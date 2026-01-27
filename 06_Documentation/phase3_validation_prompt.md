# Phase 3: Critical Validation & Multivariate Analysis
## Claude Code Comprehensive Execution Prompt

**Date:** October 11, 2025
**Objective:** Address critical statistical issues and perform rigorous external validation
**Priority:** HIGH - Publication-blocking issues must be resolved

---

## MISSION STATEMENT

You must perform **rigorous statistical validation** to address critical concerns identified in Phase 2 review:

1. **Fix Cox proportional hazards violations** using alternative models
2. **Test independent prognostic value** of clusters beyond mutations
3. **Investigate TCGA validation failure** - normalization, sample size, mutation patterns
4. **Add TARGET-AML validation** as independent external cohort
5. **Verify classifier integrity** - check for data leakage
6. **Evaluate alternative clustering solutions** (k=3, k=4, k=5)

**Scientific Standards Required:**
- All models must report assumptions checks
- Multiple testing correction applied throughout
- Honest reporting of negative results
- Transparent handling of violations

---

## PART 1: FIX COX PROPORTIONAL HAZARDS VIOLATIONS

### **Background**
Phase 2 identified violations of proportional hazards assumption:
```
Global test: p = 0.0002 (VIOLATED)
Cluster: p = 0.013 (violation)
Age: p = 0.002 (violation)
```

This invalidates all survival conclusions. Must use alternative approaches.

### **Task 1.1: Stratified Cox Regression**

**Objective:** Treat cluster as stratification variable instead of covariate

**Method:**
```r
library(survival)
library(survminer)

# Load survival data
survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")

# Stratified Cox model - do NOT include cluster as covariate
fit_stratified <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ strata(cluster_assignment) + AGE + SEX,
                        data = survival_data)

# Extract stratified results
summary(fit_stratified)

# Test if survival curves differ (log-rank test)
survdiff_result <- survdiff(Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment,
                            data = survival_data)

# Plot stratified survival
fit_km_stratified <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment,
                             data = survival_data)

pdf("04_Figures/11_Survival_Reanalysis/stratified_kaplan_meier.pdf", width=8, height=6)
ggsurvplot(fit_km_stratified,
           data = survival_data,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           legend.labs = c("Proliferative (C1)", "Immune-Inflammatory (C2)"),
           palette = c("#E41A1C", "#377EB8"),
           title = "Stratified Survival Analysis")
dev.off()
```

**Output:**
- Survival curves comparison without PH assumption
- Log-rank p-value (assumption-free test)
- File: `stratified_survival_results.csv`

---

### **Task 1.2: Time-Varying Coefficient Model**

**Objective:** Model cluster effect as changing over time

**Method:**
```r
# Time-varying coefficient using time transformation
fit_timevarying <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ 
                         cluster_assignment + tt(cluster_assignment) + AGE + SEX,
                         data = survival_data,
                         tt = function(x, t, ...) {
                           # Model cluster effect as function of log(time)
                           x * log(t + 1)  # +1 to avoid log(0)
                         })

summary(fit_timevarying)

# Test if time-varying term is significant
# If significant: effect truly changes over time
# If not: PH violation may be minor

# Create time-varying hazard ratio plot
time_points <- seq(1, max(survival_data$OS_MONTHS, na.rm=TRUE), by=1)
base_coef <- coef(fit_timevarying)["cluster_assignment"]
tv_coef <- coef(fit_timevarying)["tt(cluster_assignment)"]

hr_over_time <- exp(base_coef + tv_coef * log(time_points + 1))

pdf("04_Figures/11_Survival_Reanalysis/hazard_ratio_over_time.pdf", width=8, height=6)
plot(time_points, hr_over_time, type="l", lwd=2,
     xlab="Time (months)", ylab="Hazard Ratio (C2 vs C1)",
     main="Time-Varying Hazard Ratio",
     ylim=c(0.5, 2.5))
abline(h=1, lty=2, col="red")
dev.off()
```

**Output:**
- Time-varying HR estimates
- Plot showing how cluster effect changes over time
- File: `time_varying_cox_results.csv`

---

### **Task 1.3: Landmark Analysis**

**Objective:** Test survival differences from specific timepoints

**Method:**
```r
# Landmark analysis at 6, 12, 18, 24 months
landmarks <- c(6, 12, 18, 24)

landmark_results <- lapply(landmarks, function(lm) {
  # Subset to patients alive at landmark
  landmark_data <- survival_data %>%
    filter(OS_MONTHS > lm | (OS_MONTHS <= lm & OS_STATUS == 0))
  
  if(nrow(landmark_data) < 50) return(NULL)  # Need sufficient samples
  
  # Conditional survival from landmark
  landmark_data$landmark_time <- landmark_data$OS_MONTHS - lm
  landmark_data$landmark_time[landmark_data$landmark_time < 0] <- 0
  
  # Cox model from landmark
  fit <- coxph(Surv(landmark_time, OS_STATUS) ~ cluster_assignment + AGE + SEX,
               data = landmark_data)
  
  # Log-rank test
  lr_test <- survdiff(Surv(landmark_time, OS_STATUS) ~ cluster_assignment,
                      data = landmark_data)
  
  data.frame(
    landmark_month = lm,
    n_samples = nrow(landmark_data),
    HR = exp(coef(fit)["cluster_assignment"]),
    HR_lower = exp(confint(fit)["cluster_assignment", 1]),
    HR_upper = exp(confint(fit)["cluster_assignment", 2]),
    cox_pval = summary(fit)$coefficients["cluster_assignment", "Pr(>|z|)"],
    logrank_pval = 1 - pchisq(lr_test$chisq, df=1)
  )
}) %>% bind_rows()

write.csv(landmark_results, "03_Results/11_Survival_Reanalysis/landmark_analysis_results.csv")

# Plot landmark HRs
pdf("04_Figures/11_Survival_Reanalysis/landmark_hazard_ratios.pdf", width=8, height=6)
ggplot(landmark_results, aes(x=landmark_month, y=HR)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=HR_lower, ymax=HR_upper), width=1) +
  geom_hline(yintercept=1, linetype="dashed", color="red") +
  labs(x="Landmark Time (months)", y="Hazard Ratio (C2 vs C1)",
       title="Landmark Analysis: Conditional Survival from Different Timepoints") +
  theme_bw()
dev.off()
```

**Output:**
- Conditional survival at multiple landmarks
- Shows if cluster effect is early vs late
- File: `landmark_analysis_results.csv`

---

### **Task 1.4: Restricted Mean Survival Time (RMST)**

**Objective:** Non-parametric alternative to Cox regression (no PH assumption needed)

**Method:**
```r
library(survRM2)

# RMST at 24 months (2 years) - clinically meaningful timepoint
# Remove samples with missing cluster
complete_data <- survival_data %>% filter(!is.na(cluster_assignment))

# Convert cluster to 0/1
complete_data$arm <- ifelse(complete_data$cluster_assignment == "Cluster_1", 0, 1)

# RMST analysis
rmst_result <- rmst2(time = complete_data$OS_MONTHS,
                     status = complete_data$OS_STATUS,
                     arm = complete_data$arm,
                     tau = 24)  # Restriction time at 24 months

# Extract results
rmst_summary <- data.frame(
  Group = c("Cluster 1", "Cluster 2", "Difference"),
  RMST = c(rmst_result$RMST.arm0$rmst[1],
           rmst_result$RMST.arm1$rmst[1],
           rmst_result$unadjusted.result[1,1]),
  SE = c(rmst_result$RMST.arm0$rmst[2],
         rmst_result$RMST.arm1$rmst[2],
         rmst_result$unadjusted.result[1,2]),
  Lower_CI = c(rmst_result$RMST.arm0$rmst[3],
               rmst_result$RMST.arm1$rmst[3],
               rmst_result$unadjusted.result[1,3]),
  Upper_CI = c(rmst_result$RMST.arm0$rmst[4],
               rmst_result$RMST.arm1$rmst[4],
               rmst_result$unadjusted.result[1,4]),
  P_value = c(NA, NA, rmst_result$unadjusted.result[1,5])
)

write.csv(rmst_summary, "03_Results/11_Survival_Reanalysis/rmst_analysis_24months.csv")

# Plot RMST
pdf("04_Figures/11_Survival_Reanalysis/rmst_comparison.pdf", width=8, height=6)
plot(rmst_result, xlab="Months", ylab="Probability",
     col=c("#E41A1C", "#377EB8"),
     main="Restricted Mean Survival Time Analysis (24 months)")
dev.off()

# Also test at 36 months if sufficient follow-up
if(max(complete_data$OS_MONTHS, na.rm=TRUE) >= 36) {
  rmst_36 <- rmst2(time = complete_data$OS_MONTHS,
                   status = complete_data$OS_STATUS,
                   arm = complete_data$arm,
                   tau = 36)
  
  rmst_36_summary <- data.frame(
    Group = c("Cluster 1", "Cluster 2", "Difference"),
    RMST = c(rmst_36$RMST.arm0$rmst[1],
             rmst_36$RMST.arm1$rmst[1],
             rmst_36$unadjusted.result[1,1]),
    P_value = c(NA, NA, rmst_36$unadjusted.result[1,5])
  )
  write.csv(rmst_36_summary, "03_Results/11_Survival_Reanalysis/rmst_analysis_36months.csv")
}
```

**Output:**
- RMST at 24 months (clinically relevant timepoint)
- Direct comparison of survival time without PH assumption
- File: `rmst_analysis_24months.csv`

---

## PART 2: MULTIVARIATE ANALYSIS - TEST INDEPENDENT PROGNOSTIC VALUE

### **Critical Question**
Does cluster provide independent prognostic value beyond mutations?

### **Task 2.1: Load and Merge Data**

**Method:**
```r
# Load all necessary data
survival_data <- read.csv("03_Results/08_Survival_Analysis/survival_data_with_clusters.csv")
mutation_matrix <- readRDS("03_Results/10_Mutations/mutation_matrix.rds")
clinical_data <- readRDS("01_Data/BeatAML/cohort1_clinical.rds")

# Merge datasets
analysis_data <- survival_data %>%
  left_join(mutation_matrix, by="lab_id") %>%
  left_join(clinical_data %>% select(lab_id, ELN_RISK = eln2017_risk), by="lab_id")

# Create complete case dataset
complete_data <- analysis_data %>%
  filter(!is.na(OS_MONTHS), !is.na(OS_STATUS), !is.na(cluster_assignment))

# Save merged dataset
write.csv(complete_data, "03_Results/12_Multivariate_Analysis/merged_analysis_data.csv")

cat("Complete cases for multivariate analysis:", nrow(complete_data), "\n")
```

---

### **Task 2.2: Sequential Model Building**

**Objective:** Test if cluster adds value beyond mutations and clinical factors

**Method:**
```r
# Define models in order of complexity

# Model 1: Cluster only (baseline)
model1 <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ cluster_assignment,
                data = complete_data)

# Model 2: Key mutations only
model2 <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ NPM1 + TP53 + RUNX1 + DNMT3A + FLT3,
                data = complete_data)

# Model 3: Clinical factors only
model3 <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ AGE + SEX,
                data = complete_data)

# Model 4: Cluster + Mutations (KEY MODEL)
model4 <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ 
                cluster_assignment + NPM1 + TP53 + RUNX1 + DNMT3A + FLT3,
                data = complete_data)

# Model 5: Cluster + Clinical
model5 <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ 
                cluster_assignment + AGE + SEX,
                data = complete_data)

# Model 6: Mutations + Clinical
model6 <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ 
                NPM1 + TP53 + RUNX1 + DNMT3A + FLT3 + AGE + SEX,
                data = complete_data)

# Model 7: FULL MODEL - Cluster + Mutations + Clinical
model7 <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ 
                cluster_assignment + NPM1 + TP53 + RUNX1 + DNMT3A + FLT3 + AGE + SEX,
                data = complete_data)

# Model 8: ELN risk only (if available)
if("ELN_RISK" %in% colnames(complete_data) && sum(!is.na(complete_data$ELN_RISK)) > 100) {
  model8 <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ ELN_RISK,
                  data = complete_data %>% filter(!is.na(ELN_RISK)))
  
  # Model 9: ELN + Cluster (does cluster add value beyond ELN?)
  model9 <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ ELN_RISK + cluster_assignment,
                  data = complete_data %>% filter(!is.na(ELN_RISK)))
}

# Extract model performance
extract_model_perf <- function(model, name) {
  s <- summary(model)
  data.frame(
    Model = name,
    N = model$n,
    Events = model$nevent,
    Concordance = s$concordance[1],
    Concordance_SE = s$concordance[2],
    AIC = AIC(model),
    BIC = BIC(model),
    LogLik = model$loglik[2],
    LR_test = s$logtest[1],
    LR_pval = s$logtest[3]
  )
}

model_comparison <- rbind(
  extract_model_perf(model1, "1. Cluster only"),
  extract_model_perf(model2, "2. Mutations only"),
  extract_model_perf(model3, "3. Clinical only"),
  extract_model_perf(model4, "4. Cluster + Mutations"),
  extract_model_perf(model5, "5. Cluster + Clinical"),
  extract_model_perf(model6, "6. Mutations + Clinical"),
  extract_model_perf(model7, "7. FULL MODEL")
)

if(exists("model8")) {
  model_comparison <- rbind(
    model_comparison,
    extract_model_perf(model8, "8. ELN only"),
    extract_model_perf(model9, "9. ELN + Cluster")
  )
}

write.csv(model_comparison, "03_Results/12_Multivariate_Analysis/model_comparison.csv")
```

---

### **Task 2.3: Likelihood Ratio Tests - Nested Model Comparisons**

**Method:**
```r
# Test if cluster adds value beyond mutations
lrt_cluster_vs_mutations <- anova(model2, model4)

# Test if cluster adds value beyond clinical
lrt_cluster_vs_clinical <- anova(model3, model5)

# Test if cluster adds value to full model
lrt_full <- anova(model6, model7)

# Compile LRT results
lrt_results <- data.frame(
  Comparison = c(
    "Mutations only vs Mutations + Cluster",
    "Clinical only vs Clinical + Cluster",
    "Mutations + Clinical vs FULL (+ Cluster)"
  ),
  Model1 = c("Mutations only", "Clinical only", "Mutations + Clinical"),
  Model2 = c("Mutations + Cluster", "Clinical + Cluster", "FULL"),
  Chi_square = c(
    lrt_cluster_vs_mutations$Chisq[2],
    lrt_cluster_vs_clinical$Chisq[2],
    lrt_full$Chisq[2]
  ),
  DF = c(
    lrt_cluster_vs_mutations$Df[2],
    lrt_cluster_vs_clinical$Df[2],
    lrt_full$Chisq[2]
  ),
  P_value = c(
    lrt_cluster_vs_mutations$`Pr(>|Chi|)`[2],
    lrt_cluster_vs_clinical$`Pr(>|Chi|)`[2],
    lrt_full$`Pr(>|Chi|)`[2]
  )
)

write.csv(lrt_results, "03_Results/12_Multivariate_Analysis/likelihood_ratio_tests.csv")

# CRITICAL INTERPRETATION
if(lrt_full$`Pr(>|Chi|)`[2] < 0.05) {
  cat("\n*** CLUSTER PROVIDES INDEPENDENT PROGNOSTIC VALUE ***\n")
  cat("Cluster is significant beyond mutations and clinical factors\n")
  cat("p =", lrt_full$`Pr(>|Chi|)`[2], "\n\n")
} else {
  cat("\n*** CLUSTER DOES NOT ADD INDEPENDENT VALUE ***\n")
  cat("Cluster is redundant with mutation information\n")
  cat("p =", lrt_full$`Pr(>|Chi|)`[2], "\n\n")
}
```

---

### **Task 2.4: Extract Coefficients from Full Model**

**Method:**
```r
# Get full model results
full_results <- as.data.frame(summary(model7)$coefficients)
full_results$Variable <- rownames(full_results)
full_results <- full_results %>%
  mutate(
    HR = exp(coef),
    HR_lower = exp(coef - 1.96*`se(coef)`),
    HR_upper = exp(coef + 1.96*`se(coef)`),
    Significant = `Pr(>|z|)` < 0.05
  ) %>%
  select(Variable, coef, HR, HR_lower, HR_upper, `Pr(>|z|)`, Significant) %>%
  arrange(`Pr(>|z|)`)

write.csv(full_results, "03_Results/12_Multivariate_Analysis/full_model_coefficients.csv")

# Forest plot of full model
pdf("04_Figures/12_Multivariate_Analysis/forest_plot_full_model.pdf", width=10, height=6)
forest_data <- full_results %>%
  mutate(Variable = factor(Variable, levels = rev(Variable)))

ggplot(forest_data, aes(x=HR, y=Variable)) +
  geom_point(size=3) +
  geom_errorbarh(aes(xmin=HR_lower, xmax=HR_upper), height=0.2) +
  geom_vline(xintercept=1, linetype="dashed", color="red") +
  scale_x_log10() +
  labs(x="Hazard Ratio (log scale)", y="",
       title="Multivariate Cox Regression - Full Model") +
  theme_bw()
dev.off()
```

---

### **Task 2.5: NPM1-DNMT3A Co-mutation Analysis**

**Objective:** Test if co-mutations drive cluster membership

**Method:**
```r
# Create co-mutation categories
complete_data <- complete_data %>%
  mutate(
    NPM1_DNMT3A_status = case_when(
      NPM1 == 1 & DNMT3A == 1 ~ "Both mutated",
      NPM1 == 1 & DNMT3A == 0 ~ "NPM1 only",
      NPM1 == 0 & DNMT3A == 1 ~ "DNMT3A only",
      TRUE ~ "Neither"
    ),
    # Triple mutant (common in AML)
    Triple_mutant = ifelse(NPM1 == 1 & DNMT3A == 1 & FLT3 == 1, "Triple", "Other")
  )

# Test enrichment of co-mutations by cluster
comut_by_cluster <- table(complete_data$NPM1_DNMT3A_status, complete_data$cluster_assignment)
fisher_comut <- fisher.test(comut_by_cluster, simulate.p.value=TRUE)

# Save results
comut_results <- as.data.frame.matrix(comut_by_cluster)
comut_results$NPM1_DNMT3A_status <- rownames(comut_results)
write.csv(comut_results, "03_Results/12_Multivariate_Analysis/npm1_dnmt3a_comutation.csv")

# Survival by co-mutation status
fit_comut <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ NPM1_DNMT3A_status,
                     data = complete_data)

pdf("04_Figures/12_Multivariate_Analysis/survival_by_comutation.pdf", width=10, height=6)
ggsurvplot(fit_comut,
           data = complete_data,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           title = "Survival by NPM1-DNMT3A Co-mutation Status")
dev.off()

# Cox model with co-mutation
model_comut <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ NPM1_DNMT3A_status + cluster_assignment,
                     data = complete_data)
summary(model_comut)
```

---

## PART 3: TCGA VALIDATION - DETAILED INVESTIGATION

### **Task 3.1: Recount TCGA Samples and Check Normalization**

**Objective:** Verify sample size and compare normalization methods

**Method:**
```r
library(TCGAbiolinks)
library(SummarizedExperiment)

# Query ALL available TCGA-LAML samples
query_laml <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Check total available samples
cat("Total TCGA-LAML RNA-seq samples available:", nrow(query_laml$results[[1]]), "\n")

# Download if not already done
# GDCdownload(query_laml)

# Prepare data with multiple normalization methods
tcga_se <- GDCprepare(query_laml, summarizedExperiment = TRUE)

# Extract counts
tcga_counts <- assay(tcga_se, "unstranded")
tcga_clinical <- colData(tcga_se)

# Method 1: DESeq2 normalization (what BeatAML used)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = tcga_counts,
                               colData = tcga_clinical,
                               design = ~ 1)
dds <- estimateSizeFactors(dds)
tcga_norm_deseq2 <- counts(dds, normalized=TRUE)

# Method 2: TMM normalization (alternative)
library(edgeR)
tcga_dge <- DGEList(counts=tcga_counts)
tcga_dge <- calcNormFactors(tcga_dge, method="TMM")
tcga_norm_tmm <- cpm(tcga_dge, log=TRUE, prior.count=1)

# Method 3: TPM (if available)
# Check what BeatAML actually used
beatAML_expr <- readRDS("01_Data/BeatAML/expression_batch_corrected.rds")

# Compare data distributions
comparison_summary <- data.frame(
  Method = c("BeatAML", "TCGA DESeq2", "TCGA TMM", "TCGA raw"),
  Min = c(min(beatAML_expr, na.rm=TRUE),
          min(log2(tcga_norm_deseq2 + 1), na.rm=TRUE),
          min(tcga_norm_tmm, na.rm=TRUE),
          min(log2(tcga_counts + 1), na.rm=TRUE)),
  Max = c(max(beatAML_expr, na.rm=TRUE),
          max(log2(tcga_norm_deseq2 + 1), na.rm=TRUE),
          max(tcga_norm_tmm, na.rm=TRUE),
          max(log2(tcga_counts + 1), na.rm=TRUE)),
  Mean = c(mean(beatAML_expr, na.rm=TRUE),
           mean(log2(tcga_norm_deseq2 + 1), na.rm=TRUE),
           mean(tcga_norm_tmm, na.rm=TRUE),
           mean(log2(tcga_counts + 1), na.rm=TRUE)),
  SD = c(sd(beatAML_expr, na.rm=TRUE),
         sd(log2(tcga_norm_deseq2 + 1), na.rm=TRUE),
         sd(tcga_norm_tmm, na.rm=TRUE),
         sd(log2(tcga_counts + 1), na.rm=TRUE))
)

write.csv(comparison_summary, "03_Results/17_TCGA_Validation/normalization_comparison.csv")

cat("\n=== NORMALIZATION COMPARISON ===\n")
print(comparison_summary)

# Check if BeatAML used batch correction
if(min(beatAML_expr, na.rm=TRUE) < 0) {
  cat("\nBeatAML has negative values - batch corrected with ComBat\n")
  cat("TCGA needs same treatment!\n")
}
```

---

### **Task 3.2: Apply Same Normalization as BeatAML**

**Method:**
```r
library(sva)

# Load BeatAML processing script to see exact normalization
# Replicate the EXACT same process

# Step 1: DESeq2 normalization
tcga_log2 <- log2(tcga_norm_deseq2 + 1)

# Step 2: Check if batch correction needed
# BeatAML used ComBat on cohort1 vs cohort2
# TCGA is single batch, so just normalize

# Step 3: Gene-wise Z-score to match BeatAML distribution
common_genes <- intersect(rownames(tcga_log2), rownames(beatAML_expr))
cat("Common genes between TCGA and BeatAML:", length(common_genes), "\n")

tcga_matched <- tcga_log2[common_genes, ]
beatAML_matched <- beatAML_expr[common_genes, ]

# For each gene, match TCGA to BeatAML mean and SD
tcga_normalized <- tcga_matched
for(gene in common_genes) {
  # Z-score TCGA
  tcga_z <- (tcga_matched[gene, ] - mean(tcga_matched[gene, ])) / sd(tcga_matched[gene, ])
  
  # Scale to BeatAML distribution
  tcga_normalized[gene, ] <- tcga_z * sd(beatAML_matched[gene, ]) + mean(beatAML_matched[gene, ])
}

# Save properly normalized TCGA data
saveRDS(tcga_normalized, "01_Data/TCGA_LAML/tcga_normalized_matched_to_beatAML.rds")

# Verify matching
verification <- data.frame(
  Gene = common_genes[1:10],
  BeatAML_mean = colMeans(beatAML_matched[1:10, ]),
  TCGA_mean = colMeans(tcga_normalized[1:10, ]),
  BeatAML_sd = apply(beatAML_matched[1:10, ], 1, sd),
  TCGA_sd = apply(tcga_normalized[1:10, ], 1, sd)
)

print(verification)
```

---

### **Task 3.3: Re-apply Classifier with Proper Normalization**

**Method:**
```r
# Load classifier
rf_classifier <- readRDS("03_Results/15_Gene_Signature/final_rf_classifier.rds")
gene_signature <- read.csv("03_Results/15_Gene_Signature/50_gene_signature.csv")

# Prepare TCGA data for prediction
tcga_for_prediction <- as.data.frame(t(tcga_normalized))
tcga_for_prediction <- tcga_for_prediction[, gene_signature$Gene[gene_signature$Gene %in% colnames(tcga_for_prediction)]]

# Handle missing genes
missing_genes <- setdiff(gene_signature$Gene, colnames(tcga_for_prediction))
if(length(missing_genes) > 0) {
  cat("Missing genes in TCGA:", length(missing_genes), "\n")
  cat("Missing:", paste(missing_genes, collapse=", "), "\n")
  
  # Add missing genes as zeros (or median imputation)
  for(gene in missing_genes) {
    tcga_for_prediction[[gene]] <- 0  # Conservative imputation
  }
}

# Predict
tcga_predictions <- predict(rf_classifier, newdata = tcga_for_prediction, type="prob")
tcga_classes <- predict(rf_classifier, newdata = tcga_for_prediction)

# Get sample IDs
tcga_sample_ids <- colnames(tcga_normalized)

tcga_results <- data.frame(
  Sample = tcga_sample_ids,
  Predicted_Cluster = tcga_classes,
  Prob_C1 = tcga_predictions[, "Cluster_1"],
  Prob_C2 = tcga_predictions[, "Cluster_2"],
  Confidence = apply(tcga_predictions, 1, max)
)

# Compare to previous TCGA results
prev_results <- read.csv("03_Results/17_TCGA_Validation/tcga_sample_predictions.csv")

comparison <- data.frame(
  Previous_n = nrow(prev_results),
  Current_n = nrow(tcga_results),
  Previous_C1_pct = mean(prev_results$Predicted_Cluster == "Cluster_1") * 100,
  Current_C1_pct = mean(tcga_results$Predicted_Cluster == "Cluster_1") * 100,
  Previous_confidence = mean(prev_results$Confidence),
  Current_confidence = mean(tcga_results$Confidence)
)

write.csv(comparison, "03_Results/17_TCGA_Validation/before_after_normalization.csv")
write.csv(tcga_results, "03_Results/17_TCGA_Validation/tcga_predictions_corrected_norm.csv")
```

---

### **Task 3.4: TCGA Survival Analysis with All Available Samples**

**Method:**
```r
# Merge with clinical data
tcga_clinical_data <- as.data.frame(colData(tcga_se)) %>%
  select(barcode, vital_status, days_to_death, days_to_last_follow_up, 
         age_at_diagnosis, gender, any_other_relevant_vars)

# Create survival object
tcga_surv_data <- tcga_results %>%
  left_join(tcga_clinical_data, by=c("Sample"="barcode")) %>%
  mutate(
    OS_MONTHS = ifelse(!is.na(days_to_death),
                       days_to_death / 30.44,
                       days_to_last_follow_up / 30.44),
    OS_STATUS = ifelse(vital_status == "Dead", 1, 0)
  ) %>%
  filter(!is.na(OS_MONTHS), !is.na(OS_STATUS))

cat("TCGA samples with survival data:", nrow(tcga_surv_data), "\n")
cat("TCGA death events:", sum(tcga_surv_data$OS_STATUS), "\n")

# Survival analysis
fit_tcga <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ Predicted_Cluster,
                    data = tcga_surv_data)

survdiff_tcga <- survdiff(Surv(OS_MONTHS, OS_STATUS) ~ Predicted_Cluster,
                          data = tcga_surv_data)

cox_tcga <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ Predicted_Cluster,
                  data = tcga_surv_data)

# Save results
tcga_survival_summary <- data.frame(
  N_total = nrow(tcga_surv_data),
  N_events = sum(tcga_surv_data$OS_STATUS),
  N_C1 = sum(tcga_surv_data$Predicted_Cluster == "Cluster_1"),
  N_C2 = sum(tcga_surv_data$Predicted_Cluster == "Cluster_2"),
  Median_survival_C1 = summary(fit_tcga)$table["Predicted_Cluster=Cluster_1", "median"],
  Median_survival_C2 = summary(fit_tcga)$table["Predicted_Cluster=Cluster_2", "median"],
  HR = exp(coef(cox_tcga)),
  HR_lower = exp(confint(cox_tcga)[1]),
  HR_upper = exp(confint(cox_tcga)[2]),
  Cox_pval = summary(cox_tcga)$coefficients[,"Pr(>|z|)"],
  LogRank_pval = 1 - pchisq(survdiff_tcga$chisq, df=1)
)

write.csv(tcga_survival_summary, "03_Results/17_TCGA_Validation/tcga_survival_corrected.csv")

# Plot
pdf("04_Figures/16_TCGA_Validation/tcga_km_corrected_normalization.pdf", width=10, height=8)
ggsurvplot(fit_tcga,
           data = tcga_surv_data,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           legend.labs = c("Proliferative (C1)", "Immune-Inflammatory (C2)"),
           title = "TCGA-LAML External Validation (Corrected Normalization)")
dev.off()
```

---

### **Task 3.5: TCGA Mutation Validation**

**Objective:** Check if TCGA clusters show same mutation enrichment

**Method:**
```r
# Query TCGA mutations
query_mut <- GDCquery(
  project = "TCGA-LAML",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

# GDCdownload(query_mut)
tcga_maf <- GDCprepare(query_mut)

# Create mutation matrix for key genes
key_genes <- c("NPM1", "DNMT3A", "TP53", "FLT3", "RUNX1", "ASXL1", "IDH1", "IDH2")

tcga_mutation_matrix <- tcga_maf %>%
  filter(Hugo_Symbol %in% key_genes,
         Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation", 
                                       "Frame_Shift_Del", "Frame_Shift_Ins",
                                       "Splice_Site")) %>%
  select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  distinct() %>%
  mutate(Mutated = 1) %>%
  spread(Hugo_Symbol, Mutated, fill = 0)

# Merge with cluster assignments
tcga_mut_cluster <- tcga_results %>%
  left_join(tcga_mutation_matrix, by=c("Sample"="Tumor_Sample_Barcode"))

# Test enrichment
enrichment_tcga <- lapply(key_genes, function(gene) {
  if(!gene %in% colnames(tcga_mut_cluster)) return(NULL)
  
  contingency <- table(tcga_mut_cluster$Predicted_Cluster, tcga_mut_cluster[[gene]])
  fisher_result <- fisher.test(contingency)
  
  data.frame(
    Gene = gene,
    C1_freq = mean(tcga_mut_cluster[tcga_mut_cluster$Predicted_Cluster=="Cluster_1", gene], na.rm=TRUE),
    C2_freq = mean(tcga_mut_cluster[tcga_mut_cluster$Predicted_Cluster=="Cluster_2", gene], na.rm=TRUE),
    OR = fisher_result$estimate,
    P_value = fisher_result$p.value
  )
}) %>% bind_rows()

enrichment_tcga$FDR <- p.adjust(enrichment_tcga$P_value, method="BH")

write.csv(enrichment_tcga, "03_Results/17_TCGA_Validation/tcga_mutation_enrichment.csv")

# Compare to BeatAML
beatAML_enrichment <- read.csv("03_Results/10_Mutations/mutation_enrichment_by_cluster.csv")

comparison_mut <- enrichment_tcga %>%
  left_join(beatAML_enrichment %>% select(Gene, BeatAML_OR=Odds_Ratio, BeatAML_pval=P_Value),
            by="Gene") %>%
  mutate(
    Direction_matches = sign(log(OR)) == sign(log(BeatAML_OR)),
    Both_significant = P_value < 0.05 & BeatAML_pval < 0.05
  )

write.csv(comparison_mut, "03_Results/17_TCGA_Validation/beatAML_vs_TCGA_mutations.csv")

cat("\n=== MUTATION VALIDATION IN TCGA ===\n")
print(comparison_mut)
```

---

### **Task 3.6: Power Calculation for TCGA**

**Method:**
```r
library(powerSurvEpi)

# What HR can we detect with TCGA sample size?
tcga_n <- nrow(tcga_surv_data)
tcga_events <- sum(tcga_surv_data$OS_STATUS)
tcga_c1_n <- sum(tcga_surv_data$Predicted_Cluster == "Cluster_1")
tcga_c2_n <- sum(tcga_surv_data$Predicted_Cluster == "Cluster_2")

# Test range of HRs
hr_range <- seq(1.2, 3.0, by=0.1)
power_results <- data.frame(
  HR = hr_range,
  Power = sapply(hr_range, function(hr) {
    powerCT(nE = tcga_c2_n, nC = tcga_c1_n,
            pE = 0.6, pC = 0.6,  # approximate event rates
            RR = hr,
            alpha = 0.05)
  })
)

# Find minimum detectable HR with 80% power
min_detectable_hr <- power_results$HR[which(power_results$Power >= 0.80)[1]]

power_summary <- data.frame(
  TCGA_total_n = tcga_n,
  TCGA_events = tcga_events,
  BeatAML_observed_HR = 1.22,  # From Phase 2
  TCGA_observed_HR = exp(coef(cox_tcga)),
  Min_detectable_HR_80pct = min_detectable_hr,
  Power_for_observed_HR = power_results$Power[which.min(abs(power_results$HR - 1.22))]
)

write.csv(power_summary, "03_Results/17_TCGA_Validation/tcga_power_calculation.csv")

cat("\n=== TCGA POWER ANALYSIS ===\n")
cat("Minimum detectable HR with 80% power:", min_detectable_hr, "\n")
cat("Power to detect BeatAML HR (1.22):", 
    power_summary$Power_for_observed_HR * 100, "%\n")

if(power_summary$Power_for_observed_HR < 0.5) {
  cat("\nTCGA IS UNDERPOWERED to detect the effect size observed in BeatAML\n")
  cat("This explains the non-significant result\n")
}
```

---

## PART 4: TARGET-AML EXTERNAL VALIDATION

### **Task 4.1: Download TARGET-AML Data**

**Method:**
```r
# Query TARGET-AML project
query_target <- GDCquery(
  project = "TARGET-AML",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

cat("TARGET-AML samples available:", nrow(query_target$results[[1]]), "\n")

# Download
GDCdownload(query_target)

# Prepare
target_se <- GDCprepare(query_target, summarizedExperiment = TRUE)

# Extract data
target_counts <- assay(target_se, "unstranded")
target_clinical <- colData(target_se)

# Normalize using DESeq2 (match BeatAML)
library(DESeq2)
dds_target <- DESeqDataSetFromMatrix(countData = target_counts,
                                     colData = target_clinical,
                                     design = ~ 1)
dds_target <- estimateSizeFactors(dds_target)
target_norm <- counts(dds_target, normalized=TRUE)
target_log2 <- log2(target_norm + 1)

# Match to BeatAML distribution (same as TCGA)
common_genes_target <- intersect(rownames(target_log2), rownames(beatAML_expr))
cat("Common genes with BeatAML:", length(common_genes_target), "\n")

target_normalized <- target_log2[common_genes_target, ]
for(gene in common_genes_target) {
  target_z <- (target_normalized[gene, ] - mean(target_normalized[gene, ])) / sd(target_normalized[gene, ])
  target_normalized[gene, ] <- target_z * sd(beatAML_expr[gene, ]) + mean(beatAML_expr[gene, ])
}

saveRDS(target_normalized, "01_Data/TARGET_AML/target_normalized_matched.rds")
```

---

### **Task 4.2: Apply Classifier to TARGET-AML**

**Method:**
```r
# Prepare for prediction
target_for_pred <- as.data.frame(t(target_normalized))
target_for_pred <- target_for_pred[, gene_signature$Gene[gene_signature$Gene %in% colnames(target_for_pred)]]

# Handle missing genes
missing_target <- setdiff(gene_signature$Gene, colnames(target_for_pred))
if(length(missing_target) > 0) {
  for(gene in missing_target) {
    target_for_pred[[gene]] <- 0
  }
}

# Predict
target_predictions <- predict(rf_classifier, newdata = target_for_pred, type="prob")
target_classes <- predict(rf_classifier, newdata = target_for_pred)

target_results <- data.frame(
  Sample = colnames(target_normalized),
  Predicted_Cluster = target_classes,
  Prob_C1 = target_predictions[, "Cluster_1"],
  Prob_C2 = target_predictions[, "Cluster_2"],
  Confidence = apply(target_predictions, 1, max)
)

# Cluster proportions
target_proportions <- table(target_results$Predicted_Cluster)
cat("\n=== TARGET-AML CLUSTER PROPORTIONS ===\n")
print(prop.table(target_proportions))

write.csv(target_results, "03_Results/18_TARGET_Validation/target_predictions.csv")
```

---

### **Task 4.3: TARGET-AML Survival Analysis**

**Method:**
```r
# Extract clinical data
target_clinical_df <- as.data.frame(target_clinical) %>%
  select(barcode, vital_status, days_to_death, days_to_last_follow_up,
         age_at_diagnosis, gender)

# Merge with predictions
target_surv_data <- target_results %>%
  left_join(target_clinical_df, by=c("Sample"="barcode")) %>%
  mutate(
    OS_MONTHS = ifelse(!is.na(days_to_death),
                       days_to_death / 30.44,
                       days_to_last_follow_up / 30.44),
    OS_STATUS = ifelse(vital_status == "Dead", 1, 0),
    AGE = age_at_diagnosis / 365.25  # Convert to years
  ) %>%
  filter(!is.na(OS_MONTHS), !is.na(OS_STATUS))

cat("TARGET samples with survival data:", nrow(target_surv_data), "\n")
cat("Death events:", sum(target_surv_data$OS_STATUS), "\n")

# Survival analysis
fit_target <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ Predicted_Cluster,
                      data = target_surv_data)

survdiff_target <- survdiff(Surv(OS_MONTHS, OS_STATUS) ~ Predicted_Cluster,
                            data = target_surv_data)

cox_target <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ Predicted_Cluster + AGE + gender,
                    data = target_surv_data)

# Results summary
target_survival_summary <- data.frame(
  Cohort = "TARGET-AML",
  N_total = nrow(target_surv_data),
  N_events = sum(target_surv_data$OS_STATUS),
  N_C1 = sum(target_surv_data$Predicted_Cluster == "Cluster_1"),
  N_C2 = sum(target_surv_data$Predicted_Cluster == "Cluster_2"),
  Median_C1 = summary(fit_target)$table["Predicted_Cluster=Cluster_1", "median"],
  Median_C2 = summary(fit_target)$table["Predicted_Cluster=Cluster_2", "median"],
  HR = exp(coef(cox_target)["Predicted_ClusterCluster_2"]),
  P_value = summary(cox_target)$coefficients["Predicted_ClusterCluster_2", "Pr(>|z|)"]
)

write.csv(target_survival_summary, "03_Results/18_TARGET_Validation/target_survival.csv")

# Plot
pdf("04_Figures/17_TARGET_Validation/target_kaplan_meier.pdf", width=10, height=8)
ggsurvplot(fit_target,
           data = target_surv_data,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE,
           legend.labs = c("Proliferative (C1)", "Immune-Inflammatory (C2)"),
           title = "TARGET-AML External Validation")
dev.off()
```

---

### **Task 4.4: TARGET Mutation Validation**

**Method:**
```r
# Query TARGET mutations
query_target_mut <- GDCquery(
  project = "TARGET-AML",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation"
)

GDCdownload(query_target_mut)
target_maf <- GDCprepare(query_target_mut)

# Create mutation matrix
target_mut_matrix <- target_maf %>%
  filter(Hugo_Symbol %in% key_genes,
         Variant_Classification %in% c("Missense_Mutation", "Nonsense_Mutation",
                                       "Frame_Shift_Del", "Frame_Shift_Ins",
                                       "Splice_Site")) %>%
  select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  distinct() %>%
  mutate(Mutated = 1) %>%
  spread(Hugo_Symbol, Mutated, fill = 0)

# Enrichment analysis
target_mut_cluster <- target_results %>%
  left_join(target_mut_matrix, by=c("Sample"="Tumor_Sample_Barcode"))

enrichment_target <- lapply(key_genes, function(gene) {
  if(!gene %in% colnames(target_mut_cluster)) return(NULL)
  
  contingency <- table(target_mut_cluster$Predicted_Cluster, target_mut_cluster[[gene]])
  fisher_result <- fisher.test(contingency)
  
  data.frame(
    Gene = gene,
    C1_freq = mean(target_mut_cluster[target_mut_cluster$Predicted_Cluster=="Cluster_1", gene], na.rm=TRUE),
    C2_freq = mean(target_mut_cluster[target_mut_cluster$Predicted_Cluster=="Cluster_2", gene], na.rm=TRUE),
    OR = fisher_result$estimate,
    P_value = fisher_result$p.value
  )
}) %>% bind_rows()

write.csv(enrichment_target, "03_Results/18_TARGET_Validation/target_mutation_enrichment.csv")
```

---

## PART 5: META-ANALYSIS ACROSS COHORTS

### **Task 5.1: Combine Results from All Cohorts**

**Method:**
```r
# Combine survival results
beatAML_summary <- data.frame(
  Cohort = "BeatAML",
  N = nrow(survival_data),
  N_events = sum(survival_data$OS_STATUS),
  HR = 1.22,  # From Phase 1
  P_value = 0.001  # Approximate
)

meta_analysis_data <- rbind(
  beatAML_summary,
  tcga_survival_summary %>% mutate(Cohort = "TCGA-LAML") %>% 
    select(Cohort, N=N_total, N_events, HR, P_value=Cox_pval),
  target_survival_summary %>% select(Cohort, N=N_total, N_events, HR, P_value)
)

# Meta-analysis using random effects
library(meta)
meta_result <- metagen(
  TE = log(meta_analysis_data$HR),
  seTE = sqrt(1/meta_analysis_data$N_events),  # Approximate SE
  studlab = meta_analysis_data$Cohort,
  sm = "HR",
  random = TRUE,
  method.tau = "DL"
)

# Forest plot
pdf("04_Figures/18_Meta_Analysis/forest_plot_all_cohorts.pdf", width=10, height=6)
forest(meta_result,
       leftcols = c("studlab", "n", "event"),
       leftlabs = c("Cohort", "N", "Events"),
       xlab = "Hazard Ratio (C2 vs C1)")
dev.off()

# Save meta-analysis results
meta_summary <- data.frame(
  Pooled_HR = exp(meta_result$TE.random),
  Lower_CI = exp(meta_result$lower.random),
  Upper_CI = exp(meta_result$upper.random),
  P_value = meta_result$pval.random,
  I2 = meta_result$I2,
  Heterogeneity_p = meta_result$pval.Q
)

write.csv(meta_summary, "03_Results/19_Meta_Analysis/pooled_survival_meta_analysis.csv")
```

---

## PART 6: VERIFY CLASSIFIER INTEGRITY

### **Task 6.1: Check for Data Leakage**

**Objective:** Verify feature selection was done ONLY on training data

**Method:**
```r
# Load original script
original_script <- readLines("02_Scripts/Phase2_Validation/11_develop_minimal_gene_signature.R")

# Check workflow
cat("=== CHECKING CLASSIFIER WORKFLOW ===\n")
cat("Searching for train/test split location...\n\n")

# Find where data split occurs
split_line <- grep("createDataPartition|sample|train_indices", original_script)
cat("Data split at lines:", split_line, "\n")

# Find where LASSO/RF are fit
lasso_line <- grep("cv.glmnet|glmnet", original_script)
rf_line <- grep("randomForest|train.*rf", original_script)

cat("LASSO fit at lines:", lasso_line, "\n")
cat("RF fit at lines:", rf_line, "\n")

# CRITICAL CHECK: Is split BEFORE feature selection?
if(min(split_line) < min(c(lasso_line, rf_line))) {
  cat("\nâœ… CORRECT: Train/test split BEFORE feature selection\n")
  cat("No data leakage detected\n")
} else {
  cat("\nâš ï¸ WARNING: Feature selection may have used full dataset\n")
  cat("MUST RE-RUN with proper nested cross-validation\n")
}

# If leakage detected, create corrected script
if(min(split_line) >= min(c(lasso_line, rf_line))) {
  # Generate corrected analysis
  source("02_Scripts/Phase2_Validation/11_corrected_gene_signature.R")
}
```

---

### **Task 6.2: Nested Cross-Validation (if needed)**

**Method:**
```r
# Only run if data leakage was detected
library(caret)
library(glmnet)
library(randomForest)

# Load data
expression_data <- readRDS("01_Data/BeatAML/expression_batch_corrected.rds")
cluster_assignments <- read.csv("03_Results/02_Clustering/cluster_assignments.csv")

# Merge
data_for_ml <- as.data.frame(t(expression_data)) %>%
  rownames_to_column("lab_id") %>%
  inner_join(cluster_assignments, by="lab_id")

# STEP 1: Split FIRST
set.seed(42)
train_indices <- createDataPartition(data_for_ml$cluster_assignment, p=0.7, list=FALSE)
train_data <- data_for_ml[train_indices, ]
test_data <- data_for_ml[-train_indices, ]

# STEP 2: Feature selection ONLY on training data
train_matrix <- as.matrix(train_data %>% select(-lab_id, -cluster_assignment))
train_labels <- train_data$cluster_assignment

# LASSO on training only
cv_lasso <- cv.glmnet(train_matrix, train_labels, 
                      family="binomial", alpha=1, nfolds=10)
lasso_coef <- coef(cv_lasso, s="lambda.min")
lasso_genes <- rownames(lasso_coef)[which(lasso_coef != 0)][-1]  # Remove intercept

# Random Forest on training only
rf_model <- randomForest(cluster_assignment ~ ., 
                         data = train_data %>% select(-lab_id),
                         importance = TRUE, ntree=500)
rf_importance <- importance(rf_model)
rf_genes <- rownames(rf_importance)[order(-rf_importance[,"MeanDecreaseGini"])][1:50]

# STEP 3: Test on held-out test set
final_genes <- unique(c(lasso_genes, rf_genes))[1:50]
test_matrix <- test_data %>% select(all_of(final_genes))

# Train final classifier on training data only
final_rf <- randomForest(x = train_data %>% select(all_of(final_genes)),
                         y = train_data$cluster_assignment,
                         ntree=500)

# Predict on test set
test_predictions <- predict(final_rf, newdata = test_matrix)
test_probs <- predict(final_rf, newdata = test_matrix, type="prob")

# Evaluate
cm <- confusionMatrix(test_predictions, test_data$cluster_assignment)

corrected_performance <- data.frame(
  Method = "Nested CV (Corrected)",
  Accuracy = cm$overall["Accuracy"],
  Sensitivity = cm$byClass["Sensitivity"],
  Specificity = cm$byClass["Specificity"],
  AUC = auc(roc(test_data$cluster_assignment, test_probs[, 2]))
)

write.csv(corrected_performance, "03_Results/15_Gene_Signature/corrected_classifier_performance.csv")

cat("\n=== CORRECTED CLASSIFIER PERFORMANCE ===\n")
print(corrected_performance)
```

---

## PART 7: ALTERNATIVE CLUSTERING SOLUTIONS

### **Task 7.1: Evaluate k=3, k=4, k=5**

**Method:**
```r
library(cluster)
library(factoextra)

# Load expression data
expr_scaled <- readRDS("01_Data/BeatAML/expression_batch_corrected.rds")

# Test k=2 to k=6
k_range <- 2:6

# Silhouette analysis
silhouette_scores <- sapply(k_range, function(k) {
  set.seed(42)
  km <- kmeans(t(expr_scaled), centers=k, nstart=25, iter.max=100)
  sil <- silhouette(km$cluster, dist(t(expr_scaled)))
  mean(sil[, 3])
})

# Gap statistic
gap_stat <- clusGap(t(expr_scaled), FUN=kmeans, nstart=25, K.max=6, B=50)

# Plot silhouette
pdf("04_Figures/19_Alternative_Clustering/silhouette_by_k.pdf", width=8, height=6)
plot(k_range, silhouette_scores, type="b", pch=19,
     xlab="Number of Clusters (k)", ylab="Average Silhouette Width",
     main="Optimal k Selection - Silhouette Method")
abline(v=which.max(silhouette_scores)+1, lty=2, col="red")
dev.off()

# Plot gap statistic
pdf("04_Figures/19_Alternative_Clustering/gap_statistic.pdf", width=8, height=6)
fviz_gap_stat(gap_stat)
dev.off()

# For each k, perform clustering and survival analysis
for(k in 3:5) {
  set.seed(42)
  km_k <- kmeans(t(expr_scaled), centers=k, nstart=25, iter.max=100)
  
  # Assign clusters
  cluster_k <- data.frame(
    lab_id = colnames(expr_scaled),
    cluster = paste0("Cluster_", km_k$cluster)
  )
  
  # Merge with survival
  surv_k <- survival_data %>%
    inner_join(cluster_k, by="lab_id") %>%
    filter(!is.na(OS_MONTHS), !is.na(OS_STATUS))
  
  # Survival analysis
  fit_k <- survfit(Surv(OS_MONTHS, OS_STATUS) ~ cluster, data = surv_k)
  survdiff_k <- survdiff(Surv(OS_MONTHS, OS_STATUS) ~ cluster, data = surv_k)
  
  # Plot
  pdf(paste0("04_Figures/19_Alternative_Clustering/km_k", k, "_survival.pdf"), 
      width=10, height=8)
  print(ggsurvplot(fit_k, data = surv_k, pval = TRUE, risk.table = TRUE,
                   title = paste0("Survival Analysis - k=", k)))
  dev.off()
  
  # Save cluster assignments
  write.csv(cluster_k, paste0("03_Results/20_Alternative_k/cluster_assignments_k", k, ".csv"))
  
  # Mutation enrichment for each cluster
  if(k == 3) {
    # Detailed analysis for k=3
    mut_k3 <- cluster_k %>%
      inner_join(mutation_matrix, by="lab_id")
    
    # Test NPM1 enrichment across 3 clusters
    npm1_by_k3 <- table(mut_k3$cluster, mut_k3$NPM1)
    fisher_k3 <- fisher.test(npm1_by_k3, simulate.p.value=TRUE)
    
    write.csv(as.data.frame.matrix(npm1_by_k3),
              "03_Results/20_Alternative_k/npm1_enrichment_k3.csv")
  }
}

# Summary comparison
clustering_summary <- data.frame(
  k = k_range,
  Silhouette = silhouette_scores,
  Gap_stat = gap_stat$Tab[,"gap"],
  Optimal_by_silhouette = k_range == k_range[which.max(silhouette_scores)],
  Optimal_by_gap = k_range == with(gap_stat, maxSE(Tab[,"gap"], Tab[,"SE.sim"]))
)

write.csv(clustering_summary, "03_Results/20_Alternative_k/clustering_evaluation_summary.csv")
```

---

## FINAL OUTPUTS REQUIRED

### **Directory Structure**
```
03_Results/
â"œâ"€â"€ 11_Survival_Reanalysis/
â"‚   â"œâ"€â"€ stratified_survival_results.csv
â"‚   â"œâ"€â"€ time_varying_cox_results.csv
â"‚   â"œâ"€â"€ landmark_analysis_results.csv
â"‚   â""â"€â"€ rmst_analysis_24months.csv
â"œâ"€â"€ 12_Multivariate_Analysis/
â"‚   â"œâ"€â"€ model_comparison.csv
â"‚   â"œâ"€â"€ likelihood_ratio_tests.csv
â"‚   â"œâ"€â"€ full_model_coefficients.csv
â"‚   â""â"€â"€ npm1_dnmt3a_comutation.csv
â"œâ"€â"€ 17_TCGA_Validation/
â"‚   â"œâ"€â"€ normalization_comparison.csv
â"‚   â"œâ"€â"€ before_after_normalization.csv
â"‚   â"œâ"€â"€ tcga_survival_corrected.csv
â"‚   â"œâ"€â"€ tcga_mutation_enrichment.csv
â"‚   â""â"€â"€ tcga_power_calculation.csv
â"œâ"€â"€ 18_TARGET_Validation/
â"‚   â"œâ"€â"€ target_predictions.csv
â"‚   â"œâ"€â"€ target_survival.csv
â"‚   â""â"€â"€ target_mutation_enrichment.csv
â"œâ"€â"€ 19_Meta_Analysis/
â"‚   â""â"€â"€ pooled_survival_meta_analysis.csv
â""â"€â"€ 20_Alternative_k/
    â"œâ"€â"€ clustering_evaluation_summary.csv
    â"œâ"€â"€ cluster_assignments_k3.csv
    â"œâ"€â"€ cluster_assignments_k4.csv
    â""â"€â"€ cluster_assignments_k5.csv

04_Figures/
â"œâ"€â"€ 11_Survival_Reanalysis/
â"œâ"€â"€ 12_Multivariate_Analysis/
â"œâ"€â"€ 16_TCGA_Validation/
â"œâ"€â"€ 17_TARGET_Validation/
â"œâ"€â"€ 18_Meta_Analysis/
â""â"€â"€ 19_Alternative_Clustering/
```

---

## CRITICAL SUCCESS CRITERIA

At the end of this analysis, you MUST be able to answer:

1. **Does cluster provide independent prognostic value beyond mutations?**
   - Check LRT p-value for Model 6 vs Model 7
   - If p<0.05: YES, publishable claim
   - If p>0.05: NO, need to pivot manuscript

2. **Why did TCGA validation fail?**
   - Sample size (underpowered)?
   - Normalization mismatch (fixed)?
   - Biological differences (check mutations)?
   - Treatment era (documented)?

3. **Does TARGET-AML validate the findings?**
   - If yes: Strong evidence for generalizability
   - If no + TCGA no: Major concern, need third cohort

4. **Is k=2 the optimal clustering?**
   - Check silhouette and gap statistic
   - If k=3 or k=4 better: Re-analyze with optimal k

5. **Is the classifier unbiased?**
   - Check for data leakage
   - Report honest performance estimates

---

## REPORTING REQUIREMENTS

For EACH analysis, report:
âœ… Sample sizes used
âœ… Assumptions checked
âœ… Multiple testing correction applied
âœ… Effect sizes with confidence intervals
âœ… P-values (uncorrected and FDR-corrected)
âœ… Negative results transparently

**DO NOT:**
❌ Hide non-significant results
❌ Claim causation from correlation
❌ Overstate generalizability if validations fail
❌ Use "explains" unless tested in multivariate model

---

## EXECUTION PRIORITY

**Week 1 (CRITICAL):**
1. Part 1: Cox model fixes
2. Part 2: Multivariate analysis
3. Part 3: TCGA investigation

**Week 2 (HIGH PRIORITY):**
4. Part 4: TARGET validation
5. Part 5: Meta-analysis
6. Part 6: Classifier verification

**Week 3 (IMPORTANT):**
7. Part 7: Alternative clustering

---

## BEGIN EXECUTION

Start with Part 1, Task 1.1 and proceed sequentially through all parts. Report progress after completing each major section.

**Expected runtime:** 20-30 hours of computation
**Priority:** Publication-blocking issues - complete within 2 weeks

Good luck! 🚀