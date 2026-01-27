#!/usr/bin/env python3
"""
Failed Paper Recovery Script
Analyzes failed downloads and provides alternative access strategies
"""

import os
import csv
from pathlib import Path
from typing import Dict, List
import webbrowser

# Priority ranking for paper categories (1 = highest priority)
CATEGORY_PRIORITY = {
    "Beat_AML_Consortium": 1,  # Essential - your primary data source
    "Drug_Response_Prediction": 2,  # Core to your analysis
    "AML_Molecular_Subtyping": 2,  # Core to your analysis
    "Multi-Omics_Integration": 3,  # Important for methods
    "Mutation_Expression_Studies": 3,  # Important for validation
    "Clinical_Prognostic_Models": 4,  # Useful for context
    "Computational_Methods": 4,  # Useful for methods
    "Clinical_Trials": 5,  # Background information
}

# Individual paper priorities (overrides category priority if specified)
PAPER_PRIORITY = {
    "Tyner_2018": 1,  # Original Beat AML - MUST HAVE
    "Bottomly_2022": 1,  # Latest Beat AML integration - MUST HAVE
    "Cheng_2022": 1,  # Transcriptomic subtypes - MUST HAVE
    "Trac_2023": 2,  # MDREAM model - HIGH PRIORITY
    "Zeng_2022": 2,  # Cellular hierarchy - HIGH PRIORITY
    "Pino_2024": 2,  # Proteogenomics - HIGH PRIORITY
    "TCGA_2013": 2,  # Foundational - HIGH PRIORITY
    "Severens_2024": 2,  # Multi-cohort transcriptomics - HIGH PRIORITY
    "Burd_2020": 2,  # Beat AML Master Trial - HIGH PRIORITY
    "Argelaguet_2020": 3,  # MOFA+ integration
    "Lee_2018": 3,  # MERGE method
}


class PaperRecovery:
    """Analyze failed downloads and provide recovery strategies"""
    
    def __init__(self, base_dir: str = "AML_Literature"):
        self.base_dir = Path(base_dir)
        self.failed_papers = []
        self.successful_papers = []
        
    def analyze_downloads(self) -> Dict:
        """Analyze what was downloaded and what failed"""
        
        # Get all successful downloads
        for pdf_file in self.base_dir.rglob("*.pdf"):
            self.successful_papers.append({
                'filename': pdf_file.name,
                'category': pdf_file.parent.name,
                'path': str(pdf_file)
            })
        
        # Parse failed downloads
        failed_file = self.base_dir / "failed_downloads.txt"
        if failed_file.exists():
            with open(failed_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        title = parts[0]
                        doi = parts[1].replace('DOI:', '').strip() if 'DOI:' in parts[1] else None
                        pmid = parts[2].replace('PMID:', '').strip() if 'PMID:' in parts[2] else None
                        
                        # Determine category from directory structure
                        category = "Unknown"
                        for cat in CATEGORY_PRIORITY.keys():
                            if (self.base_dir / cat).exists():
                                category = cat
                                break
                        
                        self.failed_papers.append({
                            'title': title,
                            'doi': doi,
                            'pmid': pmid,
                            'category': category,
                            'priority': self._get_priority(title, category)
                        })
        
        # Sort failed papers by priority
        self.failed_papers.sort(key=lambda x: x['priority'])
        
        return {
            'total': len(self.successful_papers) + len(self.failed_papers),
            'successful': len(self.successful_papers),
            'failed': len(self.failed_papers),
            'success_rate': len(self.successful_papers) / (len(self.successful_papers) + len(self.failed_papers)) * 100
        }
    
    def _get_priority(self, title: str, category: str) -> int:
        """Get priority for a paper"""
        # Check individual paper priority first
        for paper_name, priority in PAPER_PRIORITY.items():
            if paper_name in title:
                return priority
        
        # Fall back to category priority
        return CATEGORY_PRIORITY.get(category, 5)
    
    def generate_recovery_report(self, output_file: str = "paper_recovery_guide.md"):
        """Generate a detailed recovery guide"""
        
        stats = self.analyze_downloads()
        
        output_path = self.base_dir / output_file
        with open(output_path, 'w') as f:
            f.write("# Paper Recovery Guide\n\n")
            f.write("## Summary\n\n")
            f.write(f"- **Total papers**: {stats['total']}\n")
            f.write(f"- **Successfully downloaded**: {stats['successful']}\n")
            f.write(f"- **Failed downloads**: {stats['failed']}\n")
            f.write(f"- **Success rate**: {stats['success_rate']:.1f}%\n\n")
            
            f.write("---\n\n")
            f.write("## Priority Papers to Recover\n\n")
            
            # Group by priority
            by_priority = {}
            for paper in self.failed_papers:
                priority = paper['priority']
                if priority not in by_priority:
                    by_priority[priority] = []
                by_priority[priority].append(paper)
            
            priority_labels = {
                1: "🔴 CRITICAL (Must Have)",
                2: "🟠 HIGH PRIORITY",
                3: "🟡 MEDIUM PRIORITY",
                4: "🟢 LOW PRIORITY",
                5: "⚪ OPTIONAL"
            }
            
            for priority in sorted(by_priority.keys()):
                papers = by_priority[priority]
                f.write(f"### {priority_labels.get(priority, f'Priority {priority}')}\n\n")
                f.write(f"**{len(papers)} papers**\n\n")
                
                for paper in papers:
                    f.write(f"#### {paper['title']}\n\n")
                    
                    if paper['doi']:
                        f.write(f"- **DOI**: {paper['doi']}\n")
                        f.write(f"- **DOI Link**: https://doi.org/{paper['doi']}\n")
                    
                    if paper['pmid']:
                        f.write(f"- **PMID**: {paper['pmid']}\n")
                        f.write(f"- **PubMed**: https://pubmed.ncbi.nlm.nih.gov/{paper['pmid']}/\n")
                    
                    f.write("\n**Alternative Access Methods:**\n\n")
                    
                    # Method 1: Sci-Hub
                    if paper['doi']:
                        f.write(f"1. **Sci-Hub** (check local laws): https://sci-hub.se/{paper['doi']}\n")
                    
                    # Method 2: Google Scholar
                    search_query = paper['title'].replace('_', ' ')
                    gs_url = f"https://scholar.google.com/scholar?q={search_query.replace(' ', '+')}"
                    f.write(f"2. **Google Scholar** (may have free PDF): {gs_url}\n")
                    
                    # Method 3: ResearchGate
                    f.write(f"3. **ResearchGate**: Search for paper and request from authors\n")
                    
                    # Method 4: PubMed Central
                    if paper['pmid']:
                        pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/?term={paper['pmid']}"
                        f.write(f"4. **PMC Search**: {pmc_url}\n")
                    
                    # Method 5: Institutional access
                    f.write(f"5. **University Library**: Use institutional VPN/access\n")
                    
                    f.write("\n---\n\n")
            
            # Add manual download checklist
            f.write("\n## Manual Download Checklist\n\n")
            f.write("Copy this checklist to track your progress:\n\n")
            
            for priority in sorted(by_priority.keys()):
                if priority <= 2:  # Only include critical and high priority
                    papers = by_priority[priority]
                    f.write(f"\n### {priority_labels.get(priority, f'Priority {priority}')}\n\n")
                    for paper in papers:
                        f.write(f"- [ ] {paper['title']}\n")
            
            # Add batch download instructions
            f.write("\n\n## Batch Download Strategies\n\n")
            f.write("### Option 1: Enable Sci-Hub in Download Script\n\n")
            f.write("```python\n")
            f.write("# Edit paper_download_script.py line ~380\n")
            f.write("USE_SCIHUB = True  # Enable Sci-Hub fallback\n")
            f.write("```\n\n")
            f.write("Then re-run: `python paper_download_script.py`\n\n")
            
            f.write("### Option 2: University Library Access\n\n")
            f.write("If you have institutional access:\n")
            f.write("1. Connect to university VPN\n")
            f.write("2. Visit library website\n")
            f.write("3. Use DOI links provided above\n")
            f.write("4. Download PDFs manually\n\n")
            
            f.write("### Option 3: Request from Authors\n\n")
            f.write("Many authors provide PDFs upon request:\n")
            f.write("1. Find paper on Google Scholar\n")
            f.write("2. Click author profile\n")
            f.write("3. Send polite email requesting PDF\n")
            f.write("4. Most respond within 24-48 hours\n\n")
            
            f.write("### Option 4: Inter-Library Loan (ILL)\n\n")
            f.write("Free service at most universities:\n")
            f.write("1. Request through library website\n")
            f.write("2. Provide DOI/PMID\n")
            f.write("3. Usually receive within 1-3 days\n\n")
            
        print(f"\n✓ Recovery guide saved to: {output_path}")
        return str(output_path)
    
    def generate_direct_links_html(self, output_file: str = "quick_access_links.html"):
        """Generate HTML file with clickable links for quick access"""
        
        output_path = self.base_dir / output_file
        
        with open(output_path, 'w') as f:
            f.write("""<!DOCTYPE html>
<html>
<head>
    <title>AML Papers - Quick Access Links</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            max-width: 1200px;
            margin: 20px auto;
            padding: 20px;
            background: #f5f5f5;
        }
        .priority-section {
            background: white;
            padding: 20px;
            margin: 20px 0;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .priority-1 { border-left: 5px solid #dc3545; }
        .priority-2 { border-left: 5px solid #fd7e14; }
        .priority-3 { border-left: 5px solid #ffc107; }
        .priority-4 { border-left: 5px solid #28a745; }
        .priority-5 { border-left: 5px solid #6c757d; }
        .paper {
            margin: 15px 0;
            padding: 15px;
            background: #f8f9fa;
            border-radius: 4px;
        }
        .paper-title {
            font-size: 16px;
            font-weight: bold;
            color: #333;
            margin-bottom: 8px;
        }
        .links {
            display: flex;
            gap: 10px;
            flex-wrap: wrap;
        }
        .link-btn {
            padding: 6px 12px;
            background: #007bff;
            color: white;
            text-decoration: none;
            border-radius: 4px;
            font-size: 14px;
        }
        .link-btn:hover {
            background: #0056b3;
        }
        .scihub { background: #6f42c1; }
        .scihub:hover { background: #5a32a3; }
        .scholar { background: #28a745; }
        .scholar:hover { background: #218838; }
        h1 { color: #333; }
        h2 { color: #555; margin-top: 30px; }
        .stats {
            background: #e7f3ff;
            padding: 15px;
            border-radius: 8px;
            margin: 20px 0;
        }
    </style>
</head>
<body>
    <h1>🔬 AML Literature - Quick Access Links</h1>
    
    <div class="stats">
        <strong>Status:</strong> 
""")
            
            stats = self.analyze_downloads()
            f.write(f"{stats['successful']} downloaded, {stats['failed']} failed ({stats['success_rate']:.1f}% success rate)\n")
            f.write("    </div>\n\n")
            
            # Group by priority
            by_priority = {}
            for paper in self.failed_papers:
                priority = paper['priority']
                if priority not in by_priority:
                    by_priority[priority] = []
                by_priority[priority].append(paper)
            
            priority_labels = {
                1: "🔴 CRITICAL (Must Have)",
                2: "🟠 HIGH PRIORITY",
                3: "🟡 MEDIUM PRIORITY",
                4: "🟢 LOW PRIORITY",
                5: "⚪ OPTIONAL"
            }
            
            for priority in sorted(by_priority.keys()):
                papers = by_priority[priority]
                f.write(f'    <div class="priority-section priority-{priority}">\n')
                f.write(f'        <h2>{priority_labels.get(priority, f"Priority {priority}")} ({len(papers)} papers)</h2>\n\n')
                
                for paper in papers:
                    f.write('        <div class="paper">\n')
                    f.write(f'            <div class="paper-title">{paper["title"]}</div>\n')
                    f.write('            <div class="links">\n')
                    
                    if paper['doi']:
                        # DOI link
                        f.write(f'                <a href="https://doi.org/{paper["doi"]}" class="link-btn" target="_blank">Publisher</a>\n')
                        # Sci-Hub
                        f.write(f'                <a href="https://sci-hub.se/{paper["doi"]}" class="link-btn scihub" target="_blank">Sci-Hub</a>\n')
                    
                    if paper['pmid']:
                        # PubMed
                        f.write(f'                <a href="https://pubmed.ncbi.nlm.nih.gov/{paper["pmid"]}/" class="link-btn" target="_blank">PubMed</a>\n')
                        # PMC
                        f.write(f'                <a href="https://www.ncbi.nlm.nih.gov/pmc/?term={paper["pmid"]}" class="link-btn" target="_blank">PMC</a>\n')
                    
                    # Google Scholar
                    search_query = paper['title'].replace('_', ' ').replace(' ', '+')
                    f.write(f'                <a href="https://scholar.google.com/scholar?q={search_query}" class="link-btn scholar" target="_blank">Google Scholar</a>\n')
                    
                    f.write('            </div>\n')
                    f.write('        </div>\n\n')
                
                f.write('    </div>\n\n')
            
            f.write("""
</body>
</html>
""")
        
        print(f"✓ Quick access HTML saved to: {output_path}")
        print(f"  Open in browser to access all links easily")
        return str(output_path)
    
    def create_scihub_batch_script(self, output_file: str = "scihub_batch_download.py"):
        """Create a script to batch download from Sci-Hub"""
        
        output_path = self.base_dir / output_file
        
        with open(output_path, 'w') as f:
            f.write('''#!/usr/bin/env python3
"""
Sci-Hub Batch Download Script
⚠️  WARNING: Check your local laws regarding Sci-Hub usage
"""

import requests
import time
import re
from pathlib import Path

# Failed papers with DOIs
PAPERS = [
''')
            
            # Add papers with DOIs
            for paper in self.failed_papers:
                if paper['doi'] and paper['priority'] <= 2:  # Only critical and high priority
                    f.write(f'    {{"title": "{paper["title"]}", "doi": "{paper["doi"]}", "category": "{paper["category"]}"}},\n')
            
            f.write(''']

def download_from_scihub(doi, title, category, scihub_url="https://sci-hub.se"):
    """Download paper from Sci-Hub"""
    
    print(f"\\nDownloading: {title}")
    print(f"DOI: {doi}")
    
    try:
        # Request page
        url = f"{scihub_url}/{doi}"
        response = requests.get(url, timeout=30)
        
        if response.status_code == 200:
            # Find PDF link
            match = re.search(rb'(?:src|href)=["\\']([^"\\']*\\.pdf[^"\\']*)["\\'']', response.content)
            if match:
                pdf_url = match.group(1).decode('utf-8')
                if not pdf_url.startswith('http'):
                    pdf_url = scihub_url + '/' + pdf_url.lstrip('/')
                
                # Download PDF
                pdf_response = requests.get(pdf_url, timeout=30)
                if pdf_response.status_code == 200 and pdf_response.content[:4] == b'%PDF':
                    # Save
                    category_dir = Path("AML_Literature") / category
                    category_dir.mkdir(parents=True, exist_ok=True)
                    
                    safe_title = title.replace('/', '_').replace('\\\\', '_')
                    filepath = category_dir / f"{safe_title}__{doi.replace('/', '_')}.pdf"
                    
                    with open(filepath, 'wb') as f:
                        f.write(pdf_response.content)
                    
                    print(f"✓ SUCCESS: Saved to {filepath.name}")
                    return True
        
        print(f"✗ FAILED: Could not download")
        return False
        
    except Exception as e:
        print(f"✗ ERROR: {e}")
        return False

def main():
    print("=" * 70)
    print("SCI-HUB BATCH DOWNLOAD")
    print("=" * 70)
    print()
    print("⚠️  WARNING: Sci-Hub usage may violate copyright laws.")
    print("    Check your local laws before proceeding.")
    print()
    
    response = input("Continue? (yes/no): ")
    if response.lower() != 'yes':
        print("Exiting...")
        return
    
    print(f"\\nAttempting to download {len(PAPERS)} papers...\\n")
    
    success = 0
    failed = 0
    
    for paper in PAPERS:
        result = download_from_scihub(paper['doi'], paper['title'], paper['category'])
        if result:
            success += 1
        else:
            failed += 1
        
        time.sleep(3)  # Be polite to servers
    
    print(f"\\n{'=' * 70}")
    print(f"DOWNLOAD COMPLETE")
    print(f"{'=' * 70}")
    print(f"Success: {success}")
    print(f"Failed: {failed}")
    print(f"{'=' * 70}\\n")

if __name__ == "__main__":
    main()
''')
        
        # Make executable
        output_path.chmod(0o755)
        
        print(f"✓ Sci-Hub batch script saved to: {output_path}")
        print(f"  Run with: python {output_path.name}")
        print(f"  ⚠️  WARNING: Check local laws before using")
        return str(output_path)


def main():
    """Main execution"""
    print("=" * 70)
    print("PAPER RECOVERY & PRIORITIZATION TOOL")
    print("=" * 70)
    print()
    
    recovery = PaperRecovery()
    
    print("Analyzing downloads...")
    stats = recovery.analyze_downloads()
    
    print(f"\n{'=' * 70}")
    print("DOWNLOAD ANALYSIS")
    print(f"{'=' * 70}")
    print(f"Total papers:      {stats['total']}")
    print(f"Successfully downloaded: {stats['successful']} ({stats['success_rate']:.1f}%)")
    print(f"Failed downloads:  {stats['failed']}")
    print(f"{'=' * 70}\n")
    
    # Generate recovery guide
    print("Generating recovery guide...")
    recovery.generate_recovery_report()
    
    # Generate HTML quick access
    print("\nGenerating quick access links...")
    html_file = recovery.generate_direct_links_html()
    
    # Generate Sci-Hub batch script
    print("\nGenerating Sci-Hub batch script...")
    recovery.create_scihub_batch_script()
    
    print(f"\n{'=' * 70}")
    print("RECOVERY TOOLS GENERATED")
    print(f"{'=' * 70}")
    print("\n📄 Files created:")
    print("  1. paper_recovery_guide.md - Detailed recovery strategies")
    print("  2. quick_access_links.html - Clickable links to all papers")
    print("  3. scihub_batch_download.py - Automated Sci-Hub download")
    print()
    print("🎯 Next steps:")
    print("  1. Open quick_access_links.html in browser for easy access")
    print("  2. Focus on Priority 1 & 2 papers (marked 🔴 and 🟠)")
    print("  3. Try Sci-Hub batch script for critical papers")
    print("  4. Request remaining papers from authors via ResearchGate")
    print()
    print("💡 Tip: You already have 23 papers (37%). With Priority 1-2 papers,")
    print("       you'll have enough to proceed with your analysis!")
    print(f"{'=' * 70}\n")
    
    # Open HTML file in browser
    response = input("Open quick access links in browser? (y/n): ")
    if response.lower() == 'y':
        import webbrowser
        webbrowser.open('file://' + str(Path(html_file).absolute()))


if __name__ == "__main__":
    main()