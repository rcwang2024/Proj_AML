#!/usr/bin/env python3
"""
AML Literature Download Script
Downloads papers from DOI/PMID using multiple sources
"""

import requests
import time
import os
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import logging
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('paper_download.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class PaperDownloader:
    """Download academic papers from multiple sources"""
    
    def __init__(self, output_dir: str = "downloaded_papers", email: str = "your.email@example.com"):
        """
        Initialize downloader
        
        Args:
            output_dir: Directory to save papers
            email: Your email for Unpaywall API (required for polite usage)
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.email = email
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        })
        
        # Statistics
        self.stats = {
            'total': 0,
            'success': 0,
            'failed': 0,
            'already_exists': 0
        }
    
    def sanitize_filename(self, filename: str) -> str:
        """Remove invalid characters from filename"""
        # Remove invalid characters
        filename = re.sub(r'[<>:"/\\|?*]', '_', filename)
        # Limit length
        if len(filename) > 200:
            filename = filename[:200]
        return filename
    
    def download_from_unpaywall(self, doi: str) -> Optional[bytes]:
        """
        Download PDF from Unpaywall API (open access only)
        
        Args:
            doi: DOI of the paper
            
        Returns:
            PDF content as bytes or None
        """
        try:
            # Query Unpaywall API
            url = f"https://api.unpaywall.org/v2/{doi}?email={self.email}"
            response = self.session.get(url, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                
                # Check if open access PDF is available
                if data.get('is_oa') and data.get('best_oa_location'):
                    pdf_url = data['best_oa_location'].get('url_for_pdf')
                    if pdf_url:
                        logger.info(f"Found OA PDF via Unpaywall: {pdf_url}")
                        pdf_response = self.session.get(pdf_url, timeout=30)
                        if pdf_response.status_code == 200:
                            return pdf_response.content
            
            return None
            
        except Exception as e:
            logger.debug(f"Unpaywall error for {doi}: {e}")
            return None
    
    def download_from_pmc(self, pmid: str) -> Optional[bytes]:
        """
        Download PDF from PubMed Central
        
        Args:
            pmid: PubMed ID
            
        Returns:
            PDF content as bytes or None
        """
        try:
            # First, get PMC ID from PMID
            url = f"https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids={pmid}&format=json"
            response = self.session.get(url, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                records = data.get('records', [])
                
                if records and 'pmcid' in records[0]:
                    pmcid = records[0]['pmcid']
                    
                    # Try to download PDF from PMC
                    pdf_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/"
                    logger.info(f"Trying PMC PDF: {pdf_url}")
                    
                    pdf_response = self.session.get(pdf_url, timeout=30)
                    if pdf_response.status_code == 200 and pdf_response.content[:4] == b'%PDF':
                        return pdf_response.content
            
            return None
            
        except Exception as e:
            logger.debug(f"PMC error for PMID {pmid}: {e}")
            return None
    
    def download_from_doi_redirect(self, doi: str) -> Optional[bytes]:
        """
        Try to download directly from DOI redirect
        
        Args:
            doi: DOI of the paper
            
        Returns:
            PDF content as bytes or None
        """
        try:
            # Some publishers allow direct PDF download via doi.org redirect
            url = f"https://doi.org/{doi}"
            response = self.session.get(url, timeout=10, allow_redirects=True)
            
            # Check if we got redirected to a PDF
            if response.status_code == 200:
                content_type = response.headers.get('Content-Type', '')
                if 'pdf' in content_type.lower():
                    return response.content
                
                # Try adding /pdf to URL
                final_url = response.url
                if not final_url.endswith('.pdf'):
                    pdf_url = final_url.rstrip('/') + '.pdf'
                    pdf_response = self.session.get(pdf_url, timeout=30)
                    if pdf_response.status_code == 200 and pdf_response.content[:4] == b'%PDF':
                        return pdf_response.content
            
            return None
            
        except Exception as e:
            logger.debug(f"DOI redirect error for {doi}: {e}")
            return None
    
    def download_from_scihub(self, doi: str, scihub_url: str = "https://sci-hub.se") -> Optional[bytes]:
        """
        Download from Sci-Hub (use as last resort, check your local laws)
        
        Args:
            doi: DOI of the paper
            scihub_url: Sci-Hub mirror URL
            
        Returns:
            PDF content as bytes or None
        """
        try:
            url = f"{scihub_url}/{doi}"
            logger.info(f"Trying Sci-Hub: {url}")
            
            response = self.session.get(url, timeout=30)
            if response.status_code == 200:
                # Parse HTML to find PDF link
                # Sci-Hub embeds PDF or provides download link
                if b'<embed' in response.content or b'<iframe' in response.content:
                    # Extract PDF URL from embed/iframe
                    match = re.search(rb'(?:src|href)=["\']([^"\']*\.pdf[^"\']*)["\']', response.content)
                    if match:
                        pdf_url = match.group(1).decode('utf-8')
                        if not pdf_url.startswith('http'):
                            pdf_url = scihub_url + '/' + pdf_url.lstrip('/')
                        
                        pdf_response = self.session.get(pdf_url, timeout=30)
                        if pdf_response.status_code == 200 and pdf_response.content[:4] == b'%PDF':
                            return pdf_response.content
            
            return None
            
        except Exception as e:
            logger.debug(f"Sci-Hub error for {doi}: {e}")
            return None
    
    def download_paper(self, doi: str = None, pmid: str = None, 
                      title: str = "paper", use_scihub: bool = False) -> bool:
        """
        Download paper using multiple methods
        
        Args:
            doi: DOI of the paper
            pmid: PubMed ID
            title: Title or identifier for filename
            use_scihub: Whether to use Sci-Hub as fallback
            
        Returns:
            True if download successful, False otherwise
        """
        self.stats['total'] += 1
        
        # Create safe filename
        safe_title = self.sanitize_filename(title)
        if doi:
            safe_title += f"__{doi.replace('/', '_')}"
        elif pmid:
            safe_title += f"__PMID{pmid}"
        
        filepath = self.output_dir / f"{safe_title}.pdf"
        
        # Check if already exists
        if filepath.exists():
            logger.info(f"✓ Already exists: {filepath.name}")
            self.stats['already_exists'] += 1
            return True
        
        logger.info(f"\n{'='*60}")
        logger.info(f"Downloading: {title}")
        if doi:
            logger.info(f"DOI: {doi}")
        if pmid:
            logger.info(f"PMID: {pmid}")
        logger.info(f"{'='*60}")
        
        pdf_content = None
        source = None
        
        # Try methods in order of preference
        # 1. Unpaywall (open access, legal)
        if doi and not pdf_content:
            logger.info("Trying Unpaywall...")
            pdf_content = self.download_from_unpaywall(doi)
            if pdf_content:
                source = "Unpaywall"
        
        # 2. PubMed Central (open access, legal)
        if pmid and not pdf_content:
            logger.info("Trying PubMed Central...")
            pdf_content = self.download_from_pmc(pmid)
            if pdf_content:
                source = "PMC"
        
        # 3. Direct DOI redirect
        if doi and not pdf_content:
            logger.info("Trying DOI redirect...")
            pdf_content = self.download_from_doi_redirect(doi)
            if pdf_content:
                source = "DOI_redirect"
        
        # 4. Sci-Hub (optional, use with caution)
        if doi and use_scihub and not pdf_content:
            logger.warning("Trying Sci-Hub (check your local laws)...")
            pdf_content = self.download_from_scihub(doi)
            if pdf_content:
                source = "Sci-Hub"
        
        # Save if successful
        if pdf_content:
            try:
                with open(filepath, 'wb') as f:
                    f.write(pdf_content)
                logger.info(f"✓ SUCCESS: Saved to {filepath.name} (via {source})")
                self.stats['success'] += 1
                return True
            except Exception as e:
                logger.error(f"✗ Error saving file: {e}")
        else:
            logger.warning(f"✗ FAILED: Could not download {title}")
            # Save DOI/PMID to failed list
            with open(self.output_dir / "failed_downloads.txt", 'a') as f:
                f.write(f"{title}\tDOI:{doi}\tPMID:{pmid}\n")
        
        self.stats['failed'] += 1
        time.sleep(2)  # Be polite to servers
        return False
    
    def print_summary(self):
        """Print download summary"""
        logger.info(f"\n{'='*60}")
        logger.info("DOWNLOAD SUMMARY")
        logger.info(f"{'='*60}")
        logger.info(f"Total papers:      {self.stats['total']}")
        logger.info(f"Successfully downloaded: {self.stats['success']}")
        logger.info(f"Already existed:   {self.stats['already_exists']}")
        logger.info(f"Failed:            {self.stats['failed']}")
        logger.info(f"{'='*60}\n")


def parse_paper_list() -> Dict[str, List[Dict]]:
    """Parse the paper list from the document"""
    
    papers = {
        "Beat_AML_Consortium": [
            {"title": "Tyner_2018", "doi": "10.1038/s41586-018-0623-z", "pmid": "30333627"},
            {"title": "Bottomly_2022", "doi": "10.1016/j.ccell.2022.07.002", "pmid": "35868306"},
            {"title": "Burd_2020", "doi": "10.1038/s41591-020-1089-8", "pmid": "33169022"},
            {"title": "Zhang_2019", "doi": "10.1038/s41467-018-08263-x", "pmid": "30651561"},
            {"title": "Joshi_2021", "doi": "10.1016/j.ccell.2021.06.003", "pmid": "34171263"},
            {"title": "Romine_2021", "doi": "10.1158/2643-3230.BCD-21-0012", "pmid": "34568834"},
            {"title": "Yang_2022", "doi": "10.1182/blood.2021011354", "pmid": "34482403"},
            {"title": "Romine_2023", "doi": "10.3389/fonc.2023.1192829", "pmid": "37359033"},
            {"title": "Drusbosky_2019", "doi": "10.1016/j.leukres.2018.11.010", "pmid": "30642575"},
            {"title": "Lachowiez_2023", "doi": "10.1182/bloodadvances.2022009010", "pmid": "36512707"},
        ],
        
        "Multi-Omics_Integration": [
            {"title": "Argelaguet_2018", "doi": "10.15252/msb.20178124", "pmid": "29925568"},
            {"title": "Argelaguet_2020", "doi": "10.1186/s13059-020-02015-1", "pmid": "32393369"},
            {"title": "Cantini_2021", "doi": "10.1038/s41467-020-20430-7", "pmid": "33402679"},
            {"title": "Lee_2018", "doi": "10.1038/s41467-017-02465-5", "pmid": "29295995"},
            {"title": "Sanders_2024", "doi": "10.1038/s41588-024-01999-x", "pmid": "39558088"},
            {"title": "Wang_2014", "doi": "10.1038/nmeth.2810", "pmid": "24464287"},
            {"title": "Wang_2024", "doi": "10.1016/j.heliyon.2024.e37155", "pmid": "39286085"},
            {"title": "Nicora_2020", "doi": "10.3389/fonc.2020.01030", "pmid": "32850447"},
            {"title": "Ahn_2021", "doi": "10.14348/molcells.2021.0042", "pmid": "34226278"},
            {"title": "Lazopoulou_2024", "doi": "10.1016/j.compbiomed.2024.108598", "pmid": "38901255"},
        ],
        
        "AML_Molecular_Subtyping": [
            {"title": "TCGA_2013", "doi": "10.1056/NEJMoa1301689", "pmid": "23634996"},
            {"title": "Papaemmanuil_2016", "doi": "10.1056/NEJMoa1516192", "pmid": "27276561"},
            {"title": "Valk_2004", "doi": "10.1056/NEJMoa040465", "pmid": "15084694"},
            {"title": "Bullinger_2004", "doi": "10.1056/NEJMoa031046", "pmid": "15084693"},
            {"title": "Cheng_2022", "doi": "10.1073/pnas.2211429119", "pmid": "36442087"},
            {"title": "Severens_2024", "doi": "10.1038/s41375-024-02137-6", "pmid": "38360865"},
            {"title": "Mou_2021", "doi": "10.1002/ajh.26141", "pmid": "33625756"},
            {"title": "Wang_2023", "doi": "10.1038/s41408-023-00836-4", "pmid": "37193681"},
        ],
        
        "Mutation_Expression_Studies": [
            {"title": "Russler-Germain_2014", "doi": None, "pmid": "24656771"},
            {"title": "Falini_2020", "doi": "10.1182/blood.2020006324", "pmid": "32810237"},
            {"title": "Figueroa_2010", "doi": "10.1016/j.ccr.2010.11.015", "pmid": "21130701"},
            {"title": "Palau_2023", "doi": "10.1038/s41375-023-01972-3", "pmid": "37452120"},
            {"title": "Daver_2019", "doi": "10.1038/s41375-018-0357-9", "pmid": "30622285"},
            {"title": "Rasmussen_2021", "doi": "10.1038/s41467-021-26093-2", "pmid": "34667156"},
            {"title": "Ley_2013", "doi": "10.1056/NEJMoa1301689", "pmid": "23634996"},
            {"title": "Loghavi_2014", "doi": "10.1186/s13045-014-0074-4", "pmid": "25260825"},
        ],
        
        "Drug_Response_Prediction": [
            {"title": "Trac_2023", "doi": "10.1038/s41698-023-00374-z", "pmid": None},
            {"title": "Pino_2024", "doi": "10.1016/j.xcrm.2023.101359", "pmid": "38232702"},
            {"title": "Lee_2018_DRP", "doi": "10.1038/s41467-017-02465-5", "pmid": "29295995"},
            {"title": "Short_2024", "doi": "10.1200/JCO.23.01911", "pmid": "38277619"},
            {"title": "Pollyea_2022", "doi": "10.1158/1078-0432.CCR-21-4119", "pmid": "35046058"},
            {"title": "Lachowiez_2023_DRP", "doi": "10.1158/2643-3230.BCD-22-0205", "pmid": "37102976"},
            {"title": "Zeng_2022", "doi": "10.1038/s41591-022-01819-x", "pmid": "35618838"},
            {"title": "Wang_2025", "doi": "10.1002/advs.202506447", "pmid": None},
        ],
        
        "Clinical_Prognostic_Models": [
            {"title": "Tazi_2022", "doi": "10.1038/s41467-022-32103-8", "pmid": "35948539"},
            {"title": "Rausch_2023", "doi": "10.1038/s41375-023-01884-2", "pmid": None},
            {"title": "Mrozek_2023", "doi": "10.1038/s41375-023-01846-8", "pmid": None},
            {"title": "Kong_2022", "doi": "10.1038/s41375-022-01662-6", "pmid": None},
            {"title": "Lai_2022", "doi": "10.1186/s12911-022-01791-z", "pmid": None},
        ],
        
        "Computational_Methods": [
            {"title": "Zhang_2020", "doi": "10.1093/nargab/lqaa078", "pmid": None},
            {"title": "Behdenna_2023", "doi": "10.1186/s12859-023-05578-5", "pmid": "38057672"},
            {"title": "Briere_2021", "doi": "10.1186/s12859-021-04279-1", "pmid": "34229602"},
            {"title": "Leng_2022", "doi": "10.1186/s13059-022-02739-2", "pmid": "35915444"},
            {"title": "Soneson_2013", "doi": "10.1186/1471-2105-14-91", "pmid": "23497356"},
        ],
        
        "Clinical_Trials": [
            {"title": "Perl_2019", "doi": "10.1056/NEJMoa1902688", "pmid": "31665578"},
            {"title": "Stone_2017", "doi": "10.1056/NEJMoa1614359", "pmid": "28644114"},
            {"title": "DiNardo_2018", "doi": "10.1056/NEJMoa1716984", "pmid": "29860938"},
            {"title": "Stein_2017", "doi": "10.1182/blood-2017-04-779405", "pmid": "28588020"},
            {"title": "DiNardo_2020", "doi": "10.1056/NEJMoa2012971", "pmid": "32786187"},
            {"title": "Daver_2022", "doi": "10.1200/JCO.22.00602", "pmid": "35849791"},
            {"title": "Burd_2020_Trial", "doi": "10.1038/s41591-020-1089-8", "pmid": "33169022"},
            {"title": "Dohner_2021", "doi": "10.1038/s41571-021-00509-w", "pmid": "34006997"},
        ],
    }
    
    return papers


def main():
    """Main execution function"""
    
    print("=" * 70)
    print("AML LITERATURE DOWNLOAD SCRIPT")
    print("=" * 70)
    print()
    
    # Configuration
    OUTPUT_DIR = "AML_Literature"
    YOUR_EMAIL = "bruim2010@gmail.com"  # CHANGE THIS!
    USE_SCIHUB = True  # Set to True to enable Sci-Hub as fallback
    
    print("Configuration:")
    print(f"  Output directory: {OUTPUT_DIR}")
    print(f"  Email (for Unpaywall): {YOUR_EMAIL}")
    print(f"  Use Sci-Hub: {USE_SCIHUB}")
    print()
    
    if YOUR_EMAIL == "bruim2010@gmail.com":
        print("⚠️  WARNING: Please set YOUR_EMAIL variable to your actual email!")
        print("   This is required for Unpaywall API access.")
        print()
    
    if USE_SCIHUB:
        print("⚠️  WARNING: Sci-Hub usage may violate copyright laws in your jurisdiction.")
        print("   Use at your own risk and check your local laws.")
        print()
    
    response = input("Continue? (y/n): ")
    if response.lower() != 'y':
        print("Exiting...")
        return
    
    print("\nStarting downloads...\n")
    
    # Initialize downloader
    downloader = PaperDownloader(output_dir=OUTPUT_DIR, email=YOUR_EMAIL)
    
    # Parse paper list
    papers_by_category = parse_paper_list()
    
    # Download papers by category
    for category, papers in papers_by_category.items():
        print(f"\n{'='*70}")
        print(f"CATEGORY: {category}")
        print(f"{'='*70}\n")
        
        # Create category subdirectory
        category_dir = downloader.output_dir / category
        category_dir.mkdir(exist_ok=True)
        
        # Temporarily change output directory
        original_dir = downloader.output_dir
        downloader.output_dir = category_dir
        
        # Download each paper in category
        for paper in papers:
            downloader.download_paper(
                doi=paper.get('doi'),
                pmid=paper.get('pmid'),
                title=paper['title'],
                use_scihub=USE_SCIHUB
            )
            time.sleep(1)  # Be polite
        
        # Restore original output directory
        downloader.output_dir = original_dir
    
    # Print summary
    downloader.print_summary()
    
    print(f"\nDownloads completed!")
    print(f"Papers saved to: {OUTPUT_DIR}/")
    print(f"Failed downloads listed in: {OUTPUT_DIR}/failed_downloads.txt")
    print(f"Log file: paper_download.log")


if __name__ == "__main__":
    main()