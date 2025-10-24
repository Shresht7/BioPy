"""
Data Collection Script
"""

# Library
import requests
import gzip
import io
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import pandas as pd
import time
import sys

if len(sys.argv) < 2:
    raise RuntimeError("Please provide an email to use for NCBI Entrez")
email = sys.argv[1]

# ------------------
Entrez.email = email
# ------------------

class DataCollector:
    """Collect and prepare training data"""

    def __init__(self, output_dir='data'):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

    def download_card_database(self):
        """Download CARD resistance gene database"""
        print("Downloading CARD database...")

        # CARD nucleotide sequences URL
        url = "https://card.mcmaster.ca/latest/data"

        print("""
            Manual Download Required:
            1. Visit: https://card.mcmaster.ca/download
            2. Click "Download" under "Nucleotide FASTA"
            3. Save as 'data/card_database.fasta'
            
            Alternative: Use wget/curl:
            wget https://card.mcmaster.ca/latest/data -O data/card_database.fasta  
        """)

    def filter_resistance_genes(self, input_fasta='data/card_database.fasta', output_fasta='data/resistance_genes.fasta', resistance_types=None):
        """Filter CARD database for specific resistance genes"""

        if resistance_types is None:
            # Focus on major antibiotic resistance classes
            resistance_types = [
                "BLA", # 'beta-lactam',
                "AMG", # 'aminoglycoside',
                "FLO", # 'fluoroquinolone',
                "TET", # 'tetracycline',
                "MAC", # 'macrolide',
                "SLF", # 'sulfonamide',
                "TMP", # 'trimethoprim'
            ]

        print(f"Filtering resistance genes from {input_fasta}...")
        
        filtered_records = []
        total = 0
        
        try:
            for record in SeqIO.parse(input_fasta, 'fasta'):
                total += 1
                description = record.description.lower()

                # Check if any resistance type is in description
                if any(rt.lower() in description for rt in resistance_types):
                    filtered_records.append(record)

            print(f"Found {len(filtered_records)} resistance genes out of {total} total")

            # Write filtered sequences
            SeqIO.write(filtered_records, output_fasta, 'fasta')
            print('Saved to {output_fasta}')

            return len(filtered_records)
        
        except FileNotFoundError:
            print(f'Error: {input_fasta} not found. Please download CARD database first!')
            return 0
        
    def download_reference_genome(self, accession, organism_name):
        """Download a reference genome from NCBI"""
        print(f"Downloading {organism_name} genome (accession: {accession}...)")

        try:
            # Fetch genome
            handle = Entrez.efetch(db='nucleotide', id=accession, rettype='fasta', retmode='text')
            genome = handle.read()
            handle.close()

            # Save genome
            output_file = self.output_dir / f"{organism_name.replace(' ', '_')}_genome.fasta"
            with open(output_file, 'w') as file:
                file.write(genome)

            print(f"Saved to {output_file}")

            return output_file
        except Exception as e:
            print(f"Error downloading genome of {organism_name}: {e}")
            return None

    def extract_housekeeping_genes(self, genome_file, organism_name):
        """Extract housekeeping genes from a genome"""

        # Common housekeeping genes (universal markers)
        housekeeping_genes = {
            'recA': 'DNA recombination/repair',
            'gyrB': 'DNA gyrase subunit B',
            'rpoB': 'RNA polymerase beta',
            'gyrA': 'DNA gyrase subunit A',
            'dnaK': 'Chaperone protein',
            'groEL': 'Chaperonin',
            'atpD': 'ATP synthase',
            'infB': 'Translation initiation factor'
        }

        print(f"\nExtracting housekeeping genes from {organism_name}...")
        
        native_genes = []

        try:
            # In practice, you'd need gene annotations (GenBank format)
            # This is a simplified example
            for record in SeqIO.parse(genome_file, 'fasta'):
                # Extract first 8 genes as proxies (in real case, use GFF/GenBank)
                seq_str = str(record.seq)

                # Extract genes at regular intervals (simplified approach)
                gene_length = 1000
                for i, (gene_name, description) in enumerate(housekeeping_genes.items()):
                    start = i * gene_length * 2
                    end = start + gene_length

                    if end < len(seq_str):
                        gene_seq = seq_str[start:end]

                        # Create record
                        gene_record = SeqRecord(
                            Seq(gene_seq),
                            id=f"{organism_name}_{gene_name}",
                            description=f"{description} [{organism_name}]"
                        )
                        native_genes.append(gene_record)

                print(f"Extracted {len(native_genes)} housekeeping genes")
                return native_genes 
            
        except Exception as e:
            print(f"Error extracting genes: {e}")
            return []
        
    def create_training_dataset(self):
        """Create complete training dataset"""
        
        print("\n" + "="*70)
        print("CREATING TRAINING DATASET")
        print("="*70)
        
        # Step 1: Download and filter resistance genes
        print("\n[Step 1] Processing resistance genes...")
        self.download_card_database()
        
        if Path('data/card_database.fasta').exists():
            n_resistance = self.filter_resistance_genes()
        else:
            n_resistance = 0
            print("Please download CARD database manually")
        
        # Step 2: Download reference genomes
        print("\n[Step 2] Downloading reference genomes...")
        
        reference_genomes = [
            ('NC_000913.3', 'Escherichia_coli_K12'),
            ('NC_003197.2', 'Salmonella_enterica'),
            ('NC_009648.1', 'Klebsiella_pneumoniae'),
            ('NC_002516.2', 'Pseudomonas_aeruginosa')
        ]
        
        all_native_genes = []
        
        for accession, organism in reference_genomes:
            genome_file = self.download_reference_genome(accession, organism)
            if genome_file:
                native_genes = self.extract_housekeeping_genes(genome_file, organism)
                if native_genes:
                    all_native_genes.extend(native_genes)
            
            # Be nice to NCBI servers
            time.sleep(1)
        
        # Save native genes
        if all_native_genes:
            output_file = self.output_dir / 'native_genes.fasta'
            SeqIO.write(all_native_genes, output_file, "fasta")
            print(f"\n[Step 3] Saved {len(all_native_genes)} native genes to {output_file}")
        
        # Summary
        print("\n" + "="*70)
        print("DATASET SUMMARY")
        print("="*70)
        print(f"Resistance genes (transferred): {n_resistance}")
        print(f"Native genes: {len(all_native_genes)}")
        print(f"Total samples: {n_resistance + len(all_native_genes)}")
        print("="*70)
        
        return n_resistance, len(all_native_genes)

if __name__ == "__main__":
    collector = DataCollector()
    collector.create_training_dataset()
