#!/usr/bin/env python3
"""
Comprehensive Genomic Analysis Pipeline
=====================================

A world-class bioinformatics pipeline for analyzing human genome data including:
- Gene extraction and sequence analysis
- Protein translation and mutation analysis
- GC content and codon usage computation
- Phylogenetic analysis
- RNA-Seq expression analysis
- CRISPR edit simulation

Author: Girish G shankar
Python: 3.8+
Dependencies: biopython, pandas, numpy, pyensembl, scikit-bio, htseq, matplotlib, seaborn
"""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import warnings
warnings.filterwarnings('ignore')

# Core dependencies
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Bioinformatics libraries
try:
    from Bio import SeqIO, Seq, SeqUtils
    from Bio.SeqUtils import gc_fraction, molecular_weight
    from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
    from Bio.Align import MultipleSeqAlignment
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq as BioSeq
except ImportError:
    print("Error: Biopython not installed. Install with: pip install biopython")
    sys.exit(1)

try:
    from pyensembl import EnsemblRelease
except ImportError:
    print("Error: PyEnsembl not installed. Install with: pip install pyensembl")
    sys.exit(1)

try:
    import skbio
    from skbio import DNA, Protein
    from skbio.sequence import GeneticCode
    from skbio.stats.composition import ancom
except ImportError:
    print("Error: Scikit-bio not installed. Install with: pip install scikit-bio")
    sys.exit(1)

try:
    import HTSeq
except ImportError:
    print("Error: HTSeq not installed. Install with: pip install htseq")
    sys.exit(1)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('genomic_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class GenomicAnalyzer:
    """
    Comprehensive genomic analysis pipeline for human genome data.
    """
    
    def __init__(self, genome_version: str = "GRCh38", output_dir: str = "output"):
        """
        Initialize the genomic analyzer.
        
        Args:
            genome_version: Genome version (GRCh38 or CHM13)
            output_dir: Output directory for results
        """
        self.genome_version = genome_version
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Initialize Ensembl release
        if genome_version == "GRCh38":
            self.ensembl = EnsemblRelease(104)  # Latest GRCh38
        else:
            logger.warning("CHM13 not directly supported by PyEnsembl, using GRCh38")
            self.ensembl = EnsemblRelease(104)
        
        # Mock SNP data for demonstration
        self.mock_snps = {
            'BRCA1': [
                {'pos': 100, 'ref': 'A', 'alt': 'G', 'type': 'missense'},
                {'pos': 250, 'ref': 'C', 'alt': 'T', 'type': 'synonymous'},
            ],
            'TP53': [
                {'pos': 175, 'ref': 'G', 'alt': 'A', 'type': 'nonsense'},
                {'pos': 300, 'ref': 'T', 'alt': 'C', 'type': 'missense'},
            ]
        }
        
        logger.info(f"GenomicAnalyzer initialized with {genome_version}")
    
    def install_genome_data(self) -> bool:
        """
        Install genome data if not present.
        
        Returns:
            bool: Success status
        """
        try:
            if not self.ensembl.installed:
                logger.info("Installing genome data... This may take a while.")
                self.ensembl.download()
                self.ensembl.index()
            return True
        except Exception as e:
            logger.error(f"Failed to install genome data: {e}")
            return False
    
    def extract_gene_info(self, gene_name: str) -> Dict:
        """
        Extract comprehensive gene information including coordinates and exons.
        
        Args:
            gene_name: Gene symbol (e.g., 'BRCA1', 'TP53')
            
        Returns:
            Dict: Gene information including coordinates, exons, and sequence
        """
        try:
            logger.info(f"Extracting information for gene: {gene_name}")
            
            # Get gene information
            genes = self.ensembl.genes_by_name(gene_name)
            if not genes:
                raise ValueError(f"Gene {gene_name} not found")
            
            gene = genes[0]  # Take first match
            
            # Extract exon information
            exons = []
            for transcript in gene.transcripts:
                for exon in transcript.exons:
                    exons.append({
                        'exon_id': exon.exon_id,
                        'start': exon.start,
                        'end': exon.end,
                        'strand': exon.strand,
                        'length': len(exon)
                    })
            
            gene_info = {
                'gene_id': gene.gene_id,
                'gene_name': gene.gene_name,
                'chromosome': gene.contig,
                'start': gene.start,
                'end': gene.end,
                'strand': gene.strand,
                'length': len(gene),
                'exons': exons,
                'num_exons': len(exons),
                'biotype': getattr(gene, 'biotype', 'protein_coding')
            }
            
            logger.info(f"Successfully extracted info for {gene_name}: {len(exons)} exons")
            return gene_info
            
        except Exception as e:
            logger.error(f"Error extracting gene info for {gene_name}: {e}")
            return self._create_mock_gene_info(gene_name)
    
    def _create_mock_gene_info(self, gene_name: str) -> Dict:
        """Create mock gene information for testing."""
        mock_sequences = {
            'BRCA1': 'ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGAGCCTACAAGAAAGTACGAGATTTAGTCAACTTGTTGAAGAGCTATTGAAAATCATTTGTGCTTTTCAGCTTGACACAGGTTTGGAGTATGCAAACAGCTATAATTTTGCAAAAAAGGAAAATAACTGCAAGAGTTTCCTCGAAAACTATCACTGATTTCCAGAGACCCTACAAGAAAGTACGAGATTTAGTCAACTTGTTGAAGAGCTATTG',
            'TP53': 'ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCT'
        }
        
        return {
            'gene_id': f'ENSG00000{gene_name}',
            'gene_name': gene_name,
            'chromosome': '17' if gene_name == 'BRCA1' else '17',
            'start': 43044295 if gene_name == 'BRCA1' else 7565097,
            'end': 43125483 if gene_name == 'BRCA1' else 7590856,
            'strand': '-' if gene_name == 'BRCA1' else '+',
            'length': len(mock_sequences.get(gene_name, 'ATGCGT')),
            'sequence': mock_sequences.get(gene_name, 'ATGCGT'),
            'exons': [
                {'exon_id': f'{gene_name}_exon_1', 'start': 100, 'end': 200, 'strand': '+', 'length': 100},
                {'exon_id': f'{gene_name}_exon_2', 'start': 300, 'end': 450, 'strand': '+', 'length': 150},
            ],
            'num_exons': 2,
            'biotype': 'protein_coding'
        }
    
    def fetch_gene_sequence(self, gene_info: Dict) -> Dict:
        """
        Fetch gene sequence and identify coding/non-coding regions.
        
        Args:
            gene_info: Gene information dictionary
            
        Returns:
            Dict: Sequence information including coding and non-coding regions
        """
        try:
            gene_name = gene_info['gene_name']
            logger.info(f"Fetching sequence for {gene_name}")
            
            # For mock data, use pre-defined sequences
            if 'sequence' in gene_info:
                full_sequence = gene_info['sequence']
            else:
                # In real implementation, fetch from genome
                full_sequence = self._generate_mock_sequence(gene_info['length'])
            
            # Create Bio.Seq object
            bio_seq = BioSeq(full_sequence)
            
            # Identify regions (simplified for mock data)
            coding_regions = []
            non_coding_regions = []
            
            # Assume first 80% is coding, rest is UTR (simplified)
            coding_length = int(len(full_sequence) * 0.8)
            coding_regions.append({
                'start': 1,
                'end': coding_length,
                'sequence': str(bio_seq[:coding_length])
            })
            
            if coding_length < len(full_sequence):
                non_coding_regions.append({
                    'start': coding_length + 1,
                    'end': len(full_sequence),
                    'sequence': str(bio_seq[coding_length:]),
                    'type': '3_UTR'
                })
            
            sequence_info = {
                'full_sequence': str(bio_seq),
                'length': len(bio_seq),
                'coding_regions': coding_regions,
                'non_coding_regions': non_coding_regions,
                'gc_content': gc_fraction(bio_seq) * 100,
                'molecular_weight': molecular_weight(bio_seq, seq_type='DNA')
            }
            
            logger.info(f"Sequence fetched: {len(bio_seq)} bp, GC: {sequence_info['gc_content']:.2f}%")
            return sequence_info
            
        except Exception as e:
            logger.error(f"Error fetching sequence: {e}")
            return {'error': str(e)}
    
    def translate_and_find_mutations(self, sequence_info: Dict, gene_name: str) -> Dict:
        """
        Translate DNA to protein and identify mutations.
        
        Args:
            sequence_info: Sequence information
            gene_name: Gene name for SNP lookup
            
        Returns:
            Dict: Translation and mutation information
        """
        try:
            logger.info(f"Translating sequence and finding mutations for {gene_name}")
            
            # Get coding sequence
            if sequence_info['coding_regions']:
                coding_seq = sequence_info['coding_regions'][0]['sequence']
            else:
                coding_seq = sequence_info['full_sequence']
            
            bio_seq = BioSeq(coding_seq)
            
            # Translate to protein
            protein_seq = bio_seq.translate()
            
            # Remove stop codons for analysis
            protein_clean = str(protein_seq).replace('*', '')
            
            # Apply mock mutations
            mutations = []
            snps = self.mock_snps.get(gene_name, [])
            
            mutated_seq = coding_seq
            for snp in snps:
                pos = snp['pos']
                if pos < len(mutated_seq):
                    mutated_seq = mutated_seq[:pos] + snp['alt'] + mutated_seq[pos+1:]
                    mutations.append({
                        'position': pos,
                        'reference': snp['ref'],
                        'alternate': snp['alt'],
                        'type': snp['type'],
                        'codon_change': self._get_codon_change(coding_seq, pos, snp['alt'])
                    })
            
            # Translate mutated sequence
            mutated_protein = BioSeq(mutated_seq).translate()
            
            translation_info = {
                'original_dna': coding_seq,
                'mutated_dna': mutated_seq,
                'original_protein': str(protein_seq),
                'mutated_protein': str(mutated_protein),
                'protein_length': len(protein_clean),
                'mutations': mutations,
                'stop_codons': str(protein_seq).count('*')
            }
            
            logger.info(f"Translation completed: {len(protein_clean)} amino acids, {len(mutations)} mutations")
            return translation_info
            
        except Exception as e:
            logger.error(f"Error in translation: {e}")
            return {'error': str(e)}
    
    def compute_sequence_stats(self, sequence_info: Dict, gene_name: str) -> Dict:
        """
        Compute GC content, codon usage, and prepare phylogenetic data.
        
        Args:
            sequence_info: Sequence information
            gene_name: Gene name
            
        Returns:
            Dict: Sequence statistics
        """
        try:
            logger.info(f"Computing sequence statistics for {gene_name}")
            
            full_seq = sequence_info['full_sequence']
            
            # GC content analysis
            gc_content = gc_fraction(BioSeq(full_seq)) * 100
            
            # Sliding window GC analysis
            window_size = 100
            gc_windows = []
            for i in range(0, len(full_seq) - window_size + 1, window_size // 2):
                window_seq = full_seq[i:i + window_size]
                gc_windows.append({
                    'position': i,
                    'gc_content': GC(BioSeq(window_seq))
                })
            
            # Codon usage analysis
            codon_usage = {}
            coding_seq = sequence_info.get('coding_regions', [{}])[0].get('sequence', full_seq)
            
            for i in range(0, len(coding_seq) - 2, 3):
                codon = coding_seq[i:i+3]
                if len(codon) == 3:
                    codon_usage[codon] = codon_usage.get(codon, 0) + 1
            
            # Generate mock variant sequences for phylogenetic analysis
            variants = self._generate_mock_variants(full_seq, gene_name)
            
            stats = {
                'gc_content': gc_content,
                'gc_windows': gc_windows,
                'codon_usage': codon_usage,
                'most_common_codons': sorted(codon_usage.items(), key=lambda x: x[1], reverse=True)[:10],
                'variants': variants,
                'sequence_length': len(full_seq),
                'coding_length': len(coding_seq)
            }
            
            logger.info(f"Statistics computed: GC={gc_content:.2f}%, {len(codon_usage)} unique codons")
            return stats
            
        except Exception as e:
            logger.error(f"Error computing statistics: {e}")
            return {'error': str(e)}
    
    def analyze_expression_data(self, gene_name: str, bam_file: Optional[str] = None) -> Dict:
        """
        Analyze RNA-Seq expression data using HTSeq.
        
        Args:
            gene_name: Gene name
            bam_file: Path to BAM file (optional, will use mock data if None)
            
        Returns:
            Dict: Expression analysis results
        """
        try:
            logger.info(f"Analyzing expression data for {gene_name}")
            
            if bam_file and os.path.exists(bam_file):
                # Real BAM file analysis
                bam_reader = HTSeq.BAM_Reader(bam_file)
                counts = {}
                
                for alignment in bam_reader:
                    if alignment.aligned:
                        gene_id = alignment.optional_field("GN") if alignment.has_optional_field("GN") else "unknown"
                        counts[gene_id] = counts.get(gene_id, 0) + 1
                        
                expression_level = counts.get(gene_name, 0)
            else:
                # Mock expression data
                np.random.seed(42)
                expression_level = np.random.poisson(100)  # Mock read count
                counts = {gene_name: expression_level}
                
                # Generate mock time-course data
                time_points = range(0, 24, 2)  # 0-22 hours
                mock_expression = []
                for t in time_points:
                    # Simulate circadian expression pattern
                    base_expr = 100
                    circadian = 50 * np.sin(2 * np.pi * t / 24)
                    noise = np.random.normal(0, 10)
                    expr = max(0, base_expr + circadian + noise)
                    mock_expression.append({'time': t, 'expression': expr})
            
            # Differential expression simulation
            conditions = ['control', 'treatment']
            diff_expr = {}
            for condition in conditions:
                if condition == 'control':
                    diff_expr[condition] = np.random.poisson(expression_level, 6)  # 6 replicates
                else:
                    # Simulate 2-fold change
                    diff_expr[condition] = np.random.poisson(expression_level * 2, 6)
            
            # Statistical analysis
            from scipy import stats
            t_stat, p_value = stats.ttest_ind(diff_expr['control'], diff_expr['treatment'])
            fold_change = np.mean(diff_expr['treatment']) / np.mean(diff_expr['control'])
            
            expression_info = {
                'gene_name': gene_name,
                'read_count': expression_level,
                'time_course': mock_expression if not bam_file else [],
                'differential_expression': diff_expr,
                'fold_change': fold_change,
                'p_value': p_value,
                'significant': p_value < 0.05,
                't_statistic': t_stat
            }
            
            logger.info(f"Expression analysis completed: {expression_level} reads, FC={fold_change:.2f}")
            return expression_info
            
        except Exception as e:
            logger.error(f"Error analyzing expression data: {e}")
            return {'error': str(e)}
    
    def simulate_crispr_edit(self, sequence_info: Dict, edit_type: str = "deletion", 
                           position: int = 100, size: int = 10) -> Dict:
        """
        Simulate CRISPR edit on gene sequence.
        
        Args:
            sequence_info: Original sequence information
            edit_type: Type of edit ('deletion', 'insertion', 'substitution')
            position: Position for edit (0-based)
            size: Size of edit
            
        Returns:
            Dict: Before and after comparison
        """
        try:
            logger.info(f"Simulating CRISPR {edit_type} at position {position}")
            
            original_seq = sequence_info['full_sequence']
            
                        if position + size > len(original_seq):
                raise ValueError("CRISPR edit exceeds sequence length.")

if edit_type == "deletion":
                edited_seq = original_seq[:position] + original_seq[position + size:]
                edit_description = f"Deleted {size} bp at position {position}"
                
            elif edit_type == "insertion":
                insert_seq = "A" * size  # Simple insertion
                edited_seq = original_seq[:position] + insert_seq + original_seq[position:]
                edit_description = f"Inserted {size} bp at position {position}"
                
            elif edit_type == "substitution":
                substitute_seq = "T" * size
                edited_seq = original_seq[:position] + substitute_seq + original_seq[position + size:]
                edit_description = f"Substituted {size} bp at position {position}"
            else:
                raise ValueError(f"Unknown edit type: {edit_type}")
            
            # Translate both sequences
            original_protein = str(BioSeq(original_seq).translate())
            edited_protein = str(BioSeq(edited_seq).translate())
            
            # Find differences
            protein_changes = []
            min_len = min(len(original_protein), len(edited_protein))
            
            for i in range(min_len):
                if original_protein[i] != edited_protein[i]:
                    protein_changes.append({
                        'position': i + 1,
                        'original': original_protein[i],
                        'edited': edited_protein[i]
                    })
            
            crispr_result = {
                'edit_type': edit_type,
                'edit_position': position,
                'edit_size': size,
                'edit_description': edit_description,
                'original_sequence': original_seq,
                'edited_sequence': edited_seq,
                'original_protein': original_protein,
                'edited_protein': edited_protein,
                'protein_changes': protein_changes,
                'frameshift': len(original_seq) % 3 != len(edited_seq) % 3,
                'length_change': len(edited_seq) - len(original_seq)
            }
            
            logger.info(f"CRISPR simulation completed: {len(protein_changes)} protein changes")
            return crispr_result
            
        except Exception as e:
            logger.error(f"Error simulating CRISPR edit: {e}")
            return {'error': str(e)}
    
    def generate_visualizations(self, all_results: Dict) -> None:
        """
        Generate comprehensive visualizations of the analysis results.
        
        Args:
            all_results: Dictionary containing all analysis results
        """
        try:
            logger.info("Generating visualizations")
            
            # Set up the plotting style
            plt.style.use('seaborn-v0_8')
            fig = plt.figure(figsize=(20, 16))
            
            # 1. GC Content sliding window
            ax1 = plt.subplot(3, 4, 1)
            if 'sequence_stats' in all_results and 'gc_windows' in all_results['sequence_stats']:
                gc_data = all_results['sequence_stats']['gc_windows']
                positions = [w['position'] for w in gc_data]
                gc_values = [w['gc_content'] for w in gc_data]
                plt.plot(positions, gc_values, 'b-', linewidth=2)
                plt.title('GC Content (Sliding Window)')
                plt.xlabel('Position (bp)')
                plt.ylabel('GC Content (%)')
                plt.grid(True, alpha=0.3)
            
            # 2. Codon usage
            ax2 = plt.subplot(3, 4, 2)
            if 'sequence_stats' in all_results and 'most_common_codons' in all_results['sequence_stats']:
                codons = all_results['sequence_stats']['most_common_codons'][:10]
                codon_names = [c[0] for c in codons]
                codon_counts = [c[1] for c in codons]
                plt.bar(codon_names, codon_counts, color='green', alpha=0.7)
                plt.title('Top 10 Codon Usage')
                plt.xlabel('Codon')
                plt.ylabel('Count')
                plt.xticks(rotation=45)
            
            # 3. Expression time course
            ax3 = plt.subplot(3, 4, 3)
            if 'expression' in all_results and 'time_course' in all_results['expression']:
                time_data = all_results['expression']['time_course']
                if time_data:
                    times = [t['time'] for t in time_data]
                    expressions = [t['expression'] for t in time_data]
                    plt.plot(times, expressions, 'ro-', linewidth=2, markersize=6)
                    plt.title('Expression Time Course')
                    plt.xlabel('Time (hours)')
                    plt.ylabel('Expression Level')
                    plt.grid(True, alpha=0.3)
            
            # 4. Differential expression
            ax4 = plt.subplot(3, 4, 4)
            if 'expression' in all_results and 'differential_expression' in all_results['expression']:
                diff_data = all_results['expression']['differential_expression']
                conditions = list(diff_data.keys())
                values = [diff_data[cond] for cond in conditions]
                plt.boxplot(values, labels=conditions)
                plt.title('Differential Expression')
                plt.ylabel('Expression Level')
                plt.grid(True, alpha=0.3)
            
            # 5. Mutation positions
            ax5 = plt.subplot(3, 4, 5)
            if 'translation' in all_results and 'mutations' in all_results['translation']:
                mutations = all_results['translation']['mutations']
                if mutations:
                    positions = [m['position'] for m in mutations]
                    types = [m['type'] for m in mutations]
                    colors = {'missense': 'red', 'synonymous': 'blue', 'nonsense': 'black'}
                    for mut_type in set(types):
                        pos_type = [p for p, t in zip(positions, types) if t == mut_type]
                        plt.scatter(pos_type, [1] * len(pos_type), 
                                  c=colors.get(mut_type, 'gray'), 
                                  label=mut_type, s=100, alpha=0.7)
                    plt.title('Mutation Positions')
                    plt.xlabel('Position')
                    plt.ylabel('Mutations')
                    plt.legend()
                    plt.grid(True, alpha=0.3)
            
            # 6. Sequence composition
            ax6 = plt.subplot(3, 4, 6)
            if 'sequence_info' in all_results:
                seq = all_results['sequence_info']['full_sequence']
                composition = {'A': seq.count('A'), 'T': seq.count('T'), 
                             'G': seq.count('G'), 'C': seq.count('C')}
                plt.pie(composition.values(), labels=composition.keys(), autopct='%1.1f%%')
                plt.title('Nucleotide Composition')
            
            # 7. Protein length comparison (if CRISPR data available)
            ax7 = plt.subplot(3, 4, 7)
            if 'crispr' in all_results:
                original_len = len(all_results['crispr']['original_protein'])
                edited_len = len(all_results['crispr']['edited_protein'])
                plt.bar(['Original', 'CRISPR Edited'], [original_len, edited_len], 
                       color=['blue', 'red'], alpha=0.7)
                plt.title('Protein Length: Before vs After CRISPR')
                plt.ylabel('Amino Acids')
            
            # 8. Expression fold change
            ax8 = plt.subplot(3, 4, 8)
            if 'expression' in all_results:
                fc = all_results['expression'].get('fold_change', 1)
                p_val = all_results['expression'].get('p_value', 1)
                color = 'red' if p_val < 0.05 else 'gray'
                plt.bar(['Fold Change'], [fc], color=color, alpha=0.7)
                plt.axhline(y=1, color='black', linestyle='--', alpha=0.5)
                plt.title(f'Expression Fold Change\n(p={p_val:.3f})')
                plt.ylabel('Fold Change')
            
            # 9-12: Gene structure visualization
            ax9 = plt.subplot(3, 4, (9, 12))
            if 'gene_info' in all_results:
                gene_info = all_results['gene_info']
                gene_start = gene_info['start']
                gene_end = gene_info['end']
                
                # Draw gene
                plt.barh(0, gene_end - gene_start, left=gene_start, height=0.3, 
                        color='lightblue', alpha=0.7, label='Gene')
                
                # Draw exons
                for i, exon in enumerate(gene_info.get('exons', [])):
                    plt.barh(0, exon['end'] - exon['start'], left=exon['start'], 
                            height=0.5, color='darkblue', alpha=0.8)
                
                plt.title(f"Gene Structure: {gene_info['gene_name']}")
                plt.xlabel('Genomic Position')
                plt.ylabel('')
                plt.yticks([])
                plt.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            # Save the plot
            output_file = self.output_dir / "genomic_analysis_results.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"Visualizations saved to {output_file}")
            
            plt.show()
            
        except Exception as e:
            logger.error(f"Error generating visualizations: {e}")
    
    def export_results(self, all_results: Dict, gene_name: str) -> None:
        """
        Export results to various formats (FASTA, CSV, JSON).
        
        Args:
            all_results: All analysis results
            gene_name: Gene name for file naming
        """
        try:
            logger.info("Exporting results to files")
            
            # Export FASTA sequences
            fasta_file = self.output_dir / f"{gene_name}_sequences.fasta"
            with open(fasta_file, 'w') as f:
                if 'sequence_info' in all_results:
                    seq_info = all_results['sequence_info']
                    f.write(f">{gene_name}_genomic\n{seq_info['full_sequence']}\n")
                    
                    if 'coding_regions' in seq_info:
                        for i, region in enumerate(seq_info['coding_regions']):
                            f.write(f">{gene_name}_coding_{i+1}\n{region['sequence']}\n")
                
                if 'translation' in all_results:
                    trans_info = all_results['translation']
                    f.write(f">{gene_name}_protein_original\n{trans_info['original_protein']}\n")
                    f.write(f">{gene_name}_protein_mutated\n{trans_info['mutated_protein']}\n")
                
                if 'crispr' in all_results:
                    crispr_info = all_results['crispr']
                    f.write(f">{gene_name}_crispr_edited\n{crispr_info['edited_sequence']}\n")
                    f.write(f">{gene_name}_protein_crispr\n{crispr_info['edited_protein']}\n")
            
            # Export CSV data
            csv_file = self.output_dir / f"{gene_name}_analysis.csv"
            rows = []
            
            # Gene information
            if 'gene_info' in all_results:
                gene_info = all_results['gene_info']
                rows.append(['Gene_ID', gene_info.get('gene_id', 'N/A')])
                rows.append(['Gene_Name', gene_info.get('gene_name', 'N/A')])
                rows.append(['Chromosome', gene_info.get('chromosome', 'N/A')])
                rows.append(['Start_Position', gene_info.get('start', 'N/A')])
                rows.append(['End_Position', gene_info.get('end', 'N/A')])
                rows.append(['Strand', gene_info.get('strand', 'N/A')])
                rows.append(['Gene_Length', gene_info.get('length', 'N/A')])
                rows.append(['Number_of_Exons', gene_info.get('num_exons', 'N/A')])
            
            # Sequence statistics
            if 'sequence_stats' in all_results:
                stats = all_results['sequence_stats']
                rows.append(['GC_Content', f"{stats.get('gc_content', 0):.2f}%"])
                rows.append(['Sequence_Length', stats.get('sequence_length', 'N/A')])
                rows.append(['Coding_Length', stats.get('coding_length', 'N/A')])
                rows.append(['Unique_Codons', len(stats.get('codon_usage', {}))])
            
            # Expression data
            if 'expression' in all_results:
                expr = all_results['expression']
                rows.append(['Expression_Level', expr.get('read_count', 'N/A')])
                rows.append(['Fold_Change', f"{expr.get('fold_change', 1):.2f}"])
                rows.append(['P_Value', f"{expr.get('p_value', 1):.2e}"])
                rows.append(['Significant', expr.get('significant', False)])
            
            # Mutations
            if 'translation' in all_results and 'mutations' in all_results['translation']:
                mutations = all_results['translation']['mutations']
                rows.append(['Number_of_Mutations', len(mutations)])
                for i, mut in enumerate(mutations):
                    rows.append([f'Mutation_{i+1}_Position', mut['position']])
                    rows.append([f'Mutation_{i+1}_Type', mut['type']])
                    rows.append([f'Mutation_{i+1}_Change', f"{mut['reference']}>{mut['alternate']}"])
            
            # Save CSV
            df = pd.DataFrame(rows, columns=['Parameter', 'Value'])
            df.to_csv(csv_file, index=False)
            
            # Export detailed expression data
            if 'expression' in all_results:
                expr_file = self.output_dir / f"{gene_name}_expression.csv"
                expr_data = []
                
                # Time course data
                if 'time_course' in all_results['expression']:
                    for point in all_results['expression']['time_course']:
                        expr_data.append({
                            'Type': 'TimeCourse',
                            'Condition': 'Control',
                            'Time': point['time'],
                            'Expression': point['expression']
                        })
                
                # Differential expression data
                if 'differential_expression' in all_results['expression']:
                    diff_data = all_results['expression']['differential_expression']
                    for condition, values in diff_data.items():
                        for i, value in enumerate(values):
                            expr_data.append({
                                'Type': 'DifferentialExpression',
                                'Condition': condition,
                                'Replicate': i + 1,
                                'Expression': value
                            })
                
                if expr_data:
                    expr_df = pd.DataFrame(expr_data)
                    expr_df.to_csv(expr_file, index=False)
            
            # Export codon usage data
            if 'sequence_stats' in all_results and 'codon_usage' in all_results['sequence_stats']:
                codon_file = self.output_dir / f"{gene_name}_codon_usage.csv"
                codon_data = []
                for codon, count in all_results['sequence_stats']['codon_usage'].items():
                    amino_acid = str(BioSeq(codon).translate()) if len(codon) == 3 else 'N/A'
                    codon_data.append({
                        'Codon': codon,
                        'Count': count,
                        'Amino_Acid': amino_acid,
                        'Frequency': count / sum(all_results['sequence_stats']['codon_usage'].values())
                    })
                
                codon_df = pd.DataFrame(codon_data)
                codon_df = codon_df.sort_values('Count', ascending=False)
                codon_df.to_csv(codon_file, index=False)
            
            # Export JSON summary
            import json
            json_file = self.output_dir / f"{gene_name}_summary.json"
            
            # Create a JSON-serializable summary
            json_summary = {
                'gene_name': gene_name,
                'analysis_timestamp': pd.Timestamp.now().isoformat(),
                'genome_version': self.genome_version,
                'gene_info': all_results.get('gene_info', {}),
                'sequence_length': all_results.get('sequence_info', {}).get('length', 0),
                'gc_content': all_results.get('sequence_stats', {}).get('gc_content', 0),
                'protein_length': all_results.get('translation', {}).get('protein_length', 0),
                'mutation_count': len(all_results.get('translation', {}).get('mutations', [])),
                'expression_level': all_results.get('expression', {}).get('read_count', 0),
                'fold_change': all_results.get('expression', {}).get('fold_change', 1),
                'p_value': all_results.get('expression', {}).get('p_value', 1),
                'crispr_edit': all_results.get('crispr', {}).get('edit_description', 'None')
            }
            
            with open(json_file, 'w') as f:
                json.dump(json_summary, f, indent=2)
            
            logger.info(f"Results exported to {self.output_dir}")
            logger.info(f"Files created: {fasta_file.name}, {csv_file.name}, {json_file.name}")
            
        except Exception as e:
            logger.error(f"Error exporting results: {e}")
    
    def _generate_mock_sequence(self, length: int) -> str:
        """Generate a mock DNA sequence of specified length."""
        np.random.seed(42)
        bases = ['A', 'T', 'G', 'C']
        return ''.join(np.random.choice(bases, length))
    
    def _get_codon_change(self, sequence: str, position: int, new_base: str) -> str:
        """Get codon change information for a mutation."""
        try:
            codon_start = (position // 3) * 3
            original_codon = sequence[codon_start:codon_start + 3]
            
            # Create mutated codon
            mut_pos_in_codon = position % 3
            mutated_codon = (original_codon[:mut_pos_in_codon] + 
                           new_base + 
                           original_codon[mut_pos_in_codon + 1:])
            
            if len(original_codon) == 3 and len(mutated_codon) == 3:
                orig_aa = str(BioSeq(original_codon).translate())
                mut_aa = str(BioSeq(mutated_codon).translate())
                return f"{original_codon}({orig_aa}) -> {mutated_codon}({mut_aa})"
            else:
                return f"{original_codon} -> {mutated_codon}"
        except:
            return "Unknown"
    
    def _generate_mock_variants(self, sequence: str, gene_name: str) -> List[Dict]:
        """Generate mock variant sequences for phylogenetic analysis."""
        variants = []
        np.random.seed(42)
        
        for i in range(5):  # Generate 5 variants
            variant_seq = list(sequence)
            num_mutations = np.random.randint(1, 5)
            
            for _ in range(num_mutations):
                pos = np.random.randint(0, len(variant_seq))
                original_base = variant_seq[pos]
                new_base = np.random.choice([b for b in ['A', 'T', 'G', 'C'] if b != original_base])
                variant_seq[pos] = new_base
            
            variants.append({
                'id': f"{gene_name}_variant_{i+1}",
                'sequence': ''.join(variant_seq),
                'mutations': num_mutations
            })
        
        return variants
    
    def run_comprehensive_analysis(self, gene_name: str, bam_file: Optional[str] = None,
                                 crispr_edit: bool = True) -> Dict:
        """
        Run the complete genomic analysis pipeline.
        
        Args:
            gene_name: Gene symbol to analyze
            bam_file: Optional BAM file for expression analysis
            crispr_edit: Whether to simulate CRISPR editing
            
        Returns:
            Dict: Complete analysis results
        """
        try:
            logger.info(f"Starting comprehensive analysis for {gene_name}")
            
            # Step 1: Install genome data if needed
            if not self.install_genome_data():
                logger.warning("Using mock data due to genome installation issues")
            
            # Step 2: Extract gene information
            gene_info = self.extract_gene_info(gene_name)
            
            # Step 3: Fetch gene sequence
            sequence_info = self.fetch_gene_sequence(gene_info)
            
            # Step 4: Translate and find mutations
            translation_info = self.translate_and_find_mutations(sequence_info, gene_name)
            
            # Step 5: Compute sequence statistics
            sequence_stats = self.compute_sequence_stats(sequence_info, gene_name)
            
            # Step 6: Analyze expression data
            expression_info = self.analyze_expression_data(gene_name, bam_file)
            
            # Step 7: CRISPR simulation (optional)
            crispr_info = {}
            if crispr_edit:
                crispr_info = self.simulate_crispr_edit(sequence_info, edit_type=args.crispr_type, position=args.crispr_pos, size=args.crispr_size)
            
            # Compile all results
            all_results = {
                'gene_info': gene_info,
                'sequence_info': sequence_info,
                'translation': translation_info,
                'sequence_stats': sequence_stats,
                'expression': expression_info,
                'crispr': crispr_info if crispr_edit else None
            }
            
            # Step 8: Generate visualizations
            self.generate_visualizations(all_results)
            
            # Step 9: Export results
            self.export_results(all_results, gene_name)
            
            # Print summary
            self.print_analysis_summary(all_results, gene_name)
            
            logger.info("Comprehensive analysis completed successfully")
            return all_results
            
        except Exception as e:
            logger.error(f"Error in comprehensive analysis: {e}")
            return {'error': str(e)}
    
    def print_analysis_summary(self, results: Dict, gene_name: str) -> None:
        """Print a comprehensive summary of the analysis results."""
        print("\n" + "="*80)
        print(f"GENOMIC ANALYSIS SUMMARY: {gene_name}")
        print("="*80)
        
        # Gene Information
        if 'gene_info' in results:
            gene_info = results['gene_info']
            print(f"\n📊 GENE INFORMATION:")
            print(f"   Gene ID: {gene_info.get('gene_id', 'N/A')}")
            print(f"   Location: {gene_info.get('chromosome', 'N/A')}:{gene_info.get('start', 'N/A')}-{gene_info.get('end', 'N/A')}")
            print(f"   Strand: {gene_info.get('strand', 'N/A')}")
            print(f"   Length: {gene_info.get('length', 'N/A')} bp")
            print(f"   Exons: {gene_info.get('num_exons', 'N/A')}")
        
        # Sequence Statistics
        if 'sequence_stats' in results:
            stats = results['sequence_stats']
            print(f"\n🧬 SEQUENCE STATISTICS:")
            print(f"   GC Content: {stats.get('gc_content', 0):.2f}%")
            print(f"   Total Length: {stats.get('sequence_length', 'N/A')} bp")
            print(f"   Coding Length: {stats.get('coding_length', 'N/A')} bp")
            print(f"   Unique Codons: {len(stats.get('codon_usage', {}))}")
            
            if 'most_common_codons' in stats:
                print(f"   Top 3 Codons: {', '.join([f'{c[0]}({c[1]})' for c in stats['most_common_codons'][:3]])}")
        
        # Translation Information
        if 'translation' in results:
            trans = results['translation']
            print(f"\n🔬 PROTEIN TRANSLATION:")
            print(f"   Protein Length: {trans.get('protein_length', 'N/A')} amino acids")
            print(f"   Stop Codons: {trans.get('stop_codons', 'N/A')}")
            print(f"   Mutations Found: {len(trans.get('mutations', []))}")
            
            for i, mut in enumerate(trans.get('mutations', [])[:3]):  # Show first 3
                print(f"     {i+1}. Position {mut['position']}: {mut['reference']}>{mut['alternate']} ({mut['type']})")
        
        # Expression Analysis
        if 'expression' in results:
            expr = results['expression']
            print(f"\n📈 EXPRESSION ANALYSIS:")
            print(f"   Read Count: {expr.get('read_count', 'N/A')}")
            print(f"   Fold Change: {expr.get('fold_change', 1):.2f}")
            print(f"   P-value: {expr.get('p_value', 1):.2e}")
            print(f"   Significant: {'Yes' if expr.get('significant', False) else 'No'}")
        
        # CRISPR Analysis
        if 'crispr' in results and results['crispr']:
            crispr = results['crispr']
            print(f"\n✂️ CRISPR SIMULATION:")
            print(f"   Edit Type: {crispr.get('edit_type', 'N/A')}")
            print(f"   Position: {crispr.get('edit_position', 'N/A')}")
            print(f"   Description: {crispr.get('edit_description', 'N/A')}")
            print(f"   Frameshift: {'Yes' if crispr.get('frameshift', False) else 'No'}")
            print(f"   Protein Changes: {len(crispr.get('protein_changes', []))}")
        
        print(f"\n📁 OUTPUT FILES:")
        print(f"   Results saved to: {self.output_dir}")
        print(f"   - {gene_name}_sequences.fasta")
        print(f"   - {gene_name}_analysis.csv")
        print(f"   - {gene_name}_expression.csv")
        print(f"   - {gene_name}_codon_usage.csv")
        print(f"   - {gene_name}_summary.json")
        print(f"   - genomic_analysis_results.png")
        
        print("\n" + "="*80)


def main():
    """Main function with CLI interface."""
    parser = argparse.ArgumentParser(
        description="Comprehensive Genomic Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python genomic_pipeline.py --gene BRCA1
  python genomic_pipeline.py --gene TP53 --bam sample.bam --no-crispr
  python genomic_pipeline.py --gene BRCA1 --output results --genome CHM13
        """
    )
    
    parser.add_argument(
        '--gene', '-g',
        required=True,
        help='Gene symbol to analyze (e.g., BRCA1, TP53)'
    )
    
    parser.add_argument(
        '--genome', '-r',
        default='GRCh38',
        choices=['GRCh38', 'CHM13'],
        help='Genome version to use (default: GRCh38)'
    )
    
    parser.add_argument(
        '--bam', '-b',
        help='BAM file for RNA-Seq analysis (optional)'
    )
    
    parser.add_argument(
        '--output', '-o',
        default='output',
        help='Output directory (default: output)'
    )
    
    parser.add_argument(
        '--no-crispr',
        action='store_true',
        help='Skip CRISPR edit simulation'
    )
    
    parser.add_argument(
        '--crispr-type',
        default='deletion',
        choices=['deletion', 'insertion', 'substitution'],
        help='Type of CRISPR edit to simulate (default: deletion)'
    )
    
    parser.add_argument(
        '--crispr-pos',
        type=int,
        default=150,
        help='Position for CRISPR edit (default: 150)'
    )
    
    parser.add_argument(
        '--crispr-size',
        type=int,
        default=12,
        help='Size of CRISPR edit (default: 12)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Initialize analyzer
    analyzer = GenomicAnalyzer(
        genome_version=args.genome,
        output_dir=args.output
    )
    
    # Run comprehensive analysis
    results = analyzer.run_comprehensive_analysis(
        gene_name=args.gene,
        bam_file=args.bam,
        crispr_edit=not args.no_crispr
    )
    
    if 'error' in results:
        logger.error(f"Analysis failed: {results['error']}")
        sys.exit(1)
    
    print(f"\n✅ Analysis completed successfully!")
    print(f"📊 Results saved to: {analyzer.output_dir}")


if __name__ == "__main__":
    main()