"""Genome utilities."""

import tempfile
import os
from pathlib import Path
from typing import Optional, Dict

from src.utils.bedtools import run_bedtools_closest


def find_nearest_gene(
    chrom: str,
    start: int,
    end: int,
    annotation_path: Path,
) -> Optional[Dict]:
    """
    Find nearest gene to a genomic interval.
    
    Args:
        chrom: Chromosome
        start: Start position
        end: End position
        annotation_path: Path to annotation file
        
    Returns:
        Dict with gene info or None if not found
    """
    # Create temporary BED file for query
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
        query_path = f.name
        f.write(f"{chrom}\t{start}\t{end}\tquery\t0\t+\n")
    
    try:
        # Convert annotation to BED if needed
        suffix = annotation_path.suffix.lower()
        if suffix in ('.gtf', '.gff', '.gff3'):
            with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
                annot_bed_path = f.name
                _convert_gtf_to_genes_bed(annotation_path, annot_bed_path)
        else:
            annot_bed_path = str(annotation_path)
        
        try:
            # Run bedtools closest
            result = run_bedtools_closest(query_path, annot_bed_path, d=True)
            
            # Parse result
            for line in result.stdout.strip().split('\n'):
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 13:
                    distance = int(parts[-1])
                    gene_attrs = parts[11] if len(parts) > 11 else ""
                    gene_name = _parse_gene_name(gene_attrs)
                    gene_id = _parse_gene_id(gene_attrs)
                    
                    return {
                        "gene_name": gene_name,
                        "gene_id": gene_id,
                        "distance": distance,
                    }
            
            return None
            
        finally:
            if annot_bed_path != str(annotation_path):
                os.unlink(annot_bed_path)
                
    finally:
        os.unlink(query_path)


def _convert_gtf_to_genes_bed(gtf_path: Path, output_path: str) -> None:
    """Convert GTF to BED with gene entries only."""
    with open(gtf_path) as f_in, open(output_path, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'gene':
                continue
            
            chrom = parts[0]
            start = int(parts[3]) - 1
            end = int(parts[4])
            strand = parts[6]
            attrs = parts[8]
            
            f_out.write(f"{chrom}\t{start}\t{end}\tgene\t0\t{strand}\t{attrs}\n")


def _parse_gene_name(attrs: str) -> Optional[str]:
    """Extract gene_name from GTF attributes."""
    for part in attrs.split(';'):
        part = part.strip()
        if 'gene_name' in part:
            if '"' in part:
                return part.split('"')[1]
            else:
                return part.split()[-1]
    return None


def _parse_gene_id(attrs: str) -> Optional[str]:
    """Extract gene_id from GTF attributes."""
    for part in attrs.split(';'):
        part = part.strip()
        if 'gene_id' in part and 'gene_name' not in part:
            if '"' in part:
                return part.split('"')[1]
            else:
                return part.split()[-1]
    return None
