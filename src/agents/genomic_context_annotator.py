"""Genomic context annotator using bedtools."""

import subprocess
import tempfile
import os
from pathlib import Path
from typing import List, Optional

from src.models.schema import OffTargetSite, AnnotatedSite, RegionType
from src.utils.bedtools import run_bedtools_intersect
from src.utils.genome import find_nearest_gene


class GenomicContextAnnotator:
    """Annotates off-target sites with genomic context."""
    
    def __init__(self, promoter_upstream: int = 2000, promoter_downstream: int = 500):
        """
        Initialize annotator.
        
        Args:
            promoter_upstream: Bases upstream of TSS to consider as promoter
            promoter_downstream: Bases downstream of TSS to consider as promoter
        """
        self.promoter_upstream = promoter_upstream
        self.promoter_downstream = promoter_downstream
    
    def annotate(
        self,
        sites: List[OffTargetSite],
        annotation_path: Path,
    ) -> List[AnnotatedSite]:
        """
        Annotate off-target sites with genomic context.
        
        Args:
            sites: List of off-target sites
            annotation_path: Path to annotation file (GTF/GFF/BED)
            
        Returns:
            List of annotated sites
        """
        # Create temporary BED file from sites
        with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
            bed_path = f.name
            for site in sites:
                f.write(f"{site.chrom}\t{site.start}\t{site.end}\t{site.name}\t"
                       f"{site.mismatch_count}\t{site.strand}\n")
        
        try:
            # Run bedtools intersect for exons
            exon_hits = self._intersect_with_feature(bed_path, annotation_path, "exon")
            
            # Run bedtools intersect for genes (to identify introns)
            gene_hits = self._intersect_with_feature(bed_path, annotation_path, "gene")
            
            # Annotate each site
            annotated = []
            for site in sites:
                ann = self._annotate_site(site, exon_hits, gene_hits, annotation_path)
                annotated.append(ann)
            
            return annotated
            
        finally:
            os.unlink(bed_path)
    
    def _intersect_with_feature(
        self,
        bed_path: str,
        annotation_path: Path,
        feature_type: str,
    ) -> dict:
        """
        Run bedtools intersect for specific feature type.
        
        Returns:
            Dict mapping site name to intersection info
        """
        # Filter annotation for feature type if GTF
        suffix = annotation_path.suffix.lower()
        
        if suffix in ('.gtf', '.gff', '.gff3'):
            # Create filtered temp file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as f:
                filtered_path = f.name
                self._convert_gtf_to_bed(annotation_path, filtered_path, feature_type)
        else:
            filtered_path = str(annotation_path)
        
        try:
            # Run intersection with full output
            result = run_bedtools_intersect(
                bed_path, 
                filtered_path,
                wa=True,  # Write original A entry
                wb=True,  # Write original B entry
            )
            
            # Parse results
            hits = {}
            for line in result.stdout.strip().split('\n'):
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 12:
                    site_name = parts[3]
                    gene_attrs = parts[11] if len(parts) > 11 else ""
                    gene_name = self._parse_gene_name(gene_attrs)
                    gene_id = self._parse_gene_id(gene_attrs)
                    hits[site_name] = {
                        "gene_name": gene_name,
                        "gene_id": gene_id,
                        "feature": feature_type,
                    }
            
            return hits
            
        finally:
            if filtered_path != str(annotation_path):
                os.unlink(filtered_path)
    
    def _convert_gtf_to_bed(
        self,
        gtf_path: Path,
        output_path: str,
        feature_type: str,
    ) -> None:
        """Convert GTF to BED format, filtering by feature type."""
        with open(gtf_path) as f_in, open(output_path, 'w') as f_out:
            for line in f_in:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                if parts[2] != feature_type:
                    continue
                
                chrom = parts[0]
                start = int(parts[3]) - 1  # GTF is 1-based, BED is 0-based
                end = int(parts[4])
                strand = parts[6]
                attrs = parts[8]
                
                f_out.write(f"{chrom}\t{start}\t{end}\t{feature_type}\t0\t{strand}\t{attrs}\n")
    
    def _annotate_site(
        self,
        site: OffTargetSite,
        exon_hits: dict,
        gene_hits: dict,
        annotation_path: Path,
    ) -> AnnotatedSite:
        """Determine region type for a single site."""
        
        # Check exon first (highest priority)
        if site.name in exon_hits:
            hit = exon_hits[site.name]
            return AnnotatedSite(
                chrom=site.chrom,
                start=site.start,
                end=site.end,
                name=site.name,
                mismatch_count=site.mismatch_count,
                strand=site.strand,
                region_type=RegionType.EXON,
                gene_name=hit.get("gene_name"),
                gene_id=hit.get("gene_id"),
            )
        
        # Check if in gene (intron)
        if site.name in gene_hits:
            hit = gene_hits[site.name]
            return AnnotatedSite(
                chrom=site.chrom,
                start=site.start,
                end=site.end,
                name=site.name,
                mismatch_count=site.mismatch_count,
                strand=site.strand,
                region_type=RegionType.INTRON,
                gene_name=hit.get("gene_name"),
                gene_id=hit.get("gene_id"),
            )
        
        # Find nearest gene for intergenic
        nearest = find_nearest_gene(
            site.chrom, site.start, site.end, annotation_path
        )
        
        # Check if in promoter region
        if nearest and nearest.get("distance", 999999) < self.promoter_upstream:
            return AnnotatedSite(
                chrom=site.chrom,
                start=site.start,
                end=site.end,
                name=site.name,
                mismatch_count=site.mismatch_count,
                strand=site.strand,
                region_type=RegionType.PROMOTER,
                gene_name=nearest.get("gene_name"),
                gene_id=nearest.get("gene_id"),
                distance_to_gene=nearest.get("distance"),
            )
        
        # Intergenic
        return AnnotatedSite(
            chrom=site.chrom,
            start=site.start,
            end=site.end,
            name=site.name,
            mismatch_count=site.mismatch_count,
            strand=site.strand,
            region_type=RegionType.INTERGENIC,
            gene_name=nearest.get("gene_name") if nearest else None,
            gene_id=nearest.get("gene_id") if nearest else None,
            distance_to_gene=nearest.get("distance") if nearest else None,
        )
    
    def _parse_gene_name(self, attrs: str) -> Optional[str]:
        """Extract gene_name from GTF attributes."""
        for part in attrs.split(';'):
            part = part.strip()
            if 'gene_name' in part:
                # Handle both quoted and unquoted values
                if '"' in part:
                    return part.split('"')[1]
                else:
                    return part.split()[-1]
        return None
    
    def _parse_gene_id(self, attrs: str) -> Optional[str]:
        """Extract gene_id from GTF attributes."""
        for part in attrs.split(';'):
            part = part.strip()
            if 'gene_id' in part and 'gene_name' not in part:
                if '"' in part:
                    return part.split('"')[1]
                else:
                    return part.split()[-1]
        return None
