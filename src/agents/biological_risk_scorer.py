"""Biological risk scorer with transparent weighted scoring."""

from pathlib import Path
from typing import List, Dict
import yaml

from src.models.schema import AnnotatedSite, ScoredSite, ScoreComponents, RiskLevel, RegionType
from src.db.gene_sets import GeneImportanceDB
from src.db.pathway_lookup import PathwayImportanceDB
from src.utils.scoring import normalize_score


class BiologicalRiskScorer:
    """
    Computes biologically informed risk scores.
    
    Risk = w_region * region_weight + w_gene * gene_importance + w_pathway * pathway_importance
    """
    
    def __init__(self, config_path: Path = None):
        """
        Initialize scorer with configuration.
        
        Args:
            config_path: Path to YAML config file. If None, uses default weights.
        """
        self.config = self._load_config(config_path)
        self.gene_db = GeneImportanceDB()
        self.pathway_db = PathwayImportanceDB()
    
    def _load_config(self, config_path: Path = None) -> Dict:
        """Load scoring configuration."""
        default_config = {
            "weights": {
                "region": 0.4,
                "gene": 0.35,
                "pathway": 0.25,
            },
            "region_priority": {
                "exon": 1.0,
                "promoter": 0.8,
                "intron": 0.4,
                "intergenic": 0.1,
                "unknown": 0.5,
            },
            "thresholds": {
                "high_risk": 0.7,
                "medium_risk": 0.4,
                "low_risk": 0.0,
            },
        }
        
        if config_path and config_path.exists():
            with open(config_path) as f:
                user_config = yaml.safe_load(f)
                if user_config:
                    # Merge user config with defaults
                    for key in default_config:
                        if key in user_config:
                            default_config[key].update(user_config[key])
        
        return default_config
    
    def score(self, sites: List[AnnotatedSite]) -> List[ScoredSite]:
        """
        Compute risk scores for annotated sites.
        
        Args:
            sites: List of annotated off-target sites
            
        Returns:
            List of scored sites with component breakdown
        """
        scored = []
        for site in sites:
            scored_site = self._score_single(site)
            scored.append(scored_site)
        
        # Sort by total score descending
        scored.sort(key=lambda x: x.total_score, reverse=True)
        return scored
    
    def _score_single(self, site: AnnotatedSite) -> ScoredSite:
        """Compute risk score for a single site."""
        
        # Get component scores
        region_score, region_raw = self._compute_region_score(site)
        gene_score, gene_raw = self._compute_gene_score(site)
        pathway_score, pathway_raw = self._compute_pathway_score(site)
        
        # Apply weights
        weights = self.config["weights"]
        weighted_region = weights["region"] * region_score
        weighted_gene = weights["gene"] * gene_score
        weighted_pathway = weights["pathway"] * pathway_score
        
        # Total score
        total = weighted_region + weighted_gene + weighted_pathway
        
        # Determine risk level
        thresholds = self.config["thresholds"]
        if total >= thresholds["high_risk"]:
            risk_level = RiskLevel.HIGH
        elif total >= thresholds["medium_risk"]:
            risk_level = RiskLevel.MEDIUM
        else:
            risk_level = RiskLevel.LOW
        
        components = ScoreComponents(
            region_score=weighted_region,
            gene_score=weighted_gene,
            pathway_score=weighted_pathway,
            region_raw=region_raw,
            gene_raw=gene_raw,
            pathway_raw=pathway_raw,
        )
        
        return ScoredSite(
            chrom=site.chrom,
            start=site.start,
            end=site.end,
            name=site.name,
            mismatch_count=site.mismatch_count,
            strand=site.strand,
            region_type=site.region_type,
            gene_name=site.gene_name,
            gene_id=site.gene_id,
            distance_to_gene=site.distance_to_gene,
            total_score=total,
            components=components,
            risk_level=risk_level,
        )
    
    def _compute_region_score(self, site: AnnotatedSite) -> tuple:
        """
        Compute region type score.
        
        Returns:
            Tuple of (weighted_score, raw_score)
        """
        priority = self.config["region_priority"]
        raw = priority.get(site.region_type.value, priority["unknown"])
        return raw, raw  # Region score is already normalized 0-1
    
    def _compute_gene_score(self, site: AnnotatedSite) -> tuple:
        """
        Compute gene importance score.
        
        Returns:
            Tuple of (normalized_score, raw_score)
        """
        if not site.gene_name:
            return 0.0, 0.0
        
        raw = self.gene_db.get_importance(site.gene_name)
        normalized = normalize_score(raw, min_val=0.0, max_val=1.0)
        return normalized, raw
    
    def _compute_pathway_score(self, site: AnnotatedSite) -> tuple:
        """
        Compute pathway importance score.
        
        Returns:
            Tuple of (normalized_score, raw_score)
        """
        if not site.gene_name:
            return 0.0, 0.0
        
        raw = self.pathway_db.get_importance(site.gene_name)
        normalized = normalize_score(raw, min_val=0.0, max_val=1.0)
        return normalized, raw
    
    def get_score_explanation(self, site: ScoredSite) -> str:
        """Generate human-readable explanation of score components."""
        lines = [
            f"Site {site.name} ({site.chrom}:{site.start}-{site.end})",
            f"  Region: {site.region_type.value} (score: {site.components.region_raw:.2f})",
            f"  Gene: {site.gene_name or 'N/A'} (importance: {site.components.gene_raw:.2f})",
            f"  Pathway: {site.components.pathway_raw:.2f}",
            f"  Total: {site.total_score:.3f} ({site.risk_level.value} risk)",
        ]
        return "\n".join(lines)
