"""Data models for CRISPR off-target analysis."""

from dataclasses import dataclass, field
from typing import Optional, Dict, Any
from enum import Enum


class RegionType(str, Enum):
    """Genomic region classification."""
    EXON = "exon"
    INTRON = "intron"
    PROMOTER = "promoter"
    INTERGENIC = "intergenic"
    UNKNOWN = "unknown"


class RiskLevel(str, Enum):
    """Risk classification levels."""
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


@dataclass
class OffTargetSite:
    """Represents a single off-target site."""
    chrom: str
    start: int
    end: int
    name: str
    mismatch_count: int
    strand: str
    
    def __post_init__(self):
        """Validate coordinates."""
        if self.start >= self.end:
            raise ValueError(f"start ({self.start}) must be < end ({self.end})")


@dataclass
class AnnotatedSite:
    """Off-target site with genomic context annotation."""
    # Original site info
    chrom: str
    start: int
    end: int
    name: str
    mismatch_count: int
    strand: str
    
    # Annotation
    region_type: RegionType
    gene_name: Optional[str] = None
    gene_id: Optional[str] = None
    distance_to_gene: Optional[int] = None  # for intergenic sites
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "name": self.name,
            "mismatch_count": self.mismatch_count,
            "strand": self.strand,
            "region_type": self.region_type.value,
            "gene_name": self.gene_name,
            "gene_id": self.gene_id,
            "distance_to_gene": self.distance_to_gene,
        }


@dataclass
class ScoreComponents:
    """Individual components of risk score."""
    region_score: float
    gene_score: float
    pathway_score: float
    
    # Raw values before weighting
    region_raw: float = field(default=0.0)
    gene_raw: float = field(default=0.0)
    pathway_raw: float = field(default=0.0)


@dataclass
class ScoredSite:
    """Annotated site with computed risk score."""
    # From AnnotatedSite
    chrom: str
    start: int
    end: int
    name: str
    mismatch_count: int
    strand: str
    region_type: RegionType
    gene_name: Optional[str] = None
    gene_id: Optional[str] = None
    distance_to_gene: Optional[int] = None
    
    # Scoring
    total_score: float = 0.0
    components: ScoreComponents = field(default_factory=lambda: ScoreComponents(0.0, 0.0, 0.0))
    risk_level: RiskLevel = RiskLevel.LOW
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "name": self.name,
            "mismatch_count": self.mismatch_count,
            "strand": self.strand,
            "region_type": self.region_type.value,
            "gene_name": self.gene_name,
            "gene_id": self.gene_id,
            "distance_to_gene": self.distance_to_gene,
            "total_score": round(self.total_score, 4),
            "risk_level": self.risk_level.value,
            "score_components": {
                "region": round(self.components.region_score, 4),
                "gene": round(self.components.gene_score, 4),
                "pathway": round(self.components.pathway_score, 4),
            },
        }


@dataclass
class PipelineResult:
    """Complete pipeline output."""
    sgRNA: str
    total_sites: int
    high_risk_count: int
    medium_risk_count: int
    low_risk_count: int
    top_hits: list  # List[ScoredSite]
    summary_text: str
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "sgRNA": self.sgRNA,
            "total_sites": self.total_sites,
            "high_risk_count": self.high_risk_count,
            "medium_risk_count": self.medium_risk_count,
            "low_risk_count": self.low_risk_count,
            "top_hits": [s.to_dict() for s in self.top_hits],
            "summary_text": self.summary_text,
        }
