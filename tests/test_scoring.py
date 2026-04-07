"""Tests for risk scoring logic."""

import pytest
from pathlib import Path

from src.models.schema import AnnotatedSite, RegionType, RiskLevel
from src.agents.biological_risk_scorer import BiologicalRiskScorer


class TestRiskScoreOrdering:
    """Test that risk scores follow expected biological priorities."""
    
    @pytest.fixture
    def scorer(self):
        """Create scorer with default config."""
        return BiologicalRiskScorer()
    
    def test_exon_higher_than_intron(self, scorer):
        """Exon hits should score higher than intron hits."""
        exon_site = AnnotatedSite(
            chrom="chr1", start=1000, end=1022,
            name="exon_test", mismatch_count=2, strand="+",
            region_type=RegionType.EXON,
            gene_name="TP53",
        )
        intron_site = AnnotatedSite(
            chrom="chr1", start=2000, end=2022,
            name="intron_test", mismatch_count=2, strand="+",
            region_type=RegionType.INTRON,
            gene_name="TP53",
        )
        
        scored = scorer.score([exon_site, intron_site])
        
        exon_score = next(s for s in scored if s.name == "exon_test").total_score
        intron_score = next(s for s in scored if s.name == "intron_test").total_score
        
        assert exon_score > intron_score, "Exon should score higher than intron"
    
    def test_promoter_higher_than_intergenic(self, scorer):
        """Promoter hits should score higher than intergenic hits."""
        promoter_site = AnnotatedSite(
            chrom="chr1", start=1000, end=1022,
            name="promoter_test", mismatch_count=2, strand="+",
            region_type=RegionType.PROMOTER,
            gene_name="TP53",
        )
        intergenic_site = AnnotatedSite(
            chrom="chr1", start=5000, end=5022,
            name="intergenic_test", mismatch_count=2, strand="+",
            region_type=RegionType.INTERGENIC,
            gene_name=None,
        )
        
        scored = scorer.score([promoter_site, intergenic_site])
        
        promoter_score = next(s for s in scored if s.name == "promoter_test").total_score
        intergenic_score = next(s for s in scored if s.name == "intergenic_test").total_score
        
        assert promoter_score > intergenic_score, "Promoter should score higher than intergenic"
    
    def test_important_gene_boosts_score(self, scorer):
        """Sites in important genes should have higher scores."""
        important_site = AnnotatedSite(
            chrom="chr1", start=1000, end=1022,
            name="important_gene", mismatch_count=2, strand="+",
            region_type=RegionType.EXON,
            gene_name="TP53",  # Very important
        )
        normal_site = AnnotatedSite(
            chrom="chr1", start=2000, end=2022,
            name="normal_gene", mismatch_count=2, strand="+",
            region_type=RegionType.EXON,
            gene_name="GENE_X",  # Not in database
        )
        
        scored = scorer.score([important_site, normal_site])
        
        important_score = next(s for s in scored if s.name == "important_gene").total_score
        normal_score = next(s for s in scored if s.name == "normal_gene").total_score
        
        assert important_score > normal_score, "Important gene should boost score"
    
    def test_region_priority_ordering(self, scorer):
        """Test that region priority follows expected order."""
        sites = [
            AnnotatedSite("chr1", 1000, 1022, "exon", 2, "+", RegionType.EXON, "GENE"),
            AnnotatedSite("chr1", 2000, 2022, "promoter", 2, "+", RegionType.PROMOTER, "GENE"),
            AnnotatedSite("chr1", 3000, 3022, "intron", 2, "+", RegionType.INTRON, "GENE"),
            AnnotatedSite("chr1", 4000, 4022, "intergenic", 2, "+", RegionType.INTERGENIC, None),
        ]
        
        scored = scorer.score(sites)
        
        scores = {s.name: s.total_score for s in scored}
        
        assert scores["exon"] > scores["promoter"], "Exon > Promoter"
        assert scores["promoter"] > scores["intron"], "Promoter > Intron"
        assert scores["intron"] > scores["intergenic"], "Intron > Intergenic"
    
    def test_risk_level_thresholds(self, scorer):
        """Test risk level classification based on thresholds."""
        # High risk: exon of important gene
        high_site = AnnotatedSite(
            chrom="chr1", start=1000, end=1022,
            name="high", mismatch_count=1, strand="+",
            region_type=RegionType.EXON,
            gene_name="TP53",
        )
        
        # Low risk: intergenic, no gene
        low_site = AnnotatedSite(
            chrom="chr1", start=5000, end=5022,
            name="low", mismatch_count=4, strand="+",
            region_type=RegionType.INTERGENIC,
            gene_name=None,
        )
        
        scored = scorer.score([high_site, low_site])
        
        high = next(s for s in scored if s.name == "high")
        low = next(s for s in scored if s.name == "low")
        
        assert high.risk_level == RiskLevel.HIGH, "TP53 exon should be high risk"
        assert low.risk_level == RiskLevel.LOW, "Intergenic should be low risk"


class TestScoreComponents:
    """Test score component breakdown."""
    
    def test_components_sum_to_total(self):
        """Weighted components should sum to total score."""
        scorer = BiologicalRiskScorer()
        
        site = AnnotatedSite(
            chrom="chr1", start=1000, end=1022,
            name="test", mismatch_count=2, strand="+",
            region_type=RegionType.EXON,
            gene_name="TP53",
        )
        
        scored = scorer.score([site])[0]
        
        component_sum = (
            scored.components.region_score +
            scored.components.gene_score +
            scored.components.pathway_score
        )
        
        assert abs(component_sum - scored.total_score) < 0.001, \
            "Components should sum to total"
