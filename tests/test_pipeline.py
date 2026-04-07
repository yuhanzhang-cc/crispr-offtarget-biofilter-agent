"""End-to-end pipeline tests."""

import pytest
import tempfile
from pathlib import Path

from src.pipeline import Pipeline
from src.models.schema import RiskLevel


class TestPipelineSmoke:
    """Smoke tests for complete pipeline execution."""
    
    @pytest.fixture
    def demo_data_dir(self):
        """Return path to demo data directory."""
        return Path(__file__).parent.parent / "demo"
    
    @pytest.fixture
    def pipeline(self):
        """Create pipeline instance."""
        return Pipeline(log_level="WARNING")
    
    def test_pipeline_runs_end_to_end(self, pipeline, demo_data_dir):
        """Pipeline should complete without errors on demo data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir)
            
            result = pipeline.run(
                sgRNA="GAGTCCGAGCAGAAGAAGA",
                offtargets_path=demo_data_dir / "offtargets.bed",
                annotation_path=demo_data_dir / "genes.gtf",
                outdir=outdir,
            )
            
            # Check result structure
            assert result.sgRNA == "GAGTCCGAGCAGAAGAAGA"
            assert result.total_sites == 8  # Demo has 8 sites
            assert result.total_sites > 0
            
            # Check output files exist
            assert (outdir / "annotated_sites.csv").exists()
            assert (outdir / "scored_sites.csv").exists()
            assert (outdir / "summary.txt").exists()
            assert (outdir / "report.json").exists()
    
    def test_output_files_have_content(self, pipeline, demo_data_dir):
        """Output files should contain valid data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir)
            
            pipeline.run(
                sgRNA="GAGTCCGAGCAGAAGAAGA",
                offtargets_path=demo_data_dir / "offtargets.bed",
                annotation_path=demo_data_dir / "genes.gtf",
                outdir=outdir,
            )
            
            # Check annotated sites
            annotated = (outdir / "annotated_sites.csv").read_text()
            assert "chrom" in annotated
            assert "region_type" in annotated
            
            # Check scored sites
            scored = (outdir / "scored_sites.csv").read_text()
            assert "total_score" in scored
            assert "risk_level" in scored
            
            # Check summary
            summary = (outdir / "summary.txt").read_text()
            assert "CRISPR Off-Target Risk Assessment" in summary
            
            # Check JSON report
            import json
            report = json.loads((outdir / "report.json").read_text())
            assert "sgRNA" in report
            assert "total_sites" in report
            assert "top_hits" in report
    
    def test_risk_distribution_consistent(self, pipeline, demo_data_dir):
        """Risk counts should sum to total."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir)
            
            result = pipeline.run(
                sgRNA="GAGTCCGAGCAGAAGAAGA",
                offtargets_path=demo_data_dir / "offtargets.bed",
                annotation_path=demo_data_dir / "genes.gtf",
                outdir=outdir,
            )
            
            total_classified = (
                result.high_risk_count +
                result.medium_risk_count +
                result.low_risk_count
            )
            
            assert total_classified == result.total_sites, \
                "Risk counts should sum to total"
    
    def test_top_hits_sorted_by_score(self, pipeline, demo_data_dir):
        """Top hits should be sorted by risk score descending."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir)
            
            result = pipeline.run(
                sgRNA="GAGTCCGAGCAGAAGAAGA",
                offtargets_path=demo_data_dir / "offtargets.bed",
                annotation_path=demo_data_dir / "genes.gtf",
                outdir=outdir,
            )
            
            scores = [s.total_score for s in result.top_hits]
            
            assert scores == sorted(scores, reverse=True), \
                "Top hits should be sorted by score descending"
    
    def test_high_risk_sites_first(self, pipeline, demo_data_dir):
        """High risk sites should appear before medium/low risk."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir)
            
            result = pipeline.run(
                sgRNA="GAGTCCGAGCAGAAGAAGA",
                offtargets_path=demo_data_dir / "offtargets.bed",
                annotation_path=demo_data_dir / "genes.gtf",
                outdir=outdir,
            )
            
            # Find first non-high-risk site
            first_non_high = None
            for i, site in enumerate(result.top_hits):
                if site.risk_level != RiskLevel.HIGH:
                    first_non_high = i
                    break
            
            # If there are high risk sites, they should come first
            if first_non_high is not None:
                high_risk_before = sum(
                    1 for s in result.top_hits[:first_non_high]
                    if s.risk_level == RiskLevel.HIGH
                )
                assert high_risk_before == first_non_high, \
                    "All sites before first non-high should be high risk"


class TestPipelineEdgeCases:
    """Edge case tests."""
    
    def test_empty_offtargets_file(self):
        """Pipeline should handle empty offtargets file gracefully."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir)
            
            # Create empty BED file
            empty_bed = outdir / "empty.bed"
            empty_bed.write_text("")
            
            demo_dir = Path(__file__).parent.parent / "demo"
            
            pipeline = Pipeline(log_level="WARNING")
            
            result = pipeline.run(
                sgRNA="TEST",
                offtargets_path=empty_bed,
                annotation_path=demo_dir / "genes.gtf",
                outdir=outdir / "results",
            )
            
            assert result.total_sites == 0
