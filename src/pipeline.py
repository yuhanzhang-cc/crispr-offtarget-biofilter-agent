"""Pipeline orchestration for CRISPR off-target analysis."""

from pathlib import Path
from typing import Optional
import logging

from src.models.schema import OffTargetSite, PipelineResult, RiskLevel
from src.agents.genomic_context_annotator import GenomicContextAnnotator
from src.agents.biological_risk_scorer import BiologicalRiskScorer
from src.agents.llm_interpreter import LLMInterpreter
from src.utils.io import load_bed_file, save_csv, save_json, save_text
from src.utils.logging_utils import setup_logging


class Pipeline:
    """Orchestrates the complete off-target analysis workflow."""
    
    def __init__(
        self,
        config_path: Optional[Path] = None,
        log_level: str = "INFO",
    ):
        """
        Initialize pipeline.
        
        Args:
            config_path: Path to YAML configuration
            log_level: Logging level
        """
        self.config_path = config_path
        self.logger = setup_logging(log_level)
        
        # Initialize agents
        self.annotator = GenomicContextAnnotator()
        self.scorer = BiologicalRiskScorer(config_path)
        self.interpreter = LLMInterpreter()
    
    def run(
        self,
        sgRNA: str,
        offtargets_path: Path,
        annotation_path: Path,
        outdir: Path,
    ) -> PipelineResult:
        """
        Run complete analysis pipeline.
        
        Args:
            sgRNA: sgRNA sequence
            offtargets_path: Path to off-target BED file
            annotation_path: Path to annotation file (GTF/GFF/BED)
            outdir: Output directory
            
        Returns:
            PipelineResult with all outputs
        """
        self.logger.info(f"Starting pipeline for sgRNA: {sgRNA}")
        self.logger.info(f"Off-targets: {offtargets_path}")
        self.logger.info(f"Annotation: {annotation_path}")
        
        # Create output directory
        outdir.mkdir(parents=True, exist_ok=True)
        
        # Step 1: Load off-target sites
        self.logger.info("Loading off-target sites...")
        sites = load_bed_file(offtargets_path)
        self.logger.info(f"Loaded {len(sites)} sites")
        
        # Step 2: Annotate genomic context
        self.logger.info("Annotating genomic context...")
        annotated = self.annotator.annotate(sites, annotation_path)
        self.logger.info("Annotation complete")
        
        # Save annotated sites
        annotated_path = outdir / "annotated_sites.csv"
        save_csv(annotated, annotated_path)
        self.logger.info(f"Saved annotated sites to {annotated_path}")
        
        # Step 3: Score biological risk
        self.logger.info("Computing biological risk scores...")
        scored = self.scorer.score(annotated)
        self.logger.info("Scoring complete")
        
        # Save scored sites
        scored_path = outdir / "scored_sites.csv"
        save_csv(scored, scored_path)
        self.logger.info(f"Saved scored sites to {scored_path}")
        
        # Step 4: Generate summary statistics
        high_risk = sum(1 for s in scored if s.risk_level == RiskLevel.HIGH)
        medium_risk = sum(1 for s in scored if s.risk_level == RiskLevel.MEDIUM)
        low_risk = sum(1 for s in scored if s.risk_level == RiskLevel.LOW)
        
        # Get top hits (all high risk + top medium if needed)
        top_hits = [s for s in scored if s.risk_level == RiskLevel.HIGH]
        if len(top_hits) < 5:
            top_hits.extend([s for s in scored if s.risk_level == RiskLevel.MEDIUM][:5-len(top_hits)])
        
        # Step 5: Generate reports
        self.logger.info("Generating reports...")
        
        # Deterministic summary (always generated)
        summary_text = self.interpreter.generate_summary_text(
            sgRNA=sgRNA,
            total_sites=len(scored),
            high_risk=high_risk,
            medium_risk=medium_risk,
            low_risk=low_risk,
            top_hits=top_hits,
        )
        
        summary_path = outdir / "summary.txt"
        save_text(summary_text, summary_path)
        self.logger.info(f"Saved summary to {summary_path}")
        
        # LLM interpretation (optional, may use mock)
        llm_summary = self.interpreter.interpret(sgRNA, top_hits)
        
        # Build result object
        result = PipelineResult(
            sgRNA=sgRNA,
            total_sites=len(scored),
            high_risk_count=high_risk,
            medium_risk_count=medium_risk,
            low_risk_count=low_risk,
            top_hits=top_hits,
            summary_text=summary_text,
        )
        
        # Save JSON report
        report_path = outdir / "report.json"
        save_json(result.to_dict(), report_path)
        self.logger.info(f"Saved JSON report to {report_path}")
        
        # Also save LLM interpretation separately
        llm_path = outdir / "llm_interpretation.txt"
        save_text(llm_summary, llm_path)
        
        self.logger.info("Pipeline complete!")
        
        return result
