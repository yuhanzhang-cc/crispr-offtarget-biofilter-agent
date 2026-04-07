"""CLI entry point for CRISPR off-target analysis."""

import argparse
import sys
from pathlib import Path

from src.pipeline import Pipeline


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="CRISPR Off-Target Risk Assessment Agent",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with demo data
  python -m src.main \\
      --sgRNA GAGTCCGAGCAGAAGAAGA \\
      --offtargets demo/offtargets.bed \\
      --annotation demo/genes.gtf \\
      --outdir results/demo

  # Run with custom config
  python -m src.main \\
      --sgRNA GAGTCCGAGCAGAAGAAGA \\
      --offtargets data/offtargets.bed \\
      --annotation data/annotation.gtf \\
      --outdir results/run1 \\
      --config config/custom.yaml
        """,
    )
    
    parser.add_argument(
        "--sgRNA",
        required=True,
        help="sgRNA sequence (e.g., GAGTCCGAGCAGAAGAAGA)",
    )
    
    parser.add_argument(
        "--offtargets",
        required=True,
        type=Path,
        help="Path to off-target sites BED file",
    )
    
    parser.add_argument(
        "--annotation",
        required=True,
        type=Path,
        help="Path to genome annotation file (GTF/GFF/BED)",
    )
    
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory for results",
    )
    
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Path to configuration YAML file (optional)",
    )
    
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level (default: INFO)",
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.offtargets.exists():
        print(f"Error: Off-targets file not found: {args.offtargets}", file=sys.stderr)
        sys.exit(1)
    
    if not args.annotation.exists():
        print(f"Error: Annotation file not found: {args.annotation}", file=sys.stderr)
        sys.exit(1)
    
    # Run pipeline
    pipeline = Pipeline(
        config_path=args.config,
        log_level=args.log_level,
    )
    
    result = pipeline.run(
        sgRNA=args.sgRNA,
        offtargets_path=args.offtargets,
        annotation_path=args.annotation,
        outdir=args.outdir,
    )
    
    # Print summary
    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print("=" * 60)
    print(f"Total sites analyzed: {result.total_sites}")
    print(f"High risk:   {result.high_risk_count}")
    print(f"Medium risk: {result.medium_risk_count}")
    print(f"Low risk:    {result.low_risk_count}")
    print(f"\nResults saved to: {args.outdir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
