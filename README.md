# CRISPR Off-Target Risk Assessment Agent

A minimal but extensible bioinformatics system for biologically informed CRISPR off-target risk assessment.

## Purpose

Given an sgRNA sequence and candidate off-target genomic coordinates, this system:
1. Annotates each off-target site with genomic context (exon/intron/promoter/intergenic)
2. Computes a transparent, biologically informed risk score
3. Generates an interpretable natural-language report

## Key Design Principles

- **Deterministic scoring**: All numerical risk calculations are performed by transparent Python code, not by LLM
- **LLM for interpretation only**: The language model generates human-readable summaries, never scores
- **Modular architecture**: Clear separation between annotation, scoring, and interpretation
- **CLI-first**: Designed for command-line usage and pipeline integration
- **Offline-capable**: No runtime internet access required

## Project Structure

```
crispr_offtarget_agent/
├── README.md                    # This file
├── requirements.txt             # Python dependencies
├── .env.example                 # Example environment variables
├── config/
│   └── default.yaml            # Scoring weights and parameters
├── data/
│   ├── raw/                    # Input data (off-target BED, annotations)
│   └── processed/              # Intermediate outputs
├── demo/                       # Synthetic demo data
├── src/
│   ├── main.py                 # CLI entry point
│   ├── pipeline.py             # Pipeline orchestration
│   ├── models/
│   │   └── schema.py          # Data models
│   ├── agents/
│   │   ├── genomic_context_annotator.py   # Agent 1: BED annotation
│   │   ├── biological_risk_scorer.py      # Agent 2: Risk scoring
│   │   └── llm_interpreter.py             # Agent 3: Report generation
│   ├── utils/
│   │   ├── io.py              # File I/O utilities
│   │   ├── bedtools.py        # Bedtools wrapper
│   │   ├── genome.py          # Genome utilities
│   │   ├── scoring.py         # Scoring utilities
│   │   └── logging_utils.py   # Logging setup
│   └── db/
│       ├── gene_sets.py       # Gene importance database
│       └── pathway_lookup.py  # Pathway importance database
└── tests/
    ├── test_scoring.py        # Risk score tests
    └── test_pipeline.py       # End-to-end tests
```

## Dependencies

### System Requirements
- Python 3.10+
- bedtools (must be installed and in PATH)

### Python Packages
```bash
pip install -r requirements.txt
```

### Installing bedtools

**Ubuntu/Debian:**
```bash
sudo apt-get install bedtools
```

**macOS:**
```bash
brew install bedtools
```

**Conda:**
```bash
conda install -c bioconda bedtools
```

## Quick Start

### Run with Demo Data

```bash
python -m src.main \
    --sgRNA GAGTCCGAGCAGAAGAAGA \
    --offtargets demo/offtargets.bed \
    --annotation demo/genes.gtf \
    --outdir results/demo_run
```

### Outputs

The pipeline produces:
- `annotated_sites.csv` - Genomic context annotations
- `scored_sites.csv` - Risk scores with component breakdown
- `summary.txt` - Human-readable summary
- `report.json` - Structured report for downstream tools

## Configuration

Edit `config/default.yaml` to customize:
- Scoring weights (region, gene importance, pathway importance)
- Region priority multipliers
- LLM prompt templates

## Placeholders for Future Expansion

The codebase includes intentional extension points:

1. **External API Adapters**: `llm_interpreter.py` supports mock and real LLM clients
2. **Database Connectors**: `db/gene_sets.py` and `db/pathway_lookup.py` use local dictionaries but can be extended to query external databases (Ensembl, KEGG, etc.)
3. **Scoring Algorithms**: `biological_risk_scorer.py` is designed to support alternative scoring methods

## Testing

```bash
# Run all tests
python -m pytest tests/

# Run specific test
python -m pytest tests/test_scoring.py -v
```

## License

MIT
