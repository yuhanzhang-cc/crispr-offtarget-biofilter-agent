"""File I/O utilities."""

from pathlib import Path
from typing import List, Dict, Any
import csv
import json

from src.models.schema import OffTargetSite, AnnotatedSite, ScoredSite


def load_bed_file(path: Path) -> List[OffTargetSite]:
    """
    Load off-target sites from BED file.
    
    Expected format: chrom start end name mismatch_count strand
    """
    sites = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 6:
                continue
            
            site = OffTargetSite(
                chrom=parts[0],
                start=int(parts[1]),
                end=int(parts[2]),
                name=parts[3],
                mismatch_count=int(parts[4]),
                strand=parts[5],
            )
            sites.append(site)
    
    return sites


def save_csv(items: List[Any], path: Path) -> None:
    """Save list of dataclass objects to CSV."""
    if not items:
        return
    
    # Get field names from first item
    if hasattr(items[0], 'to_dict'):
        dicts = [item.to_dict() for item in items]
    else:
        dicts = [item.__dict__ for item in items]
    
    if not dicts:
        return
    
    # Flatten nested dicts for CSV
    flat_dicts = [_flatten_dict(d) for d in dicts]
    
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=flat_dicts[0].keys())
        writer.writeheader()
        writer.writerows(flat_dicts)


def save_json(data: Dict[str, Any], path: Path) -> None:
    """Save data to JSON file."""
    with open(path, 'w') as f:
        json.dump(data, f, indent=2)


def save_text(text: str, path: Path) -> None:
    """Save text to file."""
    with open(path, 'w') as f:
        f.write(text)


def _flatten_dict(d: Dict, parent_key: str = '', sep: str = '_') -> Dict:
    """Flatten nested dictionary."""
    items = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(_flatten_dict(v, new_key, sep).items())
        else:
            items.append((new_key, v))
    return dict(items)
