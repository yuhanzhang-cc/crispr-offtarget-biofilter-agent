"""Bedtools wrapper utilities."""

import subprocess
from pathlib import Path
from typing import Optional


class BedtoolsError(Exception):
    """Error running bedtools."""
    pass


def run_bedtools_intersect(
    a_path: str,
    b_path: str,
    wa: bool = False,
    wb: bool = False,
    loj: bool = False,
    sorted_input: bool = False,
) -> subprocess.CompletedProcess:
    """
    Run bedtools intersect.
    
    Args:
        a_path: Path to file A
        b_path: Path to file B
        wa: Write original A entry
        wb: Write original B entry
        loj: Left outer join (report all A features)
        sorted_input: Input files are sorted (enables faster algorithm)
        
    Returns:
        CompletedProcess with stdout containing intersection results
        
    Raises:
        BedtoolsError: If bedtools fails
    """
    cmd = ["bedtools", "intersect", "-a", a_path, "-b", b_path]
    
    if wa:
        cmd.append("-wa")
    if wb:
        cmd.append("-wb")
    if loj:
        cmd.append("-loj")
    if sorted_input:
        cmd.append("-sorted")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        return result
    except subprocess.CalledProcessError as e:
        raise BedtoolsError(f"bedtools intersect failed: {e.stderr}") from e
    except FileNotFoundError:
        raise BedtoolsError(
            "bedtools not found. Please install bedtools and ensure it's in PATH."
        )


def run_bedtools_closest(
    a_path: str,
    b_path: str,
    d: bool = True,
    k: int = 1,
) -> subprocess.CompletedProcess:
    """
    Run bedtools closest to find nearest features.
    
    Args:
        a_path: Path to file A
        b_path: Path to file B
        d: Report distance
        k: Report k closest hits
        
    Returns:
        CompletedProcess with results
    """
    cmd = ["bedtools", "closest", "-a", a_path, "-b", b_path]
    
    if d:
        cmd.append("-d")
    if k > 1:
        cmd.extend(["-k", str(k)])
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
        return result
    except subprocess.CalledProcessError as e:
        raise BedtoolsError(f"bedtools closest failed: {e.stderr}") from e


def check_bedtools_installed() -> bool:
    """Check if bedtools is installed and available."""
    try:
        subprocess.run(
            ["bedtools", "--version"],
            capture_output=True,
            check=True,
        )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False
