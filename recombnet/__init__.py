"""RecombNet: Python 3 recombination screening and visual ARG generation."""

from .core import (
    ARGResult,
    BarcodeResult,
    DistanceMatrix,
    ExtendedBarcodeResult,
    GroupAnalysis,
    H1Interval,
    analyze_group,
    b1_matrix,
    breaking_points,
    compute_barcode,
    compute_extended_barcode,
    dist_list,
    load_fasta,
    load_manifest,
    segregating_sites,
    subsample,
    build_arg_graph,
    render_arg_graph,
)

__all__ = [
    "ARGResult",
    "BarcodeResult",
    "DistanceMatrix",
    "ExtendedBarcodeResult",
    "GroupAnalysis",
    "H1Interval",
    "analyze_group",
    "b1_matrix",
    "breaking_points",
    "build_arg_graph",
    "compute_barcode",
    "compute_extended_barcode",
    "dist_list",
    "load_fasta",
    "load_manifest",
    "render_arg_graph",
    "segregating_sites",
    "subsample",
]

__version__ = "1.0.0"
