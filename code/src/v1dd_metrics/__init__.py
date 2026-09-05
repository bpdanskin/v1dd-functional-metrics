"""Per-ROI functional metrics from V1DD NWB sessions.

One module per analysis family under ``families/``; see docs/index.md.
"""

from .config import DEFAULT_CONFIG, REFERENCE_CONFIG, MetricConfig
from .families.drifting_gratings import drifting_gratings_metrics
from .families.natural_images import natural_images_metrics
from .families.natural_movie import natural_movie_metrics
from .families.receptive_fields import receptive_field_metrics
from .families.roi_quality import roi_summary_metrics
from .families.surround_suppression import surround_suppression_metrics
from .schema import OUTPUT_COLUMNS, to_output_schema

__version__ = "0.1.0"

__all__ = [
    "MetricConfig", "DEFAULT_CONFIG", "REFERENCE_CONFIG",
    "OUTPUT_COLUMNS", "to_output_schema",
    "drifting_gratings_metrics", "natural_images_metrics", "natural_movie_metrics",
    "receptive_field_metrics", "roi_summary_metrics", "surround_suppression_metrics",
]
