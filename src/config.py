"""Shared configuration for the DOT1L-H3K79 co-dependency pipeline."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
DATA = ROOT / "data"
STATS = ROOT / "stats"
FIGS = ROOT / "figures"

# The seven CRISPR gene-effect profiles used as network seeds.
SEEDS = ["DOT1L", "MLLT10", "MLLT6", "MLLT1", "MLLT3", "KMT2A", "MEN1"]

# Expected DepMap release. Pipeline asserts every input carries this label.
DATASET_LABEL = "CRISPR (DepMap Public 26Q1+Score, Chronos)"

# Node admission: positive correlation with >= MIN_RECURRENCE seeds,
# or >= SINGLE_SEED_R with any one seed.
MIN_RECURRENCE = 2
SINGLE_SEED_R = 0.35
# Threshold at which the recurrence filter reaches an acceptable q-value
# (see stats/README.md); used to flag the high-confidence subset.
HIGH_CONFIDENCE_RECURRENCE = 3

# STRING
# Pin the STRING release. The unversioned host string-db.org round-robins across
# servers that can run different builds: the same query returned 505 edges from
# one and 490 from another on the same day. The versioned host is stable.
STRING_VERSION = "12.0"
STRING_API = "https://version-12-0.string-db.org/api"
STRING_SPECIES = 9606
STRING_MIN_SCORE = 0.4        # network construction
STRING_RESOLVE_1 = 0.7        # first re-resolution of an unresolved component
STRING_RESOLVE_2 = 0.9        # second re-resolution
STRING_INTER_MODULE = 0.7     # edges drawn between modules in the figure
# Identify the client to STRING, as its terms of use request.
CALLER_IDENTITY = "dot1l-h3k79-codependency-pipeline"

RANDOM_SEED = 0
CONSENSUS_RUNS = 400          # label-propagation restarts for module stability
