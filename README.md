# elembio-spatialdata-analysis

Analysis tools for Element Biosciences AVITI24 spatial data, built on [SpatialData](https://spatialdata.scverse.org/) and [AnnData](https://anndata.readthedocs.io/).

A run is configured with some combination of assays: the 3′ transcriptome (`ThreePrimeUntargeted`), OPS / specialized panels (`SpecializedTargeted`), and imaging (CellPaint morphology, optionally with IF protein). Which are present depends on how the run was set up. The earlier Transcript and Protein barcoding assays are deprecated and remain supported only for reading older datasets. See [docs/spatialdata_object_overview.md](docs/spatialdata_object_overview.md) for the full object structure.

## Installation

Requires Python 3.12 or later on Windows, macOS, or Linux. See
[docs/installation.md](docs/installation.md) for installing Python and for uv or conda
setups.

<!-- TODO(pypi): replace the git+https URL with `pip install elembio-spatialdata-analysis` once published to PyPI. -->

```bash
python3 -m venv .venv                # Windows: py -m venv .venv
source .venv/bin/activate            # Windows: .venv\Scripts\Activate.ps1
pip install "git+https://github.com/Elembio/elembio-spatialdata-analysis.git"
pip install jupyterlab               # optional, for the tutorials
```

Check it worked:

```bash
python -c "import elembio_spatialdata_analysis as ebsd; print(ebsd.__version__)"
```

See **[docs/installation.md](docs/installation.md)** for step-by-step instructions
on every OS, with uv or conda, updating, and troubleshooting.

## Quick start

```python
from elembio_spatialdata_analysis.load import Loader

loader = Loader("path/to/run.zarr.zip")
print(loader.summary())
adata = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
```

The [tutorials](tutorials/README.md) walk through loading a real dataset. They are not
installed with the package; [docs/installation.md](docs/installation.md#get-the-tutorials)
explains how to get them.

## Package layout

```
src/elembio_spatialdata_analysis/
├── load/
│   ├── backends.py       # RunInfo ABC + concrete backends (local dir / zip / S3 / in-memory)
│   ├── loader.py         # Multi-run Loader (load_tables / load_points / export_tables)
│   └── inputs/           # Typed parsers for metadata/inputs/{bases2fastq,cells2stats,
│                         #   fastq2parquet_3prime,run,witness_hashes.json}
├── pp/                   # Preprocessing: filtering, normalization
├── qc/
│   ├── gates.py          # YAML-driven traffic-light evaluation
│   ├── config/           # qc_gates.yml — default traffic-light thresholds
│   └── aggregate.py      # Per-batch_type and per-well rollups
├── tl/                   # Tools: morphology channel selection, constants
└── vis/
    ├── filter/               # Per-well histograms / violins / imaging feature plots
    └── spatial/              # Well-image plotting (channels, outlines, table fill)
```

## Documentation

- [Installation](docs/installation.md): every OS, pip / uv / conda, troubleshooting.
- [Tutorials](tutorials/README.md): start with `GettingStarted.ipynb`.
- [SpatialData object overview](docs/spatialdata_object_overview.md) and
  [AnnData object overview](docs/anndata_object_overview.md): what's inside a run
  and where each field lives, including cell coordinates.
- Module specifications, installed alongside the code and readable by coding
  agents: [`load`](src/elembio_spatialdata_analysis/load/SPECIFICATION.md),
  [`pp.filter`](src/elembio_spatialdata_analysis/pp/filter/SPECIFICATION.md),
  [`pp.normalize`](src/elembio_spatialdata_analysis/pp/normalize/SPECIFICATION.md),
  [`vis.filter`](src/elembio_spatialdata_analysis/vis/filter/SPECIFICATION.md).

## License

Licensed under the [Apache License 2.0](LICENSE). Copyright © 2026 Element Biosciences, Inc.
