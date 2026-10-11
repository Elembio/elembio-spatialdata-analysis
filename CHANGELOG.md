# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [0.18.1]

### Changed
- Tutorials stream the public S3 `.zarr.zip` by default, so you can run them without
  downloading the dataset first.

## [0.18.0] - 2026-10-02

First public release.

- `load`: read Element Biosciences SpatialData objects (local directory, `.zarr.zip`,
  or S3) into AnnData tables and Polars points, for one or many runs, with typed
  parsers for the pipeline QC files.
- `pp`: modality-aware cell and gene QC filtering (3′ transcriptome, OPS, imaging,
  and legacy Transcript / Protein) and count normalization.
- `qc`: run-level pass / warn / fail gates and per-well / per-batch summaries.
- `tl`: imaging morphology helpers and modality-name utilities.
- `vis`: per-well QC plots and well-image / single-cell plotting.
- Tutorials: a getting-started notebook and two case studies on public datasets (an
  optical pooled CRISPR screen and a HeLa tunicamycin stress response).
