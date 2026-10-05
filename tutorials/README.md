# Tutorials

Hands-on notebooks for working with Element Biosciences **SpatialData objects (SDOs)**, the `.zarr.zip` stores produced by AVITI24 Cytos cell screening runs.

## Before you start

1. Install the package and JupyterLab into a Python environment. See the
   [installation guide](../docs/installation.md).
2. Get this `tutorials/` folder, keeping its layout (for example, clone the repository).
   Notebooks are not included in the pip install; see
   [Get the tutorials](../docs/installation.md#get-the-tutorials).
3. Download the example datasets you need into **`tutorials/data/`**. Each notebook's
   first section has copy-paste download commands. The datasets are large (see
   [Example datasets](#example-datasets)), so check your free disk space first.
4. With the environment active, run `jupyter lab` from `tutorials/` and open a notebook.

Notebooks read the data through relative paths: `data/…` from `GettingStarted.ipynb`,
and `../data/…` from the notebooks in `case-studies/`. You don't need to unzip the
`.zarr.zip` files.

The notebooks are saved with their outputs, so you can read the results on GitHub
before downloading anything.

## Notebooks

| Notebook | Start here? | What it covers |
|----------|:-----------:|----------------|
| [`GettingStarted.ipynb`](GettingStarted.ipynb) | ✓ | Download and load an SDO, explore the per-cell tables inside it. No analysis. |
| [`case-studies/ops-crispr-tgfbr2-tnfa.ipynb`](case-studies/ops-crispr-tgfbr2-tnfa.ipynb) | | OPS CRISPR screen: TGFBR2 knockout RNA/morphology, then TNFα protein IF. |
| [`case-studies/hela-tunicamycin-upr.ipynb`](case-studies/hela-tunicamycin-upr.ipynb) | | HeLa tunicamycin: 3′ RNA volcano and single-cell UPR score. |

## Example datasets

Public, no AWS credentials required.

**Cytos-MPS_TNFa-CRISPR-Screen** (~111 GB) — used by the OPS CRISPR case study.

```
s3://element-public-data/Cytos/elembio-multiomics/Cytos-MPS_TNFa-CRISPR-Screen/0c47c400756647deb9b3a37f502458de/spatialdata/Cytos-MPS_TNFa-CRISPR-Screen.zarr.zip
```

**Cytos-3Prime_HeLa-Tunicamycin** (~87.5 GB) — used by GettingStarted and the tunicamycin case study.

```
s3://element-public-data/Cytos/elembio-multiomics/Cytos-3Prime_HeLa-Tunicamycin/cb6f3a2ebb6c43acac4ad08b2b2a7779/spatialdata/Cytos-3Prime_HeLa-Tunicamycin.zarr.zip
```

## Reference

For field-by-field detail on the data structures, see the reference docs:

- [`../docs/spatialdata_object_overview.md`](../docs/spatialdata_object_overview.md): the SDO container (images, labels, shapes, points).
- [`../docs/anndata_object_overview.md`](../docs/anndata_object_overview.md): the per-cell AnnData tables.
