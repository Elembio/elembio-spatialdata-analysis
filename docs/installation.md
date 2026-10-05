# Installation

<!-- Source of truth for the docs-site page
informatics/multiomics-python-package/install.
Edit here, then copy to the docs site. -->

`elembio-spatialdata-analysis` is a Python package. You install it into a Python
environment and then use it from scripts or Jupyter notebooks. It works on
Windows, macOS, and Linux.

## Quick version

**Prerequisite:** Python 3.12 or later. See [Requirements](#requirements) to check your
version or install Python.

```bash
python3 -m venv .venv                # Windows: py -m venv .venv
source .venv/bin/activate            # Windows: .venv\Scripts\Activate.ps1
pip install "git+https://github.com/Elembio/elembio-spatialdata-analysis.git"
pip install jupyterlab               # optional, for the tutorials
```

<!-- TODO(pypi): once the package is on PyPI, replace the git+https install URLs in this
file with `pip install elembio-spatialdata-analysis`. -->

The rest of this page covers each step in detail, alternative installers (uv,
conda), getting the tutorials, and troubleshooting.

---

## Requirements

| | |
|---|---|
| **Python** | 3.12 or newer (3.12, 3.13, and 3.14 are tested) |
| **Operating system** | Windows, macOS, or Linux (64-bit) |
| **Disk space** | About 1.5 GB for the environment, plus space for your data |
| **Git** | Needed to install from GitHub ([download Git](https://git-scm.com/downloads)) |

To check your Python version, run `python3 --version` (Windows: `py --version`). To
install Python 3.12 or later, use the installer from
[python.org](https://www.python.org/downloads/), or use [uv](#option-b-uv) or
[conda](#option-c-conda--mamba), which install Python as part of setup.

---

## Install

Always install into a **dedicated environment**, not your system Python. That keeps
this package's dependencies from conflicting with other software. Pick one of the
three options below; they all produce the same result.

### Option A: pip + venv (built into Python)

**macOS / Linux**

```bash
mkdir elembio-analysis && cd elembio-analysis
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install "git+https://github.com/Elembio/elembio-spatialdata-analysis.git"
```

**Windows (PowerShell)**

```powershell
mkdir elembio-analysis; cd elembio-analysis
py -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install "git+https://github.com/Elembio/elembio-spatialdata-analysis.git"
```

The first install downloads roughly 1 GB of dependencies and takes a few minutes.
Each time you open a new terminal, re-activate the environment (the `activate` line)
before using the package.

### Option B: uv

[uv](https://docs.astral.sh/uv/getting-started/installation/) is a fast Python
package manager that can also download the right Python version for you.

```bash
mkdir elembio-analysis && cd elembio-analysis
uv venv --python 3.12
source .venv/bin/activate            # Windows: .venv\Scripts\activate
uv pip install "git+https://github.com/Elembio/elembio-spatialdata-analysis.git"
```

### Option C: conda / mamba

Use conda (or mamba) only to create the environment and provide Python, then install
the package with `pip` inside it. Don't install this package's dependencies with
`conda install`; mixing the two for the same packages causes conflicts.

```bash
conda create -n elembio-analysis python=3.12 git -y
conda activate elembio-analysis
pip install "git+https://github.com/Elembio/elembio-spatialdata-analysis.git"
```

(`git` is included in the environment so the GitHub install works even where Git
isn't installed system-wide.)

### Jupyter (optional, for notebooks and tutorials)

With the environment active:

```bash
pip install jupyterlab               # uv: uv pip install jupyterlab
```

---

## Verify the install

With the environment active:

```bash
python -c "import elembio_spatialdata_analysis as ebsd; print(ebsd.__version__)"
```

This should print a version number such as `0.16.0`. In your own code, the
convention is:

```python
import elembio_spatialdata_analysis as ebsd
from elembio_spatialdata_analysis.load import Loader
```

---

## Get the tutorials

The tutorial notebooks are **not** installed with the package (a pip install only
contains the library code). Download them separately:

1. Open the repository's
   [`tutorials/`](https://github.com/Elembio/elembio-spatialdata-analysis/tree/main/tutorials)
   folder and download `GettingStarted.ipynb` (open the file, then use the download
   button).
2. Put it in your working folder, for example the `elembio-analysis` folder created
   above.
3. Start Jupyter from that folder, with the environment active:

   ```bash
   jupyter lab
   ```

4. Open `GettingStarted.ipynb`. If the notebook asks for a kernel, choose the default
   Python 3 kernel. Because Jupyter was started from the active environment, it uses
   that environment.

The notebook's first section downloads the example dataset into a `data/` folder next
to the notebook. The dataset is large (about 87.5 GB), so check your free disk space
first.

If you prefer, you can clone the whole repository just to get the tutorials. Keep
using the package you installed with pip, rather than installing from the clone:

```bash
git clone --depth 1 https://github.com/Elembio/elembio-spatialdata-analysis.git
```

---

## Update or uninstall

**Update to the latest version:**

```bash
pip install --upgrade --force-reinstall --no-deps "git+https://github.com/Elembio/elembio-spatialdata-analysis.git"
```

(`--force-reinstall` is needed because pip doesn't re-fetch a Git URL whose version
number hasn't changed. `--no-deps` keeps your other packages as they are. Drop it if
an update needs newer dependencies.) With uv:
`uv pip install --reinstall-package elembio-spatialdata-analysis "git+https://github.com/Elembio/elembio-spatialdata-analysis.git"`.

**Uninstall:** `pip uninstall elembio-spatialdata-analysis`, or simply delete the
environment (the `.venv` folder, or `conda env remove -n elembio-analysis`).

---

## Troubleshooting

**`ERROR: Package 'elembio-spatialdata-analysis' requires a different Python` or
`No matching distribution found`**
Your environment's Python is older than 3.12. Recreate the environment with Python
3.12 or newer (see [Requirements](#requirements)).

**`git: command not found` / `Cannot find command 'git'`**
Install [Git](https://git-scm.com/downloads) and open a new terminal, or use the
conda option, which installs Git into the environment.

**SSL / certificate errors (`CERTIFICATE_VERIFY_FAILED`, `invalid peer certificate`)**
Common on corporate networks that inspect HTTPS traffic. Upgrade pip first
(`python -m pip install --upgrade pip`): recent pip versions use your operating
system's certificates automatically. On an older pip that you can't upgrade, add
`--use-feature=truststore`. With uv, add `--native-tls` to `uv venv` and
`uv pip install`. If that doesn't help, ask your IT team for the proxy's certificate
settings.

**Windows: `Activate.ps1 cannot be loaded because running scripts is disabled`**
Allow local scripts for your user once, then activate again:
`Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`.

**Jupyter can't import the package (`ModuleNotFoundError`)**
Jupyter is running from a different environment. Either start `jupyter lab` from a
terminal where this environment is active, or register the environment as a named
kernel and select it in the notebook:

```bash
python -m ipykernel install --user --name elembio-analysis --display-name "Python (elembio-analysis)"
```

---

## Next steps

- [Tutorials](https://github.com/Elembio/elembio-spatialdata-analysis/tree/main/tutorials): load a real run and explore its tables.
- [SpatialData object overview](spatialdata_object_overview.md) and [AnnData object overview](anndata_object_overview.md): what's inside a run and where each field lives.
