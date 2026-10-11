# Dockerfile for elembio-spatialdata-analysis
#
# Build:
#   docker build -t elembio-spatialdata-analysis .
#
# Build with embedded git commit (recommended for CI / provenance):
#   docker build --build-arg GIT_COMMIT=$(git rev-parse HEAD) \
#       -t elembio-spatialdata-analysis .
#
# Run (library image — import the package or launch a Python / Jupyter session):
#   docker run --gpus all -v /path/to/data:/data elembio-spatialdata-analysis \
#       python -c "import elembio_spatialdata_analysis as e; print(e.__version__)"

FROM nvidia/cuda:13.0.0-runtime-ubuntu24.04

ENV PYTHONUNBUFFERED=1
ENV DEBIAN_FRONTEND=noninteractive

ARG GIT_COMMIT=""
ENV ELEMBIO_SPATIALDATA_ANALYSIS_GIT_COMMIT=$GIT_COMMIT

# System deps: just enough to bootstrap uv and let it manage its own venv.
# The library does all I/O in-process (S3 / zip via s3fs / boto3 / zipfile),
# so no extra system tools are needed.
RUN apt-get update && apt-get install -y --no-install-recommends \
      python3.12 \
      python3.12-venv \
      python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Bootstrap uv. PEP 668 marks the system python as externally-managed on
# Ubuntu 24.04; --break-system-packages is safe here because uv then
# manages its own venv at /app/.venv and nothing else uses the system
# site-packages.
RUN pip install --break-system-packages uv

WORKDIR /app

# Project files. .dockerignore keeps tutorials/, tests/, docs/, ci-utils/
# and host .venv/ out of the build context.
COPY pyproject.toml uv.lock README.md LICENSE ./
COPY src/ ./src/

# Runtime deps, resolved from the lockfile only — no implicit updates.
RUN uv sync --frozen

ENV PATH="/app/.venv/bin:$PATH"

# Smoke check: the package imports in the built environment.
RUN python -c "import elembio_spatialdata_analysis as e; print(e.__package_name__, e.__version__)"

CMD ["/bin/bash"]
