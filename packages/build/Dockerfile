# Build distributions once; serving and optional workers use the same lockfile.
FROM ghcr.io/astral-sh/uv:0.7.6 AS uv
FROM rust:1.97.1-slim-bookworm AS wheels
COPY --from=uv /uv /usr/local/bin/uv
RUN apt-get update \
    && apt-get install -y --no-install-recommends python3 python3-dev python3-venv build-essential \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY pyproject.toml uv.lock ./
COPY packages ./packages
COPY pypath ./pypath
COPY deploy/build-constraints.txt ./deploy/build-constraints.txt
RUN uv export --frozen --package omnipath-api --extra build --no-dev --no-emit-workspace --no-emit-package pypath-omnipath --no-hashes -o /worker-requirements.txt \
    && uv build --package omnipath-core --wheel --build-constraints deploy/build-constraints.txt --out-dir /wheels \
    && uv build --package omnipath-resolver --wheel --build-constraints deploy/build-constraints.txt --out-dir /wheels \
    && uv build --package omnipath-build --wheel --build-constraints deploy/build-constraints.txt --out-dir /wheels \
    && uv build --package omnipath-api --wheel --build-constraints deploy/build-constraints.txt --out-dir /wheels \
    && uv build pypath --wheel --build-constraints deploy/build-constraints.txt --out-dir /wheels

FROM python:3.12-slim-bookworm AS runtime
COPY --from=uv /uv /usr/local/bin/uv
RUN apt-get update \
    && apt-get install -y --no-install-recommends curl libgomp1 \
    && rm -rf /var/lib/apt/lists/*
ENV PYTHONUNBUFFERED=1
ENV OMNIPATH_READ_ONLY=true
WORKDIR /app

FROM runtime AS worker
ENV OMNIPATH_READ_ONLY=false
ENV PYPATH_DOWNLOAD_DATADIR=/data/pypath-data
RUN apt-get update \
    && apt-get install -y --no-install-recommends git \
    && rm -rf /var/lib/apt/lists/*
COPY --from=wheels /wheels /wheels
COPY --from=wheels /worker-requirements.txt /requirements.txt
COPY --from=wheels /src/deploy/build-constraints.txt /build-constraints.txt
RUN uv pip install --system --no-cache --no-deps --build-constraints /build-constraints.txt -r /requirements.txt /wheels/*.whl \
    && rm -rf /wheels /requirements.txt /build-constraints.txt
CMD ["python", "-m", "omnipath_api.jobs.worker", "--data-root", "/data"]
