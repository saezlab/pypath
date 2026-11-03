# Build only the serving packages; resource building belongs to the worker image.
FROM ghcr.io/astral-sh/uv:0.7.6 AS uv
FROM python:3.12-slim-bookworm AS wheels
COPY --from=uv /uv /usr/local/bin/uv
WORKDIR /src
COPY pyproject.toml uv.lock ./
COPY packages ./packages
COPY deploy/build-constraints.txt ./deploy/build-constraints.txt
RUN uv export --frozen --package omnipath-api --no-dev --no-emit-workspace --no-hashes -o /api-requirements.txt \
    && uv build --package omnipath-core --wheel --build-constraints deploy/build-constraints.txt --out-dir /wheels \
    && uv build --package omnipath-api --wheel --build-constraints deploy/build-constraints.txt --out-dir /wheels

FROM python:3.12-slim-bookworm AS api
COPY --from=uv /uv /usr/local/bin/uv
RUN apt-get update \
    && apt-get install -y --no-install-recommends curl libgomp1 \
    && rm -rf /var/lib/apt/lists/*
ENV PYTHONUNBUFFERED=1 OMNIPATH_READ_ONLY=true
WORKDIR /app
COPY --from=wheels /wheels /wheels
COPY --from=wheels /api-requirements.txt /requirements.txt
COPY --from=wheels /src/deploy/build-constraints.txt /build-constraints.txt
RUN uv pip install --system --no-cache --no-deps --build-constraints /build-constraints.txt -r /requirements.txt /wheels/*.whl \
    && rm -rf /wheels /requirements.txt /build-constraints.txt
EXPOSE 8085
HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
    CMD curl -fsS http://127.0.0.1:8085/api/health >/dev/null
CMD ["python", "-m", "omnipath_api.cli", "--host", "0.0.0.0", "--port", "8085", "--data-root", "/data"]
