"""Global API settings and environment configuration."""

from __future__ import annotations

import os
from pathlib import Path
from pydantic import BaseModel, Field


class Settings(BaseModel):
    data_root: Path = Field(
        default_factory=lambda: Path(os.getenv("OMNIPATH_DATA_ROOT", "data")).resolve()
    )
    admin_secret: str = Field(
        default_factory=lambda: os.getenv("OMNIPATH_ADMIN_SECRET")
        or os.getenv("OMNIPATH_ADMIN_PASSWORD")
        or ""
    )
    root_path: str = Field(default_factory=lambda: os.getenv("OMNIPATH_ROOT_PATH", ""))
    port: int = Field(default_factory=lambda: int(os.getenv("PORT", "8085")))
    host: str = Field(default_factory=lambda: os.getenv("HOST", "0.0.0.0"))
    # Compute the explorer's default views at startup (serving deployments).
    warm_cache: bool = Field(
        default_factory=lambda: os.getenv("OMNIPATH_WARM_CACHE", "").strip().lower()
        in {"1", "true", "yes"}
    )
    read_only: bool = Field(
        default_factory=lambda: os.getenv("OMNIPATH_READ_ONLY", "true").strip().lower()
        not in {"0", "false", "no"}
    )


settings = Settings()
