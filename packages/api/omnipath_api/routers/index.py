from __future__ import annotations
from fastapi.responses import RedirectResponse

from fastapi import APIRouter

router = APIRouter()


@router.get("/", include_in_schema=False, response_model=None)
def root() -> RedirectResponse:
    return RedirectResponse(url="/docs")
