"""The subsets of the served release: what it holds and the step files that read them."""

from __future__ import annotations

from typing import Any


class SetsQueries:
    """Subset products copied into the release database (MetSigDB memberships and the
    exclusion metabolites)."""

    def subset_status(self) -> dict[str, dict[str, Any]]:
        """{product: {row_contract, partial, tables, stats, manifest_sha256}} of the release
        in scope; empty for a release built without subsets."""
        return {product: dict(info) for product, info in self._served().subsets.items()}

    def run_subset_query(
        self,
        product: str,
        name: str,
        params: dict[str, Any],
        *,
        row_contract: int,
        also: tuple[str, ...] = (),
    ) -> str | None:
        """The JSON answer of a step file over a subset product of the release in scope.

        ``SubsetUnavailable`` when the release has no such subset (or one in a row contract
        this API does not read), or lacks a product the request also needs (``also``).
        """
        status = self.subset_status()
        for needed in (product, *also):
            if needed not in status:
                raise SubsetUnavailable(
                    f"This release has no {_TITLES.get(needed, needed)} subset, so the request "
                    "cannot be answered. Build the release's subsets (omnipath-subsets build) "
                    "and its serving database again."
                )
        found = status[product]["row_contract"]
        if found != row_contract:
            raise SubsetUnavailable(
                f"The {_TITLES.get(product, product)} subset of this release has row contract "
                f"{found}; this API reads {row_contract}."
            )
        return self.run_query(name, params, None, output="json")


class SubsetUnavailable(LookupError):
    """A request needs a subset the release in scope does not have."""

    status_code = 503


_TITLES = {"metsigdb": "MetSigDB", "exclusion_metabolites": "exclusion metabolite"}
