"""Current product runner publication boundaries without server credentials."""

from copy import deepcopy

import pytest

from omnipath_subsets import runner
from omnipath_subsets.build import BuildResult
from omnipath_subsets import scientific as pipeline


class Database:
    def __init__(self):
        self.state = {
            "identity": ("release", "digest"),
            "phases": {"base": {}, "derived": {}},
            "metadata": {},
            "contents": {},
            "status": "derived",
        }
        self.owner = None
        self.connections = []

    def connect(self, _):
        connection = Connection(self)
        self.connections.append(connection)
        return connection


class Connection:
    def __init__(self, database):
        self.database = database
        self._pending = None
        self.commits = self.rollbacks = 0
        self.closed = False
        self.statements = []

    @property
    def pending(self):
        if self._pending is None:
            self._pending = deepcopy(self.database.state)
        return self._pending

    def cursor(self):
        return Cursor(self)

    def commit(self):
        if self._pending is not None:
            self.database.state = deepcopy(self._pending)
        self._pending = None
        self.commits += 1

    def rollback(self):
        self._pending = None
        self.rollbacks += 1

    def close(self):
        self.closed = True
        self._pending = None
        if self.database.owner is self:
            self.database.owner = None


class Cursor:
    def __init__(self, connection):
        self.connection = connection
        self.rows = []

    def __enter__(self):
        return self

    def __exit__(self, *_):
        pass

    def execute(self, statement, parameters=None):
        text = str(statement)
        connection = self.connection
        connection.statements.append((text, parameters))
        if "pg_try_advisory_lock" in text:
            locked = connection.database.owner is None
            if locked:
                connection.database.owner = connection
            self.rows = [(locked,)]
        elif "SELECT version,manifest_sha256" in text:
            self.rows = [connection.pending["identity"]]
        elif "SELECT phase FROM" in text:
            self.rows = [(phase,) for phase in connection.pending["phases"]]
        elif "subset_build_metadata" in text and "INSERT INTO" in text:
            product, release, digest, stats = parameters
            connection.pending["metadata"][product] = (release, digest, deepcopy(stats.adapted))
        elif "parquet_phase" in text and "INSERT INTO" in text:
            product, seconds, outcome = parameters
            connection.pending["phases"][product] = (seconds, deepcopy(outcome.adapted))
        elif "UPDATE" in text and "parquet_release" in text:
            product, seconds, products, expected, _ = parameters
            published = len(set(products) & set(connection.pending["phases"]))
            connection.pending["status"] = "complete" if published == expected else product
        else:
            self.rows = []

    def fetchone(self):
        return self.rows[0] if self.rows else None

    def fetchall(self):
        return self.rows


@pytest.fixture
def environment(monkeypatch):
    database = Database()
    monkeypatch.setattr(runner.psycopg2, "connect", database.connect)
    return database


def builder(*, failed_product=None, failure=RuntimeError, change_identity=False):
    def run(connection, product, **kwargs):
        # The imported scientific builder may commit internally. Its production
        # runner suppresses those commits until tables+publication markers agree.
        wrapped = pipeline.CommitDeferredConnection(connection)
        current = connection.pending["contents"].get(product, 0)
        connection.pending["contents"][product] = current + 1
        wrapped.commit()
        assert connection.commits == 0
        if change_identity:
            connection.pending["identity"] = ("changed", "changed-digest")
        if failed_product == product:
            raise failure("deliberate product interruption")
        return {"phase_seconds": 0.5, "result": {"rows": current + 1}}

    return run


def test_checkpoint_callbacks_observe_durable_products_metadata_and_held_importer_lock(
    environment, monkeypatch
):
    database = environment
    monkeypatch.setattr(runner, "run_product", builder())
    observed = []

    def committed(product, result):
        assert database.owner is not None and not database.owner.closed
        assert database.state["contents"][product] == 1
        assert database.state["metadata"][product][:2] == ("release", "digest")
        assert product in database.state["phases"]
        observed.append((product, tuple(result.products)))
        result.products[product]["rows"] = -100

    result = runner.build_subsets(
        "unused",
        "target",
        products=("network_views", "cosmos"),
        checkpoint_products=True,
        on_product_committed=committed,
    )
    assert isinstance(result, BuildResult)
    assert result.release == "release" and result.manifest_sha256 == "digest"
    assert observed == [("network_views", ("network_views",)), ("cosmos", ("network_views", "cosmos"))]
    assert result.products["network_views"]["rows"] == 1
    assert len(database.connections) == 3
    assert database.connections[0].statements[0][1] == ("omnipath:main-parquet:target",)
    assert all(connection.closed for connection in database.connections)
    assert database.owner is None


@pytest.mark.parametrize("failure", [RuntimeError, KeyboardInterrupt])
def test_checkpoint_failure_retains_previous_product_and_rolls_back_only_current(
    environment, monkeypatch, failure
):
    database = environment
    monkeypatch.setattr(runner, "run_product", builder(failed_product="cosmos", failure=failure))
    with pytest.raises(failure, match="deliberate"):
        runner.build_subsets(
            "unused", "target", products=("network_views", "cosmos"), checkpoint_products=True
        )
    assert database.state["contents"] == {"network_views": 1}
    assert set(database.state["metadata"]) == {"network_views"}
    assert "cosmos" not in database.state["phases"]
    assert database.connections[-1].rollbacks == 1
    assert database.owner is None


def test_default_mode_rolls_back_all_selected_products_on_later_failure(environment, monkeypatch):
    database = environment
    monkeypatch.setattr(runner, "run_product", builder(failed_product="cosmos"))
    with pytest.raises(RuntimeError, match="deliberate"):
        runner.build_subsets(
            "unused", "target", products=("network_views", "cosmos"), checkpoint_products=False
        )
    assert database.state["contents"] == {}
    assert database.state["metadata"] == {}
    assert set(database.state["phases"]) == {"base", "derived"}


def test_callback_failure_preserves_current_commit_and_does_not_start_next_product(
    environment, monkeypatch
):
    database = environment
    monkeypatch.setattr(runner, "run_product", builder())

    def failed(*_):
        raise OSError("observer persistence failed")

    with pytest.raises(OSError, match="observer"):
        runner.build_subsets(
            "unused",
            "target",
            products=("network_views", "cosmos"),
            checkpoint_products=True,
            on_product_committed=failed,
        )
    assert database.state["contents"] == {"network_views": 1}
    assert set(database.state["metadata"]) == {"network_views"}
    assert len(database.connections) == 2


def test_selected_rebuild_updates_only_explicit_product_and_keeps_other_checkpoints(
    environment, monkeypatch
):
    database = environment
    database.state["contents"] = {
        "exclusion_metabolites": 5,
        "metsigdb": 10,
        "network_views": 20,
        "cosmos": 30,
    }
    database.state["phases"].update(
        {product: "old phase" for product in database.state["contents"]}
    )
    database.state["metadata"] = {product: "old metadata" for product in database.state["contents"]}
    monkeypatch.setattr(runner, "run_product", builder())
    result = runner.build_subsets(
        "unused", "target", products=("cosmos",), checkpoint_products=True
    )
    assert result.products == {"cosmos": {"rows": 31}}
    assert database.state["contents"] == {
        "exclusion_metabolites": 5,
        "metsigdb": 10,
        "network_views": 20,
        "cosmos": 31,
    }
    assert database.state["metadata"]["metsigdb"] == "old metadata"
    assert database.state["phases"]["network_views"] == "old phase"
    assert database.state["status"] == "complete"


def test_release_is_complete_only_when_loaded_and_built_products_have_published(
    environment, monkeypatch
):
    database = environment
    monkeypatch.setattr(runner, "run_product", builder())
    runner.build_subsets("unused", "target", products=runner.BUILT_PRODUCTS, checkpoint_products=True)
    # The subsets are loaded from Parquet, not built here: until they are, it is not complete.
    assert database.state["status"] == runner.BUILT_PRODUCTS[-1]

    database.state["phases"].update({product: "loaded" for product in runner.LOADED_PRODUCTS})
    runner.build_subsets("unused", "target", products=("cosmos",), checkpoint_products=True)
    assert database.state["status"] == "complete"


def test_a_loaded_product_is_not_built_here(environment, monkeypatch):
    monkeypatch.setattr(runner, "run_product", builder())
    with pytest.raises(ValueError, match="Products"):
        runner.build_subsets("unused", "target", products=("metsigdb",), checkpoint_products=True)
    with pytest.raises(ValueError, match="loaded from the release's Parquet"):
        runner.run_products(
            "unused", "target", owner=None, identity=("release", "digest"), products=("metsigdb",)
        )


def test_exclusion_metabolites_come_before_the_products_that_may_read_them():
    assert runner.PRODUCTS[0] == "exclusion_metabolites"
    assert runner.PRODUCTS == runner.LOADED_PRODUCTS + runner.BUILT_PRODUCTS


def test_release_change_between_transactions_stops_before_next_product(environment, monkeypatch):
    database = environment
    monkeypatch.setattr(runner, "run_product", builder())

    def mutate(*_):
        database.state["identity"] = ("changed", "changed-digest")

    with pytest.raises(ValueError, match="release changed"):
        runner.build_subsets(
            "unused",
            "target",
            products=("network_views", "cosmos"),
            checkpoint_products=True,
            on_product_committed=mutate,
        )
    assert database.state["contents"] == {"network_views": 1}
    assert database.state["metadata"]["network_views"][:2] == ("release", "digest")


def test_product_mutating_release_identity_rolls_back_its_tables_and_metadata(
    environment, monkeypatch
):
    database = environment
    monkeypatch.setattr(runner, "run_product", builder(change_identity=True))
    with pytest.raises(ValueError, match="release changed"):
        runner.build_subsets("unused", "target", products=("cosmos",), checkpoint_products=True)
    assert database.state["identity"] == ("release", "digest")
    assert database.state["contents"] == {}
    assert database.state["metadata"] == {}


def test_missing_shared_checkpoint_rejects_product_before_build(environment, monkeypatch):
    database = environment
    del database.state["phases"]["derived"]
    monkeypatch.setattr(
        runner,
        "run_product",
        lambda *_args, **_kwargs: pytest.fail("Incomplete base must not build products"),
    )
    with pytest.raises(ValueError, match="shared derivations"):
        runner.build_subsets("unused", "target", products=("cosmos",), checkpoint_products=True)


def test_lock_collision_rejects_before_worker_is_created(environment):
    database = environment
    previous_owner = object()
    database.owner = previous_owner
    with pytest.raises(ValueError, match="Another migration owns"):
        runner.build_subsets("unused", "target", products=("cosmos",), checkpoint_products=True)
    assert len(database.connections) == 1
    assert database.owner is previous_owner


def test_public_api_is_the_current_runner_without_layout_discovery():
    from omnipath_subsets import build

    assert build.build_subsets is runner.build_subsets
    assert not hasattr(runner, "is_main_layout")


def test_finish_resume_skips_existing_products_and_emits_events_after_commit(
    environment, monkeypatch
):
    database = environment
    database.state["phases"]["network_views"] = {"old": True}
    database.state["contents"]["network_views"] = 10
    monkeypatch.setattr(runner, "run_product", builder())
    events = []
    owner = database.connect("unused")
    runner.acquire_schema_lock(owner, "target")

    def observer(event, **fields):
        assert database.owner is owner
        if event == "phase_committed":
            assert database.state["contents"][fields["phase"]] == 1
        events.append((event, fields["phase"]))

    try:
        result = runner.run_products(
            "unused",
            "target",
            owner=owner,
            identity=("release", "digest"),
            products=("network_views", "cosmos"),
            resume=True,
            observer=observer,
        )
        assert result.products == {"cosmos": {"rows": 1}}
        assert database.state["contents"] == {"network_views": 10, "cosmos": 1}
        assert events == [("phase_start", "cosmos"), ("phase_committed", "cosmos")]
    finally:
        owner.close()


def test_finish_observer_failure_preserves_commit_and_stops_next_product(environment, monkeypatch):
    database = environment
    monkeypatch.setattr(runner, "run_product", builder())
    owner = database.connect("unused")
    runner.acquire_schema_lock(owner, "target")

    def observer(event, **fields):
        if event == "phase_committed":
            raise OSError("observer persistence failed")

    try:
        with pytest.raises(OSError, match="persistence"):
            runner.run_products(
                "unused",
                "target",
                owner=owner,
                identity=("release", "digest"),
                products=("network_views", "cosmos"),
                resume=True,
                observer=observer,
            )
        assert database.state["contents"] == {"network_views": 1}
        assert "network_views" in database.state["phases"]
        assert database.owner is owner
    finally:
        owner.close()


@pytest.mark.parametrize(
    "options",
    [
        {"products": []},
        {"products": ["unknown"]},
        {"products": ["network_views", "network_views"]},
        {"products": "cosmos"},
        {"products": [None]},
        {"checkpoint_products": 1},
        {"checkpoint_products": None},
        {"checkpoint_products": "true"},
        {"checkpoint_products": True, "on_product_committed": "callback"},
        {"on_product_committed": lambda product, result: None},
    ],
)
def test_invalid_selection_fails_before_connecting(monkeypatch, options):
    def unexpected(*args, **kwargs):
        pytest.fail("Invalid subset selection must not connect")

    monkeypatch.setattr(runner.psycopg2, "connect", unexpected)
    with pytest.raises(ValueError):
        runner.build_subsets("unused", "test_release", **options)


@pytest.mark.parametrize("product", runner.BUILT_PRODUCTS)
@pytest.mark.parametrize("release_rows", [[], [("release", "digest"), ("other", "digest")]])
def test_unpublished_or_ambiguous_release_preserves_existing_products(
    environment, monkeypatch, product, release_rows
):
    database = environment
    database.state["contents"][product] = "published product"
    database.state["metadata"][product] = "published metadata"
    before = deepcopy(database.state)
    execute = Cursor.execute

    def unpublished(self, statement, parameters=None):
        execute(self, statement, parameters)
        if "SELECT version,manifest_sha256" in str(statement):
            self.rows = release_rows

    monkeypatch.setattr(Cursor, "execute", unpublished)
    monkeypatch.setattr(
        runner,
        "run_product",
        lambda *_args, **_kwargs: pytest.fail("Unpublished release must not build products"),
    )
    with pytest.raises(ValueError, match="exactly one loaded release"):
        runner.build_subsets("unused", "target", products=(product,), checkpoint_products=True)
    assert database.state == before
    assert len(database.connections) == 1
    assert database.owner is None
