"""Small fixture builds through the production flat-table writer."""

import tempfile
from omnipath_build.silver import SilverExtractor
from omnipath_build.writer import ParquetWriter
from tables_fixture import read_tables


def write_observations(resolver, entities, relations, payloads):
    """The written resource tables (``tables_fixture.Tables``)."""
    with tempfile.TemporaryDirectory() as directory:
        writer = ParquetWriter(directory, library_dir=resolver.library_dir)
        extractor = SilverExtractor("fixture", "test")
        extractor.entities, extractor.relations, extractor.payloads = entities, relations, payloads
        writer.append_observations(extractor, resolver)
        return read_tables(writer.close()["files"])
