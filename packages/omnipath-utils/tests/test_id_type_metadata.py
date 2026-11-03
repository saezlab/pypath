"""Every identifier type carries a CURIE prefix and a resolvable URL template.

The URL template derives from the Bioregistry CURIE prefix unless a type sets an
explicit one, so a type needs no per-type curation to gain a working link.
"""

from omnipath_utils.mapping._id_types import IdTypeRegistry


def _reg() -> IdTypeRegistry:
    return IdTypeRegistry.get()


class TestUrlPattern:
    def test_curie_typed_ids_get_a_url_template(self):
        reg = _reg()
        missing = [
            name
            for name in reg.all_names()
            if reg.curie_prefix(name) and not reg.url_pattern(name)
        ]
        assert not missing, f'no URL template for: {missing}'

    def test_derived_template_uses_the_prefix_and_placeholder(self):
        reg = _reg()
        # A type with a CURIE prefix and no explicit override resolves through
        # Bioregistry, with {$id} as the identifier placeholder.
        assert reg.url_pattern('uniprot') == (
            'https://bioregistry.io/uniprot:{$id}'
        )
        assert '{$id}' in reg.url_pattern('uniprot')

    def test_explicit_url_pattern_wins(self):
        reg = _reg()
        info = reg.info('uniprot')
        try:
            info['url_pattern'] = 'https://example.org/{$id}'
            assert reg.url_pattern('uniprot') == 'https://example.org/{$id}'
        finally:
            info.pop('url_pattern', None)

    def test_no_prefix_no_pattern(self):
        reg = _reg()
        # A hypothetical type without a CURIE prefix and no explicit pattern has
        # no URL template rather than a broken one.
        assert reg.url_pattern('does-not-exist') is None


class TestIdPattern:
    def test_declared_id_pattern_returned(self):
        import re

        reg = _reg()
        pattern = reg.id_pattern('uniprot')
        assert pattern
        assert re.match(pattern, 'P04637')
        assert not re.match(pattern, 'not-an-accession')

    def test_undeclared_id_pattern_is_none(self):
        reg = _reg()
        assert reg.id_pattern('does-not-exist') is None


class TestSample:
    def test_sampled_types_return_curie_and_url(self):
        reg = _reg()
        # Common types across entity kinds carry both a CURIE and a URL template.
        for name in ('uniprot', 'entrez', 'ensg', 'chebi'):
            assert reg.curie_prefix(name), name
            url = reg.url_pattern(name)
            assert url and url.startswith('http'), name
