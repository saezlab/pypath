"""A page saved in place of a data file is an error; an expected HTML page is not."""

from pypath.share.downloads import _is_error_page


def test_html_page_is_data_when_html_is_expected(tmp_path):
    page = tmp_path / 'slctables.html'
    page.write_text('<!DOCTYPE html><html><title>SLCtables</title></html>')
    assert not _is_error_page(page, '.html')
    assert not _is_error_page(page, 'html')
    assert not _is_error_page(page, None)


def test_html_page_in_place_of_a_table_is_an_error(tmp_path):
    table = tmp_path / 'table.tsv'
    table.write_text('<!DOCTYPE html><html><title>Error</title></html>')
    assert _is_error_page(table, None)
    assert _is_error_page(table, '.tsv')
