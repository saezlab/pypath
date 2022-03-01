from time import sleep, time
import pypath.share.beaker as bk
import beaker
import pytest



@pytest.fixture
def cache():
    cache = bk.get_cache_manager()
    yield cache

    # teardown


def test_manager(cache):
    assert isinstance(cache, beaker.cache.CacheManager)

def test_region_long_term(cache):
    # specify region and name of cache function
    @cache.region('long_term', 'test_array')
    def get_array():
        sleep(2)
        return [1, 2, 3, 4, 5]

    starttime = time()
    _ = get_array()
    duration1 = time() - starttime

    starttime = time()
    _ = get_array()
    duration2 = time() - starttime

    cache.region_invalidate(get_array, None, 'test_array')

    assert duration1 > 2 and duration2 < 1

