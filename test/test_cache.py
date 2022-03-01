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

def test_simple_cache(cache):
    # specify name space of cache function
    @cache.cache("test_cache")
    def get_array():
        sleep(2)
        return [1, 2, 3, 4, 5]

    starttime = time()
    _ = get_array()
    duration1 = time() - starttime

    starttime = time()
    _ = get_array()
    duration2 = time() - starttime

    cache.invalidate(get_array, "test_cache")

    assert duration1 > 2 and duration2 < 1



def test_multiple_params(cache):
    
    @cache.region("long_term", "test_cache")
    def get_array(l):
        sleep(2)
        array = []
        for i in range(l):
            array.append(i)
        return array

    starttime = time()
    _ = get_array(5)
    duration1 = time() - starttime

    starttime = time()
    _ = get_array(5)
    duration2 = time() - starttime

    starttime = time()
    _ = get_array(6)
    duration3 = time() - starttime

    cache.region_invalidate(get_array, None, "test_cache", 5)
    cache.region_invalidate(get_array, None, "test_cache", 6)
    # invalidate for each param separately; omission does not work


    assert duration1 > 2 and duration2 < 1 and duration3 > 2


def test_clear_entire_namespace(cache):
    """
    Clearing the namespace after the fact is not possible without 
    specifying all params of the call; therefore, a wrapper is needed
    that invalidates the cache before the single result is returned.
    """

    def get_array_invalidate(l, invalidate=False):
        @cache.region("long_term", "test_cache")
        def get_array(l):
            sleep(1)
            array = []
            for i in range(l):
                array.append(i)

            return array

        if invalidate:
            cache.region_invalidate(get_array, None, "test_cache", l)

        return get_array(l)

    
    _ = get_array_invalidate(5)
    _ = get_array_invalidate(6)

    starttime = time()
    _ = get_array_invalidate(5, invalidate = True)
    duration2 = time() - starttime

    starttime = time()
    _ = get_array_invalidate(6, invalidate = True)
    duration3 = time() - starttime

    assert duration2 > 1 and duration3 > 1