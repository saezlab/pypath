# Migrating from pypath.share.curl to download-manager

## Overview

This guide documents the process of migrating pypath input modules from using `pypath.share.curl` to the new `download-manager` module, based on the experience with `transmir.py`.

## Key Differences

### Old Approach (curl)
```python
import pypath.share.curl as curl

c = curl.Curl(
    url,
    silent=False,
    large=True,
    encoding='iso-8859-1',
)
# Access results via c.result
```

### New Approach (download-manager)
```python
from download_manager import DownloadManager

dm = DownloadManager(pkg='pypath')
desc, item, downloader, path = dm._download(url, dest=True)

if item and item.status == 3:  # Status.READY = 3
    opener = item.open(large=True, encoding='iso-8859-1')
    # Access results via opener.result
```

## Important Concepts

### 1. Download Manager with Cache

The download manager integrates with cache-manager, providing automatic caching:
- `DownloadManager(pkg='pypath')` - Creates a cache in the OS cache directory under 'pypath'
- Files are cached automatically, subsequent requests are served from cache
- Cache handles both newly downloaded and previously cached files uniformly

### 2. The `_download()` Method

Returns a tuple of 4 elements:
1. `desc` - Descriptor object (usually not needed)
2. `item` - CacheItem object (important for accessing the file)
3. `downloader` - Downloader object (not reliable for cached files)
4. `path` - File path in cache

### 3. Status Codes

Cache items have status codes (from `cache_manager._status.Status`):
- `0` - UNINITIALIZED
- `1` - WRITE
- `2` - FAILED
- `3` - READY ✓ (what we check for)
- `-1` - TRASH
- `-2` - DELETED

### 4. Opening Files

The cache item's `open()` method:
- Returns an `Opener` object (not the content directly!)
- Access content via `opener.result`
- Automatically handles decompression (gzip, etc.)
- Parameters:
  - `large=True` - Returns an iterator (memory efficient)
  - `large=False` - Returns all content at once
  - `encoding='...'` - Specifies text encoding

## Migration Steps

### 1. Update Imports

```python
# Remove
import pypath.share.curl as curl

# Add
from download_manager import DownloadManager
```

### 2. Replace curl.Curl with DownloadManager

```python
# Old
c = curl.Curl(url, silent=False, large=True, encoding='iso-8859-1')

# New
dm = DownloadManager(pkg='pypath')
desc, item, downloader, path = dm._download(url, dest=True)
```

### 3. Handle File Access

```python
# Old
for l in c.result:
    process_line(l)

# New
if item and item.status == 3:  # Status.READY
    opener = item.open(large=True, encoding='iso-8859-1')
    if opener and opener.result:
        for l in opener.result:
            process_line(l)
```

## Common Patterns

### Pattern 1: Download to Memory (Buffer)

```python
# If you need content in memory instead of file
buffer = dm.download(url, dest=False)
# Note: This returns raw content (may be compressed)
# For automatic decompression, use file-based approach (dest=True)
```

### Pattern 2: Headers and Other Options

```python
# Headers can be passed but format differs from curl
dm._download(
    url,
    dest=True,
    headers=['Accept-Charset: iso-8859-1'],  # List format
)
```

### Pattern 3: Error Handling

```python
desc, item, downloader, path = dm._download(url, dest=True)

if item and item.status == 3:  # Successfully ready
    # Process file
elif item and item.status == 2:  # Failed
    # Handle failure
```

## Benefits of Migration

1. **Automatic Caching** - Files are cached, subsequent runs are faster
2. **Automatic Decompression** - Gzip and other formats handled transparently
3. **Unified Interface** - Same code handles both new downloads and cached files
4. **Better Error Handling** - Status codes provide clear state information
5. **Memory Efficiency** - Iterator-based approach for large files

## Gotchas and Tips

1. **Don't use `downloader.open()` for cached files** - It returns None when file is from cache
2. **Always check `item.status == 3`** - Ensures file is ready and not corrupted
3. **`opener.result`, not `opener`** - The open() method returns an Opener object
4. **Encoding handling** - The encoding parameter in `item.open()` handles text decoding automatically
5. **Line endings** - Lines from the iterator are already stripped of newlines (unlike curl)

## Example: Complete Migration

### Before (curl)
```python
def some_interactions():
    url = urls.urls['some_db']['url']
    c = curl.Curl(url, silent=False, large=True, encoding='utf-8')
    
    result = []
    for l in c.result:
        l = l.strip().split('\t')
        if len(l) >= 3:
            result.append(process_line(l))
    
    return result
```

### After (download-manager)
```python
def some_interactions():
    url = urls.urls['some_db']['url']
    
    dm = DownloadManager(pkg='pypath')
    desc, item, downloader, path = dm._download(url, dest=True)
    
    result = []
    if item and item.status == 3:  # Status.READY
        opener = item.open(large=True, encoding='utf-8')
        if opener and opener.result:
            for l in opener.result:
                l = l.strip()
                if l:  # Skip empty lines
                    l = l.split('\t')
                    if len(l) >= 3:
                        result.append(process_line(l))
    
    return result
```

## Summary

The migration from curl to download-manager requires understanding that:
1. Use `_download()` with `dest=True` for file-based operations
2. Check `item.status == 3` before accessing the file
3. Use `item.open()` which returns an Opener object
4. Access content via `opener.result`
5. The system handles caching and decompression automatically

This approach provides a cleaner, more maintainable solution with built-in caching and compression handling.