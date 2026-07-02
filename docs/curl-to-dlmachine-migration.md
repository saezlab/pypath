# Migrating from pypath.share.curl.Curl to dlmachine

## Overview

The old `pypath.share.curl.Curl` class is being replaced by `dlmachine.DownloadManager` 
for all HTTP downloads. A shared `DownloadManager` instance is available in 
`pypath.share.downloads` as `dm`. For cases needing file opening/extraction, 
use `download_and_open()`.

## Quick Reference

### Basic download (small, text response)

**Before (Curl):**
```python
from pypath.share import curl
c = curl.Curl(url, silent=False)
data = c.result  # string
```

**After (dlmachine):**
```python
from pypath.share.downloads import dm
path = dm.download(url)
with open(path) as f:
    data = f.read()
```

### Download with GET parameters

**Before:**
```python
c = curl.Curl(url, get={'query': 'organism_id:9606', 'format': 'tsv'})
data = c.result
```

**After:**
```python
path = dm.download(url, query={'query': 'organism_id:9606', 'format': 'tsv'})
with open(path) as f:
    data = f.read()
```

### Large file / streaming (line-by-line)

**Before:**
```python
c = curl.Curl(url, large=True, silent=False)
for line in c.result:  # c.result is a file handle
    process(line)
```

**After:**
```python
path = dm.download(url)
with open(path) as f:
    for line in f:
        process(line)
```

### Download with POST data

**Before:**
```python
c = curl.Curl(url, post={'key': 'value'})
data = c.result
```

**After:**
```python
path = dm.download(url, post={'key': 'value'})
with open(path) as f:
    data = f.read()
```

### Custom headers

**Before:**
```python
c = curl.Curl(url, req_headers=['User-Agent: my-agent'])
```

**After:**
```python
path = dm.download(url, headers={'User-Agent': 'my-agent'})
```

Note: Curl takes a list of header strings, dlmachine takes a dict.

### Archive handling (zip, gz, tar.gz)

**Before:**
```python
c = curl.Curl(url, large=True, files_needed=['data.txt'])
# c.result is a dict for zip: {filename: file_handle}
# c.result is a file handle for gz
```

**After:**
```python
from pypath.share.downloads import download_and_open
opener = download_and_open(url, 'archive.zip', 'resource_name', 
                           needed=['data.txt'])
# opener.result is a dict for zip: {filename: file_handle}
```

### Gzip-compressed responses (e.g. UniProt with compressed=true)

**Before:**
```python
c = curl.Curl(url, get={'compressed': 'true'}, compr='gz', large=True)
for line in c.result:
    process(line)
```

**After:**
```python
import gzip
path = dm.download(url, query={'compressed': 'true'})
with gzip.open(path, 'rt') as f:
    for line in f:
        process(line)
```

Note: When the server returns gzip-compressed content (e.g. UniProt with
`compressed=true`), dlmachine saves the raw bytes. Use `gzip.open()` in
text mode (`'rt'`) to decompress on read.

### Timeout

**Before:**
```python
c = curl.Curl(url, timeout=2400)
```

**After:**
```python
# dlmachine does not have a per-request timeout parameter currently.
# Use requests timeout via config or add to descriptor.
path = dm.download(url)  # uses default timeout
```

## Parameter Mapping

| Curl parameter | dlmachine equivalent | Notes |
|---------------|---------------------|-------|
| `url` (positional) | `url` (positional) | Same |
| `get={}` | `query={}` | Renamed |
| `post={}` | `post={}` | Same |
| `req_headers=[]` | `headers={}` | List of strings to dict |
| `large=True` | N/A | dm always downloads to file; open it for streaming |
| `silent=False` | N/A | dm uses logging |
| `slow=True` | N/A | No rate limiting in dm |
| `timeout=N` | N/A | Not yet in dm |
| `cache=True` | Default | dm caches by default |
| `follow=True` | Default | dm follows redirects |
| `compr='gz'` | `gzip.open(path, 'rt')` | Decompress on read |
| `files_needed=[]` | `needed=[]` via `download_and_open` | |
| `encoding='utf-8'` | Specify when opening file | |
| `connect_timeout=N` | N/A | Not yet in dm |

## Response Access

| Curl | dlmachine | Notes |
|------|----------|-------|
| `c.result` (large=False) | `open(path).read()` | String content |
| `c.result` (large=True) | `open(path)` | File handle for iteration |
| `c.status` | N/A | Check return value (None = failed) |
| `c.download_failed` | `path is None` | |

## Common Patterns

### Iterating TSV lines:

```python
# Before
c = curl.Curl(url, large=True, silent=False)
for line in c.result:
    parts = line.strip().split('\\t')

# After
path = dm.download(url)
with open(path) as f:
    for line in f:
        parts = line.strip().split('\\t')
```

### JSON response:

```python
# Before
c = curl.Curl(url)
data = json.loads(c.result)

# After
path = dm.download(url)
with open(path) as f:
    data = json.load(f)
```

### Gzip-compressed TSV (e.g. UniProt stream):

```python
# Before
c = curl.Curl(url, get=params, compr='gz', large=True)
for line in c.result:
    parts = line.strip().split('\\t')

# After
import gzip
path = dm.download(url, query=params)
with gzip.open(path, 'rt') as f:
    for line in f:
        parts = line.strip().split('\\t')
```

## Notes

- `dm` is the shared singleton from `pypath.share.downloads`
- All downloads are cached by `cachedir` -- repeated calls return cached files
- The `download_and_open()` helper handles archive extraction (zip, gz, tar.gz)
- dlmachine uses `requests` backend by default (no pycurl needed)
- When the server returns gzip data (e.g. `compressed=true` in query), the
  file on disk is gzip-compressed; use `gzip.open(path, 'rt')` to read it
- `dm.download()` returns `None` on failure -- always check the return value
