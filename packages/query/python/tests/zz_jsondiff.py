"""Print the paths where two JSON values differ."""
def diff(a, b, path="", out=None, limit=30):
    out = [] if out is None else out
    if len(out) >= limit:
        return out
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b), key=str):
            if k not in a: out.append(f"{path}.{k}: missing in old (new {str(b[k])[:80]})")
            elif k not in b: out.append(f"{path}.{k}: missing in new (old {str(a[k])[:80]})")
            else: diff(a[k], b[k], f"{path}.{k}", out, limit)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): out.append(f"{path}: lengths {len(a)} vs {len(b)}")
        for i, (x, y) in enumerate(zip(a, b)): diff(x, y, f"{path}[{i}]", out, limit)
    elif a != b:
        out.append(f"{path}: {str(a)[:100]!r} vs {str(b)[:100]!r}")
    return out
