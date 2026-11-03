"""Synthetic subprocess workload for scheduler lifecycle and Linux cgroup tests."""

import json
import os
from pathlib import Path
import subprocess
import sys
import time
from omnipath_core.versioning import RESOURCE_FILES, SERVING_FILES

spec = json.loads(Path(sys.argv[1]).read_text())
source = spec["build"]["source"]
directory = Path(spec["result"]).parent
if source == "child":
    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"])
    (directory / "child.pid").write_text(str(child.pid))
if source == "kernel_oom":
    blocks = []
    while True:
        blocks.append(bytearray(16 * 1024 * 1024))
if source == "cancel":
    time.sleep(120)
with Path(spec["events"]).open("a", buffering=1) as events:
    events.write(json.dumps({"stage": "test", "message": f"{source}: working"}) + "\n")
    time.sleep(0.15)
if source == "fail" or (source == "retry" and directory.name == "attempt-1"):
    Path(spec["result"]).write_text(
        json.dumps(
            {"status": "failed", "error": "synthetic failure", "memory_failure": source == "retry"}
        )
    )
    sys.exit(1)
output = Path(spec["build"]["output_dir"]) / "resources" / source / spec["build"]["version"]
output.mkdir(parents=True)
manifest = {"resource": source, "version": spec["build"]["version"]}
(output / "build_manifest.json").write_text(json.dumps(manifest))
for name in RESOURCE_FILES + SERVING_FILES:
    (output / name).write_bytes(b"fixture")
Path(spec["result"]).write_text(
    json.dumps(
        {
            "status": "success",
            "elapsed_seconds": 0.15,
            "peak_rss_bytes": 64 * 1024 * 1024,
            "allowed_cpus": sorted(os.sched_getaffinity(0))
            if hasattr(os, "sched_getaffinity")
            else [],
        }
    )
)
