#!/usr/bin/env python3
"""Fetch explicitly registered tool artifacts and verify their pinned SHA256."""
import argparse
import json
from pathlib import Path
import tempfile
import urllib.request
from decisions import sha256

ROOT=Path(__file__).resolve().parents[1]

def tool_path(name):
    spec=json.loads((ROOT/'config/tools.json').read_text())[name]
    return ROOT/spec['path'],spec


def check_tool(name):
    path,spec=tool_path(name)
    if not path.is_file() or sha256(path)!=spec['sha256']:
        raise ValueError('Missing or mismatched '+name+'; run python scripts/fetch_tools.py --tool '+name)
    return path


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--tool',required=True,choices=['macs2']);a=p.parse_args()
    path,spec=tool_path(a.tool)
    if path.exists():
        print(check_tool(a.tool))
    else:
        path.parent.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.download-',dir=path.parent) as temporary:
            download=Path(temporary)/'artifact'
            with urllib.request.urlopen(spec['url'],timeout=120) as response,download.open('wb') as output:
                while True:
                    block=response.read(1024*1024)
                    if not block:break
                    output.write(block)
            if sha256(download)!=spec['sha256']:raise ValueError('Download SHA256 mismatch; artifact not installed')
            download.rename(path)
        print(check_tool(a.tool))
