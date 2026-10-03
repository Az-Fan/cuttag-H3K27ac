#!/usr/bin/env python3
"""Expose ignored/retried failures in Nextflow traces even when the launcher exits zero."""
import argparse
import collections
import csv
import json
from pathlib import Path
from decisions import sha256


def audit(run):
    traces=[run/'trace.tsv'] if (run/'trace.tsv').is_file() else sorted((run/'output/pipeline_info').glob('*trace*'))
    if len(traces)!=1:
        return {'state':'REVIEW_REQUIRED','reason':'Expected exactly one trace; absent or ambiguous','traces':list(map(str,traces))}
    with traces[0].open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    failures=[{key:row.get(key) for key in ('task_id','name','status','exit','hash')} for row in rows if row.get('status') not in ('COMPLETED','CACHED')]
    return {'state':'REVIEW_REQUIRED' if failures or not rows else 'COMPUTATIONAL_PASS','trace':str(traces[0].resolve()),'trace_sha256':sha256(traces[0]),
            'task_count':len(rows),'task_status_counts':dict(collections.Counter(row.get('status') for row in rows)),'nonpassing_tasks':failures,
            'interpretation':'Nextflow launcher success can include ignored process errors. Review failed attempts, final outcomes and missing QC before scientific acceptance.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args()
    with a.out.open('x') as f:json.dump(audit(a.run),f,indent=2);f.write('\n')
