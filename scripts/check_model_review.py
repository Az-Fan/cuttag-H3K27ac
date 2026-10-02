#!/usr/bin/env python3
"""Shared by direct R and configured Python entries; no implicit scientific acceptance."""
import argparse,sys
from decisions import require_review,known_text
p=argparse.ArgumentParser();p.add_argument('--review',required=True);p.add_argument('--counts',required=True);p.add_argument('--metadata',required=True);p.add_argument('--numerator',required=True);p.add_argument('--denominator',required=True);p.add_argument('--normalization',required=True,choices=['conventional','spikein']);p.add_argument('--spikein');p.add_argument('--alpha',type=float,required=True);p.add_argument('--lfc',type=float,required=True);p.add_argument('--min-count',type=float,required=True);a=p.parse_args()
inputs={'counts_sha256':a.counts,'metadata_sha256':a.metadata}
needed=['upstream_accepted','peaks_accepted','replicates_confirmed','simple_design_accepted']
if a.normalization=='spikein':
 if not a.spikein:p.error('spike-in table required')
 inputs['spikein_sha256']=a.spikein;needed+=['spikein_counting_accepted','spikein_calibration_accepted']
try:
 review=require_review(a.review,needed,inputs,{'numerator':a.numerator,'denominator':a.denominator,'normalization':a.normalization,'design':'~ condition','alpha':a.alpha,'abs_log2fc':a.lfc,'min_total_count':a.min_count})
 if a.normalization=='spikein':
  sp=review.get('spikein',{})
  if sp.get('equal_amount_confirmed') is not True or not known_text(sp.get('added_at')) or not known_text(sp.get('calibration_scope')):raise ValueError('Spike-in experimental scope/information incomplete')
except (ValueError,OSError,KeyError) as e:print(str(e),file=sys.stderr);sys.exit(2)
