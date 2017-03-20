#!/usr/bin/env python
import json
import sys

def dup_finder(jsonstring):
    donor = json.loads(jsonstring)
    donor_id = donor.get('donor_unique_id')
    rna_seq_tumors = donor.get('rna_seq', {}).get('alignment', {}).get('tumor', [])
    for t in rna_seq_tumors:
        top_hat_sp = t.get('tophat', {}).get('submitter_specimen_id', '')
        star_sp = t.get('star', {}).get('submitter_specimen_id', '')
        if not top_hat_sp == star_sp: print "donor: {}, specimen_id differs. top_hat: {}, star: {}".format(donor_id, top_hat_sp, star_sp)

        top_hat_sa = t.get('tophat', {}).get('submitter_sample_id', '')
        star_sa = t.get('star', {}).get('submitter_sample_id', '')
        if not top_hat_sa == star_sa: print "donor: {}, sample_id differs. top_hat: {}, star: {}".format(donor_id, top_hat_sa, star_sa)

for l in sys.stdin:
    dup_finder(l)


