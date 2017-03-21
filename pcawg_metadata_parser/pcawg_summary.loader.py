#!/usr/bin/env python


import sys
import json 
from elasticsearch1 import Elasticsearch


def init_es(es_host, es_index):
    es = Elasticsearch([ es_host ])
    es.indices.delete( es_index, ignore=[400, 404] )
    es.indices.create( es_index, ignore=400 )

    # create mappings
    with open('pcawg_summary.mapping.json', 'r') as m:
        es_mapping = m.read()
        es.indices.put_mapping(index=es_index, doc_type='donor', body=es_mapping)

    return es


def main(argv=None):
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    es_host = 'localhost:9200'
    es_index = 'pcawg_summary'
    es = init_es(es_host, es_index)

    with open('pcawg_summary.jsonl', 'r') as t:
        for entity in t:
            doc = json.loads(entity)
            es.index(index=es_index, doc_type='donor', id=doc['donor_unique_id'], \
                body=doc, timeout=90 )


if __name__ == "__main__":
    sys.exit(main())

