#!/usr/bin/env python

from elasticsearch import Elasticsearch


ES_QUERY = {
  "size": 3000,
  "sort": [
    "dcc_project_code",
    "submitter_donor_id"
  ],
  "fields": [
    "dcc_project_code",
    "submitter_donor_id",
    "wgs.normal_specimen.bwa_alignment.gnos_repo",
    "wgs.normal_specimen.bwa_alignment.gnos_id",
    "wgs.normal_specimen.bwa_alignment.files.file_name",
    "wgs.normal_specimen.bwa_alignment.files.file_size",
    "wgs.normal_specimen.bwa_alignment.files.file_md5sum",
    "wgs.tumor_specimens.bwa_alignment.gnos_repo",
    "wgs.tumor_specimens.bwa_alignment.gnos_id",
    "wgs.tumor_specimens.bwa_alignment.files.file_name",
    "wgs.tumor_specimens.bwa_alignment.files.file_size",
    "wgs.tumor_specimens.bwa_alignment.files.file_md5sum"
  ],
  "filter": {
    "terms": {
      "dcc_project_code": [
        "LICA-FR",
        "BTCA-SG",
        "CMDI-UK",
        "EOPC-DE",
        "OV-AU",
        "BOCA-UK",
        "BRCA-EU",
        "PACA-AU",
        "CLLE-ES"
      ]
    }
  },
  "aggs": {
    "project": {
      "terms": {"field": "dcc_project_code", "size": 100, "order" : { "_count" : "asc" }},
      "aggs": {
        "normal_repo": {
          "terms": {"field": "wgs.normal_specimen.bwa_alignment.gnos_repo"}
        },
        "tumor": {
          "nested": {"path": "wgs.tumor_specimens"},
          "aggs": {
            "tumor_repo": {
              "terms": {"field": "wgs.tumor_specimens.bwa_alignment.gnos_repo"}
            }
          }
        }
      }
    }
  }
}


es = Elasticsearch(['localhost:9200'])

response = es.search(index='pcawg_summary', body=ES_QUERY)

header = [
          "project_code",
          "donor_id",
          "specimen_type",
          "gnos_id",
          "file",
          "size",
          "md5sum",
          "gnos_repo",
         ]

print('\t'.join(header))

for hit in response['hits']['hits']:
    fields = hit['fields']

    project_code = fields["dcc_project_code"]
    donor_id = fields["submitter_donor_id"]
    normal_gnos_repo = ','.join(fields["wgs.normal_specimen.bwa_alignment.gnos_repo"])
    normal_gnos_id = fields["wgs.normal_specimen.bwa_alignment.gnos_id"]
    normal_bam_file = fields["wgs.normal_specimen.bwa_alignment.files.file_name"]
    normal_bam_size = fields["wgs.normal_specimen.bwa_alignment.files.file_size"]
    normal_bam_md5sum = fields["wgs.normal_specimen.bwa_alignment.files.file_md5sum"]
    tumor_gnos_repo = ','.join(fields["wgs.tumor_specimens.bwa_alignment.gnos_repo"])
    tumor_gnos_id = fields["wgs.tumor_specimens.bwa_alignment.gnos_id"]
    tumor_bam_file = fields["wgs.tumor_specimens.bwa_alignment.files.file_name"]
    tumor_bam_size = fields["wgs.tumor_specimens.bwa_alignment.files.file_size"]
    tumor_bam_md5sum = fields["wgs.tumor_specimens.bwa_alignment.files.file_md5sum"]

    for x in range(len(normal_bam_file)):
        print('\t'.join([project_code[0], donor_id[0], 'normal',
                     normal_gnos_id[0], normal_bam_file[x],
                     str(normal_bam_size[x]), normal_bam_md5sum[x],
                     normal_gnos_repo]))

    for x in range(len(tumor_bam_file)):
        print('\t'.join([project_code[0], donor_id[0], 'tumor',
                     tumor_gnos_id[ int(x / 2) ], tumor_bam_file[x],
                     str(tumor_bam_size[x]), tumor_bam_md5sum[x],
                     tumor_gnos_repo]))




