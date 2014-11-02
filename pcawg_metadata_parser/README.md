# Parser code for parsing/analyzing/aggregating GNOS XMLs

## Overview

The script takes a list of GNOS metadata XML file as input, one XML per GNOS
Analysis Object. It extracts most important information out for each AO (ie, a
BAM entry), then create its associated donor entry and/or specimen entry if one
does not yet exist. The parser produces two JSONL files, one is organized at
donor leve, the other at BAM level. These JSON docs are also pushed to
Elasticsearch for easy search/browse later.

## Dependencies

Python 3 (will move to python 2 later) and some python modules. Install them
with pip when needed.

Elasticsearch installed and up running on the same host using port 9200.

Kibana, simply download it from Elasticsearch website, then run an HTTP server
over the Kibana folder to serve the static content.

## Run the parser

```
./parse_gnos_xml.py -d <folder_with_GNOS_xmls> -r <revision_string>
```

In addition to build an ES index name as 'pancan_<revision_string>', two JSONL
files will also be created in the current directory.
