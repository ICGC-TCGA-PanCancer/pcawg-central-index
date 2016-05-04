#!/usr/bin/env python

import sys
import os
import re
import glob
import xmltodict
import json
import yaml
import copy
import logging
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from elasticsearch import Elasticsearch
from collections import OrderedDict
import datetime
import dateutil.parser
from itertools import izip
from distutils.version import LooseVersion
import csv
import shutil
from operator import itemgetter

es_queries = [
# query 0: PCAWGDATA-45_Sanger GNOS entries with study field ends with _test
{ 
    "name": "sanger_vcf_with_study_field_ends_with_test",
    "content":
              {
                 "fields":[
                     "donor_unique_id"
                     ],
                 "query":{
                    "wildcard":{
                            "variant_calling_results.sanger_variant_calling.study": "*_test"
                             }
                      }, 
                  "filter": {
                            "bool": {
                              "must": [
                                {
                                  "type": {
                                    "value": "donor"
                                  }
                                }                   
                              ]
                              # "must_not": [
                              #   {
                              #     "terms": {
                              #       "flags.is_manual_qc_failed": [
                              #         "T"
                              #       ]
                              #     }
                              #   },
                              #   {
                              #     "terms": {
                              #       "flags.is_donor_blacklisted": [
                              #         "T"
                              #       ]
                              #     }
                              #   }
                              # ]
                            }
                      },
                 "size": 10000
                }
},
# query 1: PCAWGDATA-47_donors with mismatch lane count
{
    "name": "specimens_with_mismatch_lane_count",
    "content":
            {
               "fields":[
                     "donor_unique_id"
                     ],  
               "filter":{
                  "bool":{
                     "must": [
                        {
                           "type":{
                              "value":"donor"
                           }
                        }
                      ],
                      "should": [
                        {
                           "terms":{
                              "normal_alignment_status.do_lane_count_and_bam_count_match":[
                                 "F"
                              ]
                           }
                        },
                        {
                           "terms":{
                              "normal_alignment_status.do_lane_counts_in_every_bam_entry_match":[
                                 "F"
                              ]
                           }
                        },
                        {
                          "nested": {
                            "path": "tumor_alignment_status",
                            "filter":{
                              "bool": {
                                "should": [
                                  {
                                     "terms":{
                                        "tumor_alignment_status.do_lane_count_and_bam_count_match":[
                                           "F"
                                        ]
                                     }
                                  },
                                  {
                                     "terms":{
                                        "tumor_alignment_status.do_lane_counts_in_every_bam_entry_match":[
                                           "F"
                                        ]
                                     }
                                  }                         
                                ]
                              }
                            }
                          }
                        }                                   
                      ],
                      "must_not": [
                      {
                        "terms": {
                          "gnos_repos_with_alignment_result":[
                            "https://cghub.ucsc.edu/"
                          ]
                        }
                      }
                     ]          
                    }
                  },
                  "size": 10000
        }
},
# query 2: variant calling missing input
{
      "name": "variant_callings_missing_input",
      "content":{
          "fields":[
                "donor_unique_id"
          ],
          "filter": {
                      "bool": {
                        "must": [
                          {
                            "type": {
                              "value": "donor"
                            }
                          },
                          {
                            "terms":{
                              "flags.is_bam_used_by_variant_calling_missing":[
                                "T"
                              ]
                            }
                          }
                        ],
                        "must_not": [
                          {
                            "terms": {
                              "flags.is_manual_qc_failed": [
                                "T"
                              ]
                            }
                          },
                          {
                            "terms": {
                              "flags.is_donor_blacklisted": [
                                "T"
                              ]
                            }
                          }
                        ]
                      }
                    },
                "size": 10000

      }
},
# query 3: get donors with mismatch duplicate bwa bams
{
      "name": "donors_with_mismatch_duplicate_bwa_bams",
      "content":{
           "fields":[
               "donor_unique_id"
           ],  
           "filter":{
              "bool":{
                 "must":[
                    {
                       "type":{
                          "value":"donor"
                       }
                    },          
                    {
                       "terms":{
                          "duplicated_bwa_alignment_summary.exists_mismatch_bwa_bams":[
                             "T"
                          ]
                       }
                    }                        
                  ],
                  "must_not": [
                  {
                    "terms": {
                      "flags.is_manual_qc_failed": [
                              "T"
                            ]
                          }
                      },
                  {
                    "terms": {
                      "flags.is_donor_blacklisted": [
                              "T"
                            ]
                          }
                      },
                  # {
                  #   "terms":{
                  #     "dcc_project_code": [
                  #             "LIRI-JP"
                  #     ]
                  #   }
                  # }
                 ]
                }
              },
              "size": 10000
      }
},
# query 4: get missing gnos_entry from santa_cruz_freeze 
{
     "name": "missing_gnos_entry_from_santa_cruz_freeze",
     "content":{
         "fields": ["donor_unique_id"],
         "filter":{
             "bool": {
                 "must":[
                    {
                       "type":{
                          "value":"donor"
                       }
                    },          
                    {
                       "terms":{
                          "flags.is_santa_cruz_donor":[
                             "T"
                          ]
                       }
                    }                        
                  ],
                  "must_not": [
                  {
                    "terms": {
                      "flags.is_manual_qc_failed": [
                              "T"
                            ]
                          }
                      },
                  {
                    "terms": {
                      "flags.is_donor_blacklisted": [
                              "T"
                            ]
                          }
                      }
                 ]
             }
         },
         "size": 10000
     }
},

# query 5: get specimens with mismatch effective xml md5sum
{
      "name": "specimens_with_mismatch_effective_xml_md5sum",
      "content":{
           "fields":[
               "donor_unique_id"
           ],  
           "filter":{
              "bool":{
                 "must":[
                    {
                       "type":{
                          "value":"donor"
                       }
                    },          
                    {
                       "terms":{
                          "flags.exists_xml_md5sum_mismatch":[
                             "T"
                          ]
                       }
                    }                        
                  ],
                  "must_not": [
                  {
                    "terms": {
                      "flags.is_manual_qc_failed": [
                              "T"
                            ]
                          }
                      },
                  {
                    "terms": {
                      "flags.is_donor_blacklisted": [
                              "T"
                            ]
                          }
                      }
                  ]
                }
              },
              "size": 10000
      }
},

# query 6: get donors with broad_incomplete_uploads 
{
      "name": "broad_incomplete_uploads",
      "content":{
           "fields":[
               "donor_unique_id"
           ],  
           "filter":{
              "bool":{
                 "must":[
                    {
                       "type":{
                          "value":"donor"
                       }
                    },          
                    {
                       "terms":{
                          "flags.is_broad_variant_calling_performed":[
                             "F"
                          ]
                       }
                    },
                    {
                      "bool": {
                        "should":[
                          {
                            "terms": {
                              "flags.broad.broad_file_subset_exist": [
                                "T"
                              ]
                            }
                          },
                          {
                            "terms": {
                              "flags.broad.muse_file_subset_exist": [
                                "T"
                              ]
                            }
                          },
                          {
                            "terms": {
                              "flags.broad.broad_tar_file_subset_exist": [
                                "T"
                              ]
                            }
                          },
                          {
                            "terms": {
                              "flags.broad.exist_file_subsets_mismatch": [
                                "T"
                              ]
                            }
                          }                        
                        ]
                      }
                    }                        
                  ],
                  "must_not": [
                  {
                    "terms": {
                      "flags.is_manual_qc_failed": [
                              "T"
                            ]
                          }
                      },
                  {
                    "terms": {
                      "flags.is_donor_blacklisted": [
                              "T"
                            ]
                          }
                      }
                 ]
                }
              },
              "size": 10000
      }
},

# query 7: get donors with broad_successful_uploads 
{
      "name": "broad_successful_uploads",
      "content":{
           "fields":[
               "donor_unique_id"
           ],  
           "filter":{
              "bool":{
                 "must":[
                    {
                       "type":{
                          "value":"donor"
                       }
                    },          
                    {
                       "terms":{
                          "flags.is_broad_variant_calling_performed":[
                             "T"
                          ]
                       }
                    }                        
                  ],
                  "must_not": [
                  {
                    "terms": {
                      "flags.is_manual_qc_failed": [
                              "T"
                            ]
                          }
                      },
                  {
                    "terms": {
                      "flags.is_donor_blacklisted": [
                              "T"
                            ]
                          }
                      }
                 ]
                }
              },
              "size": 10000
      }
},

# query 8: get donors exist_specimen_type_mismatch
{
    "name": "aliquots_with_mismatch_specimen_type",
    "content":{
        "fields":["donor_unique_id"],  
        "filter":{
          "bool":{
             "must":[
                {
                   "type":{
                      "value":"donor"
                   }
                }                    
              ],
              "should": [
                            {
                               "terms":{
                                  "normal_alignment_status.exist_specimen_type_mismatch":[
                                     "T"
                                  ]
                               }
                            },
                            {
                              "nested": {
                                "path": "tumor_alignment_status",
                                "filter":{
                                  "bool": {
                                    "should": [
                                      {
                                         "terms":{
                                            "tumor_alignment_status.exist_specimen_type_mismatch":[
                                               "T"
                                            ]
                                         }
                                      }                         
                                    ]
                                  }
                                }
                              }
                            }                                   
                          ],
              "must_not": [
                        {
                          "terms": {
                            "flags.is_manual_qc_failed": [
                                    "T"
                            ]
                          }
                        },
                        {
                          "terms": {
                            "flags.is_donor_blacklisted": [
                              "T"
                            ]
                          }
                        }
                      ]
                  }
          },
          "size": 10000
    }
},

# query 9: get donors exist vcf file prefix mismatch
{
      "name": "donors_exist_vcf_file_prefix_mismatch",
      "content":{
           "fields":[
               "donor_unique_id"
           ],  
           "filter":{
              "bool":{
                 "must":[
                    {
                       "type":{
                          "value":"donor"
                       }
                    },          
                    {
                       "terms":{
                          "flags.exists_vcf_file_prefix_mismatch":[
                             "T"
                          ]
                       }
                    }                        
                  ],
                  "must_not": [
                  {
                    "terms": {
                      "flags.is_manual_qc_failed": [
                              "T"
                            ]
                          }
                      },
                  {
                    "terms": {
                      "flags.is_donor_blacklisted": [
                              "T"
                            ]
                          }
                      }
                 ]
                }
              },
              "size": 10000
      }
},

# query 10: get missing gnos_entry from aug2015_release 
{
     "name": "missing_gnos_entry_from_aug2015_release",
     "content":{
         "fields": ["donor_unique_id"],
         "filter":{
             "bool": {
                 "must":[
                    {
                       "type":{
                          "value":"donor"
                       }
                    },          
                    {
                       "terms":{
                          "flags.is_aug2015_donor":[
                             "T"
                          ]
                       }
                    }                        
                  ],
                  "must_not": [
                  {
                    "terms": {
                      "flags.is_manual_qc_failed": [
                              "T"
                            ]
                          }
                      },
                  {
                    "terms": {
                      "flags.is_donor_blacklisted": [
                              "T"
                            ]
                          }
                      }
                 ]
             }
         },
         "size": 10000
     }
},

#query 11: get donors having fixed or version higher than 1.0.4 DKFZ/EMBL
{
  "name": "dkfz_embl_single_gnos_uploads",
  "content":{
       "fields":["donor_unique_id"],  
       "filter":{
          "bool":{
                "must":[         
                  {
                     "term":{
                        "flags.is_dkfz_embl_variant_calling_performed":[
                           "T"
                        ]
                     }
                  },
                  {
                     "term":{
                        "flags.is_dkfz_variant_calling_performed":[
                           "F"
                        ]
                     }
                  },
                  {
                     "term":{
                        "flags.is_embl_variant_calling_performed":[
                           "F"
                        ]
                     }
                  }                        
                ],
                "must_not": [
                        {
                          "term": {
                            "flags.is_manual_qc_failed": [
                                    "T"
                            ]
                          }
                        },
                        {
                          "term": {
                            "flags.is_donor_blacklisted": [
                              "T"
                            ]
                          }
                        }
                      ]
                  }
          },
          "size": 10000
      }
},

# query 12: ESAD-UK_GNOS_entries 
{
    "name": "ESAD-UK_unaligned_gnos_entries",
    "content": {
        "fields": "donor_unique_id",      
        "filter": {
                    "bool": {
                          "must": [
                            {
                              "type": {
                                "value": "donor"
                              }
                            },
                            {
                              "terms": {
                                "dcc_project_code": [
                                  "ESAD-UK"
                                ]
                              }
                            }, 
                            {
                              "terms": {
                                "flags.is_sanger_variant_calling_performed": [
                                  "F"
                                ]
                              }
                            }                   
                          ],
                          "must_not": [
                            {
                              "terms": {
                                "flags.is_manual_qc_failed": [
                                  "T"
                                ]
                              }
                            },
                            # {
                            #   "terms": {
                            #     "flags.is_donor_blacklisted": [
                            #       "T"
                            #     ]
                            #   }
                            # }
                          ]
                        }
                    },
          "size": 10000
        }
},

# query 13: get missing gnos_entry from oct2015_release 
{
     "name": "missing_gnos_entry_from_oct2015_release",
     "content":{
         "fields": ["donor_unique_id"],
         "filter":{
             "bool": {
                 "must":[
                    {
                       "type":{
                          "value":"donor"
                       }
                    },          
                    {
                       "terms":{
                          "flags.is_oct2015_donor":[
                             "T"
                          ]
                       }
                    }                        
                  ],
                  "must_not": [
                  {
                    "terms": {
                      "flags.is_manual_qc_failed": [
                              "T"
                            ]
                          }
                      },
                  {
                    "terms": {
                      "flags.is_donor_blacklisted": [
                              "T"
                            ]
                          }
                      }
                 ]
             }
         },
         "size": 10000
     }
},
# query 14: get sanger uploads 
{
      "name": "sanger_uploads_svfixed_status",
      "content":{
           "fields":[
               "donor_unique_id"
           ],  
           "filter":{
              "bool":{
                 "must":[
                    {
                       "type":{
                          "value":"donor"
                       }
                    },          
                    {
                       "terms":{
                          "flags.is_sanger_variant_calling_performed":[
                             "T"
                          ]
                       }
                    }                        
                  ],
                  "must_not": [
                  {
                    "terms": {
                      "flags.is_manual_qc_failed": [
                              "T"
                            ]
                          }
                      },
                  {
                    "terms": {
                      "flags.is_donor_blacklisted": [
                              "T"
                            ]
                          }
                      }
                 ]
                }
              },
              "size": 10000
      }
},
# query 15: get broad uploads 
{
      "name": "broad_uploads_snvfixed_status",
      "content":{
           "fields":[
               "donor_unique_id"
           ],  
           "filter":{
              "bool":{
                 "must":[
                    {
                       "type":{
                          "value":"donor"
                       }
                    },          
                    {
                       "terms":{
                          "flags.is_broad_variant_calling_performed":[
                             "T"
                          ]
                       }
                    }                        
                  ],
                  "must_not": [
                  {
                    "terms": {
                      "flags.is_manual_qc_failed": [
                              "T"
                            ]
                          }
                      },
                  {
                    "terms": {
                      "flags.is_donor_blacklisted": [
                              "T"
                            ]
                          }
                      }
                 ]
                }
              },
              "size": 10000
      }
},


# query 16: get missing gnos_entry from mar2016_release 
{
     "name": "missing_gnos_entry_from_mar2016_release",
     "content":{
         "fields": ["donor_unique_id"],
         "filter":{
             "bool": {
                 "must":[
                    {
                       "type":{
                          "value":"donor"
                       }
                    },          
                    {
                       "terms":{
                          "flags.is_mar2016_donor":[
                             "T"
                          ]
                       }
                    }                        
                  ],
                  "must_not": [
                  {
                    "terms": {
                      "flags.is_manual_qc_failed": [
                              "T"
                            ]
                          }
                      },
                  {
                    "terms": {
                      "flags.is_donor_blacklisted": [
                              "T"
                            ]
                          }
                      }
                 ]
             }
         },
         "size": 10000
     }
},

# query 17: get oxog uploads 
{
      "name": "variant_call_entries_with_broad_oxog_filter_applied",
      "content":{
           "fields":[
               "donor_unique_id"
           ],  
           "filter":{
              "bool":{
                 "must":[
                    {
                       "type":{
                          "value":"donor"
                       }
                    },          
                    {
                       "terms":{
                          "flags.is_oxog_variant_calling_performed":[
                             "T"
                          ]
                       }
                    }                        
                  ],
                  "must_not": [
                  {
                    "terms": {
                      "flags.is_manual_qc_failed": [
                              "T"
                            ]
                          }
                      },
                  {
                    "terms": {
                      "flags.is_donor_blacklisted": [
                              "T"
                            ]
                          }
                      }
                 ]
                }
              },
              "size": 10000
      }
},

# query 18: get minibam uploads 
{
      "name": "minibams_extracted_from_variant_regions",
      "content":{
           "fields":[
               "donor_unique_id"
           ],  
           "filter":{
              "bool":{
                 "must":[
                    {
                       "type":{
                          "value":"donor"
                       }
                    },          
                    {
                       "terms":{
                          "flags.is_minibam_variant_calling_performed":[
                             "T"
                          ]
                       }
                    }                        
                  ],
                  "must_not": [
                  {
                    "terms": {
                      "flags.is_manual_qc_failed": [
                              "T"
                            ]
                          }
                      },
                  {
                    "terms": {
                      "flags.is_donor_blacklisted": [
                              "T"
                            ]
                          }
                      }
                 ]
                }
              },
              "size": 10000
      }
},

# query 19: get gender information of pcawg donors
{
      "name": "pcawg_donors_gender_info",
      "content":{
           "fields":[
               "donor_unique_id"
           ],  
           "filter":{
              "bool":{
                 "must":[
                    {
                       "type":{
                          "value":"donor"
                       }
                    },
                    {
                      "terms":{
                        "flags.is_normal_specimen_aligned":[
                          "T"
                        ]
                      }
                    },
                    {
                      "terms":{
                        "flags.are_all_tumor_specimens_aligned":[
                          "T"
                        ]
                      }
                    },
                    # {
                    #    "terms":{
                    #       "flags.is_mar2016_donor":[
                    #          "T"
                    #       ]
                    #    }
                    # }                               
                  ],
                  "must_not": [
                  {
                    "terms": {
                      "flags.is_manual_qc_failed": [
                              "T"
                            ]
                          }
                      },
                  {
                    "terms": {
                      "flags.is_donor_blacklisted": [
                              "T"
                            ]
                          }
                      }
                 ]
                }
              },
              "size": 10000
      }
},

]


def get_donor_json(es, es_index, donor_unique_id):
    es_query_donor = {
        "query": {
            "term": {
                "donor_unique_id": donor_unique_id
            }
        }
    }
    response = es.search(index=es_index, body=es_query_donor)

    es_json = response['hits']['hits'][0]['_source']
 
    return es_json


def get_donors_list(es, es_index, es_queries, q_index):
    response = es.search(index=es_index, body=es_queries[q_index].get('content'))
    
    donors_list = []
    for p in response['hits']['hits']:
      donors_list.append(p.get('fields').get('donor_unique_id')[0])

    return donors_list 

def create_report_info(donor_unique_id, es_json, q_index, annotations):
    report_info_list = []

    report_info = OrderedDict()
    report_info['donor_unique_id'] = donor_unique_id
    report_info['submitter_donor_id'] = es_json['submitter_donor_id']
    report_info['dcc_project_code'] = es_json['dcc_project_code']
    
    # annotations = {}
    if q_index == 0:
        add_report_info_0(report_info, report_info_list, es_json)

    if q_index == 1:
        add_report_info_1(report_info, report_info_list, es_json)

    if q_index == 2:
        add_report_info_2(report_info, report_info_list, es_json)

    if q_index == 3:
        report_info['is_train2_donor'] = es_json.get('flags').get('is_train2_donor')
        report_info['is_sanger_variant_calling_performed'] = es_json.get('flags').get('is_sanger_variant_calling_performed') 
        add_report_info_3(report_info, report_info_list, es_json)

    if q_index == 4:
        flag = 'is_santa_cruz_entry' 
        add_report_info_4_10(report_info, report_info_list, es_json, flag)

    if q_index == 10:
        flag = 'is_aug2015_entry' 
        add_report_info_4_10(report_info, report_info_list, es_json, flag)

    if q_index == 13:
        flag = 'is_oct2015_entry' 
        add_report_info_4_10(report_info, report_info_list, es_json, flag)

    if q_index == 16:
        flag = 'is_mar2016_entry' 
        add_report_info_4_10(report_info, report_info_list, es_json, flag)

    if q_index == 5:
        add_report_info_5(report_info, report_info_list, es_json)

    if q_index in [6, 7]:
        add_report_info_6_7(report_info, report_info_list, es_json)

    if q_index == 8:
        add_report_info_8(report_info, report_info_list, es_json)

    if q_index == 9:
        add_report_info_9(report_info, report_info_list, es_json)

    if q_index == 11:
        add_report_info_14_15(report_info, report_info_list, es_json, 'dkfz_embl')

    if q_index == 12:
        # annotations = read_annotations(annotations, 'esad-uk_reheader_uuid', 'esad-uk_uuids.txt')
        add_report_info_12(report_info, report_info_list, es_json, annotations)

    if q_index == 14:
        add_report_info_14_15(report_info, report_info_list, es_json, 'sanger')

    if q_index == 15:
        add_report_info_14_15(report_info, report_info_list, es_json, 'broad')

    if q_index == 17:
        add_report_info_14_15(report_info, report_info_list, es_json, 'oxog')

    if q_index == 18:
        add_report_info_14_15(report_info, report_info_list, es_json, 'minibam')

    if q_index == 19:
        # annotations = read_annotations(annotations, 'gender', '../pcawg-ega-submission/annotation/donor.all_projects.release20.tsv')
        # annotations = read_annotations(annotations, 'gender_update', '../pcawg-ega-submission/annotation/donor.gender_update.release21.tsv')
        add_report_info_19(report_info, report_info_list, es_json, annotations)

    return report_info_list

def get_mapping(source):
    target = {
        "XX": "female",
        "XY": "male",
        "female": "female",
        "male": "male"
    }
    return target.get(source)

def gen_dict_extract(key, var):
    if hasattr(var,'iteritems'):
        for k, v in var.iteritems():
            if k.startswith(key) and not isinstance(v, dict):
                yield v
            if isinstance(v, dict):
                for result in gen_dict_extract(key, v):
                    yield result
            elif isinstance(v, list):
                for d in v:
                    for result in gen_dict_extract(key, d):
                        yield result


def add_report_info_19(report_info, report_info_list, es_json, annotations):
    report_info['gender'] = set()
    report_info['dcc_gender'] = annotations.get('gender').get(es_json.get('donor_unique_id')) if annotations.get('gender').get(es_json.get('donor_unique_id')) else None
    if report_info.get('dcc_gender'): report_info['gender'].add(report_info['dcc_gender'])
    # only use the dcc clinical gender data if donor in the below list
    if not es_json.get('donor_unique_id') in ['LIRI-JP::RK165', 'LIRI-JP::RK212']: 
        for vcf in ['sanger', 'dkfz_embl', 'broad', 'muse']:
            report_info[vcf+'_gender'] = None
            if es_json.get('variant_calling_results') and es_json.get('variant_calling_results').get(vcf+'_variant_calling'):
                if es_json.get('variant_calling_results').get(vcf+'_variant_calling').get('workflow_details') and \
                    es_json.get('variant_calling_results').get(vcf+'_variant_calling').get('workflow_details').get('variant_qc_metrics'):
                    v = list(gen_dict_extract('gender', es_json.get('variant_calling_results').get(vcf+'_variant_calling').get('workflow_details').get('variant_qc_metrics')))
                    report_info[vcf+'_gender'] = get_mapping(v[0]) if v else None
                    if report_info.get(vcf+'_gender'): report_info['gender'].add(report_info[vcf+'_gender'])
    report_info['exist_gender_discrepancy'] = False if len(report_info['gender']) == 1 else True
    report_info_list.append(copy.deepcopy(report_info))


def add_report_info_14_15(report_info, report_info_list, es_json, workflow_type):
    report_info['tumor_aliquot_ids'] = es_json.get('all_tumor_specimen_aliquots')
    if es_json.get('variant_calling_results') and es_json.get('variant_calling_results').get(workflow_type+'_variant_calling'):
        report_info[workflow_type+'_gnos_repo'] = es_json.get('variant_calling_results').get(workflow_type+'_variant_calling').get('gnos_repo')
        report_info[workflow_type+'_gnos_id'] = es_json.get('variant_calling_results').get(workflow_type+'_variant_calling').get('gnos_id')
        report_info['gnos_last_modified'] = es_json.get('variant_calling_results').get(workflow_type+'_variant_calling').get('gnos_last_modified')
        report_info['vcf_workflow_result_version'] = es_json.get('variant_calling_results').get(workflow_type+'_variant_calling').get('vcf_workflow_result_version') 
        report_info_list.append(copy.deepcopy(report_info)) 
    return report_info_list

def add_report_info_12(report_info, report_info_list, es_json, annotations):
    report_info['is_sanger_variant_calling_performed'] = es_json.get('flags').get('is_sanger_variant_calling_performed')
    report_info['is_dkfz_embl_variant_calling_performed'] = es_json.get('flags').get('is_dkfz_embl_variant_calling_performed')
    report_info['dkfz_embl_variant_calling_gnos_id'] = es_json.get('variant_calling_results').get('dkfz_embl_variant_calling').get('gnos_id') if es_json.get('flags').get('is_dkfz_embl_variant_calling_performed') else None
    report_info['dkfz_variant_calling_gnos_id'] = es_json.get('variant_calling_results').get('dkfz_variant_calling').get('gnos_id') if es_json.get('flags').get('is_dkfz_variant_calling_performed') else None
    report_info['embl_variant_calling_gnos_id'] = es_json.get('variant_calling_results').get('embl_variant_calling').get('gnos_id') if es_json.get('flags').get('is_embl_variant_calling_performed') else None
    report_info['is_broad_variant_calling_performed'] = es_json.get('flags').get('is_broad_variant_calling_performed')
    report_info['broad_variant_calling_gnos_id'] = es_json.get('variant_calling_results').get('broad_variant_calling').get('gnos_id') if es_json.get('flags').get('broad').get('broad_file_subset_exist') else None
    report_info['broad_tar_variant_calling_gnos_id'] = es_json.get('variant_calling_results').get('broad_tar_variant_calling').get('gnos_id') if es_json.get('flags').get('broad').get('broad_tar_file_subset_exist') else None
    report_info['muse_variant_calling_gnos_id'] = es_json.get('variant_calling_results').get('muse_variant_calling').get('gnos_id') if es_json.get('flags').get('broad').get('muse_file_subset_exist') else None
    

    if es_json.get('normal_alignment_status'):
        add_report_info_12_aliquot(es_json.get('normal_alignment_status'), report_info, report_info_list, annotations)

    if es_json.get('tumor_alignment_status'):
        for aliquot in es_json.get('tumor_alignment_status'):
            add_report_info_12_aliquot(aliquot, report_info, report_info_list, annotations)

    return report_info_list


def add_report_info_12_aliquot(aliquot, report_info, report_info_list, annotations):
    report_info['library_strategy'] = 'WGS'
    report_info['aliquot_id'] = aliquot.get('aliquot_id')
    report_info['submitter_specimen_id'] = aliquot.get('submitter_specimen_id')
    report_info['submitter_sample_id'] = aliquot.get('submitter_sample_id')
    report_info['dcc_specimen_type'] = aliquot.get('dcc_specimen_type')
    report_info['aligned'] = aliquot.get('aligned')
    report_info['total_lanes'] = aliquot.get('lane_count')

    report_info['aligned_bam_gnos_id'] = None
    report_info['aligned_bam_gnos_repo'] = []

    if aliquot.get('aligned_bam'):
        report_info['aligned_bam_gnos_id'] = aliquot.get('aligned_bam').get('gnos_id')
        report_info['aligned_bam_gnos_repo'] = aliquot.get('aligned_bam').get('gnos_repo')

    report_info['bam_with_unmappable_reads_gnos_id'] = None
    report_info['bam_with_unmappable_reads_gnos_repo'] = []
    if aliquot.get('bam_with_unmappable_reads'):
        report_info['bam_with_unmappable_reads_gnos_id'] = aliquot.get('bam_with_unmappable_reads').get('gnos_id')
        report_info['bam_with_unmappable_reads_gnos_repo'] = aliquot.get('bam_with_unmappable_reads').get('gnos_repo')  

    if aliquot.get('unaligned_bams'):
        report_info['entity_type'] = 'unaligned_bams'
        for unaligned_bams in aliquot.get('unaligned_bams'):
            report_info['unaligned_bams_gnos_id'] = unaligned_bams.get('gnos_id')
            report_info['is_reheader'] = True if unaligned_bams.get('gnos_id') in annotations.get('esad-uk_reheader_uuid') else False
            for gnos_repo in unaligned_bams.get('gnos_repo'):
                report_info['unaligned_bams_gnos_repo'] = gnos_repo
                report_info['unaligned_bams_gnos_metadata_url'] = gnos_repo + 'cghub/metadata/analysisFull/' + report_info['unaligned_bams_gnos_id']
                report_info_list.append(copy.deepcopy(report_info))  

    return report_info_list


# def add_report_info_11(report_info, report_info_list, es_json):
#     report_info['tumor_aliquot_ids'] = es_json.get('all_tumor_specimen_aliquots')
#     if es_json.get('variant_calling_results') and es_json.get('variant_calling_results').get('dkfz_embl_variant_calling'):
#         report_info['gnos_repo'] = es_json.get('variant_calling_results').get('dkfz_embl_variant_calling').get('gnos_repo')
#         report_info['gnos_id'] = es_json.get('variant_calling_results').get('dkfz_embl_variant_calling').get('gnos_id')
#         report_info['gnos_last_modified'] = es_json.get('variant_calling_results').get('dkfz_embl_variant_calling').get('gnos_last_modified')
#         report_info_list.append(copy.deepcopy(report_info))
#     return report_info_list

def add_report_info_9(report_info, report_info_list, es_json):
    report_info['tumor_aliquot_ids'] = es_json.get('all_tumor_specimen_aliquots')
    if es_json.get('variant_calling_results'):
        vcf = es_json.get('variant_calling_results')
        for workflow in ['sanger', 'embl', 'dkfz', 'dkfz_embl', 'broad', 'muse', 'broad_tar']:
            if vcf.get(workflow+'_variant_calling') and vcf.get(workflow+'_variant_calling').get('exists_' + workflow + '_file_prefix_mismatch'):
                report_info['workflow_name'] = workflow+'_variant_calling'
                report_info['gnos_repo'] = vcf.get(workflow+'_variant_calling').get('gnos_repo')[0]
                report_info['gnos_id'] = vcf.get(workflow+'_variant_calling').get('gnos_id')
                file_prefix = set()
                for f in vcf.get(workflow+'_variant_calling').get('files'):
                    file_prefix.add(f.get('file_name').split('.')[0])
                report_info['file_prefix'] = file_prefix
                report_info_list.append(copy.deepcopy(report_info))

    return report_info_list


def add_report_info_8_aliquot(aliquot, report_info, report_info_list):
    if aliquot.get('exist_specimen_type_mismatch'):
        report_info['aliquot_id'] = aliquot.get('aliquot_id')
        report_info['submitter_specimen_id'] = aliquot.get('submitter_specimen_id')
        report_info['submitter_sample_id'] = aliquot.get('submitter_sample_id')
        report_info['dcc_specimen_type'] = aliquot.get('dcc_specimen_type')   
        report_info['aligned'] = True if aliquot.get('aligned') else False
        report_info['exist_aligned_bam_specimen_type_mismatch'] = aliquot.get('exist_aligned_bam_specimen_type_mismatch')
        report_info['exist_unaligned_bam_specimen_type_mismatch'] = aliquot.get('exist_unaligned_bam_specimen_type_mismatch')
        report_info['exist_bam_with_unmappable_reads_specimen_type_mismatch'] = aliquot.get('exist_bam_with_unmappable_reads_specimen_type_mismatch')
        report_info_list.append(copy.deepcopy(report_info))

    return report_info_list


def add_report_info_8(report_info, report_info_list, es_json):
    if es_json.get('normal_alignment_status'):
        add_report_info_8_aliquot(es_json.get('normal_alignment_status'), report_info, report_info_list)

    if es_json.get('tumor_alignment_status'):
        for aliquot in es_json.get('tumor_alignment_status'):
            add_report_info_8_aliquot(aliquot, report_info, report_info_list)

    return report_info_list


def add_report_info_6_7(report_info, report_info_list, es_json):
    report_info['tumor_aliquot_ids'] = es_json.get('aligned_tumor_specimen_aliquots')
    if es_json.get('variant_calling_results'):
        vcf = es_json.get('variant_calling_results')
        for workflow_type in ['broad', 'broad_tar', 'muse']:
            report_info[workflow_type+'_gnos_repo'] = vcf.get(workflow_type+'_variant_calling').get('gnos_repo')[0] if vcf.get(workflow_type+'_variant_calling') else None
            report_info[workflow_type+'_gnos_id'] = vcf.get(workflow_type+'_variant_calling').get('gnos_id') if vcf.get(workflow_type+'_variant_calling') else None
            report_info[workflow_type+'_gnos_last_modified'] = vcf.get(workflow_type+'_variant_calling').get('gnos_last_modified')[0].split('T')[0] if vcf.get(workflow_type+'_variant_calling') else None
        report_info['is_cross_referencing_mismatch'] = es_json.get('flags').get('broad').get('exist_file_subsets_mismatch')

        report_info_list.append(copy.deepcopy(report_info))

    return report_info_list



def add_report_info_5(report_info, report_info_list, es_json):
    if es_json.get('normal_alignment_status'):
        aliquot = es_json.get('normal_alignment_status')
        add_report_info_5_aliquot(aliquot, report_info, report_info_list, 'wgs_bwa_alignment')

    if es_json.get('tumor_alignment_status'):
        for aliquot in es_json.get('tumor_alignment_status'):
            add_report_info_5_aliquot(aliquot, report_info, report_info_list, 'wgs_bwa_alignment')

    if es_json.get('variant_calling_results'):
        for k, v in es_json.get('variant_calling_results').iteritems():
            add_report_info_5_aliquot(v, report_info, report_info_list, k)

    if es_json.get('rna_seq').get('alignment'):
        for k, v in es_json.get('rna_seq').get('alignment').iteritems():
            if not v: continue
            if k=='normal':
                for key, value in v.iteritems():
                    workflow = 'rna_seq_'+key+'_alignment'
                    add_report_info_5_aliquot(value, report_info, report_info_list, workflow)
            if k=='tumor':
                for aliquot in v:
                    for key, value in aliquot.iteritems():
                        workflow = 'rna_seq_'+key+'_alignment'
                        add_report_info_5_aliquot(value, report_info, report_info_list, workflow)                    
    return report_info_list


def add_report_info_5_aliquot(aliquot, report_info, report_info_list, workflow):
    if aliquot.get('exists_xml_md5sum_mismatch'):
        report_info['aliquot_id'] = aliquot.get('aliquot_id') if workflow.endswith('alignment') else None
        report_info['dcc_specimen_type'] = aliquot.get('dcc_specimen_type') if workflow.endswith('alignment') else None
        report_info['workflow'] = workflow
        report_info['gnos_repo'] = aliquot.get('aligned_bam').get('gnos_repo') if workflow.endswith('alignment') else aliquot.get('gnos_repo')
        report_info['gnos_id'] = aliquot.get('aligned_bam').get('gnos_id') if workflow.endswith('alignment') else aliquot.get('gnos_id')
        report_info['effective_xml_md5sum'] = aliquot.get('aligned_bam').get('effective_xml_md5sum') if workflow.endswith('alignment') else aliquot.get('effective_xml_md5sum')
        report_info_list.append(copy.deepcopy(report_info))
    return report_info_list


def add_report_info_4_10(report_info, report_info_list, es_json, flag):
    if es_json.get('bam_files'):
        for bam in es_json.get('bam_files'):
            if not bam.get(flag): 
                continue
            report_info['gnos_id'] =  bam.get('bam_gnos_ao_id')
            report_info_list.append(copy.deepcopy(report_info))

    if es_json.get('variant_calling_results'):
        vcf = es_json.get('variant_calling_results')
        for workflow in ['sanger', 'dkfz_embl', 'broad', 'muse', 'broad_tar']:
            if vcf.get(workflow+'_variant_calling') and vcf.get(workflow+'_variant_calling').get(flag):
                report_info['gnos_id'] =  vcf.get(workflow+'_variant_calling').get('gnos_id')
                report_info_list.append(copy.deepcopy(report_info))

    return report_info_list 


def add_report_info_0(report_info, report_info_list, es_json):
    report_info['gnos_id'] = es_json.get('variant_calling_results').get('sanger_variant_calling').get('gnos_id')
    report_info['study'] = es_json.get('variant_calling_results').get('sanger_variant_calling').get('study')
    for gnos_repo in es_json.get('variant_calling_results').get('sanger_variant_calling').get('gnos_repo'):
        report_info['gnos_repo'] = gnos_repo
        report_info['gnos_metadata_url'] = gnos_repo + 'cghub/metadata/analysisFull/' + report_info['gnos_id']
        report_info_list.append(copy.deepcopy(report_info))


def add_report_info_1(report_info, report_info_list, es_json):
    if es_json.get('normal_alignment_status'):
        add_report_info_1_aliquot(es_json.get('normal_alignment_status'), report_info, report_info_list)

    if es_json.get('tumor_alignment_status'):
        for aliquot in es_json.get('tumor_alignment_status'):
            add_report_info_1_aliquot(aliquot, report_info, report_info_list)

    return report_info_list

def add_report_info_1_aliquot(aliquot, report_info, report_info_list):
    report_info['aliquot_id'] = aliquot.get('aliquot_id')
    report_info['submitter_specimen_id'] = aliquot.get('submitter_specimen_id')
    report_info['submitter_sample_id'] = aliquot.get('submitter_sample_id')
    report_info['dcc_specimen_type'] = aliquot.get('dcc_specimen_type')   
    report_info['aligned'] = True if aliquot.get('aligned') else False

    if not aliquot.get('do_lane_count_and_bam_count_match') or not aliquot.get('do_lane_counts_in_every_bam_entry_match'):
        report_info['number_of_bams'] = len(aliquot.get('unaligned_bams'))
        report_info['total_lanes'] = aliquot.get('lane_count')
        report_info_list.append(copy.deepcopy(report_info))

    return report_info_list

def add_report_info_2(report_info, report_info_list, es_json):
    if es_json.get('variant_calling_results'):
        vcf = es_json.get('variant_calling_results')
        report_info['normal_bam_gnos_id'] = es_json.get('normal_alignment_status').get('aligned_bam').get('gnos_id') if es_json.get('normal_alignment_status').get('aligned_bam') else None
        report_info['tumor_bam_gnos_id'] = []
        for bam in es_json.get('tumor_alignment_status'):
            if bam.get('aligned_bam'):
                report_info['tumor_bam_gnos_id'].append(bam.get('aligned_bam').get('gnos_id'))
        for workflow in ['sanger', 'embl', 'dkfz', 'dkfz_embl', 'broad', 'muse', 'broad_tar']:
            if vcf.get(workflow+'_variant_calling') and vcf.get(workflow+'_variant_calling').get('is_bam_used_by_' + workflow + '_missing'):
                report_info['workflow_name'] = vcf.get(workflow+'_variant_calling').get('workflow_details').get('variant_workflow_name')
                report_info['vcf_gnos_repo'] = vcf.get(workflow+'_variant_calling').get('gnos_repo')[0]
                report_info['vcf_gnos_id'] = vcf.get(workflow+'_variant_calling').get('gnos_id')  
                report_info['used_normal_bam_gnos_id'] = None
                report_info['used_normal_bam_gnos_url'] = None
                report_info['is_normal_bam_used_by_vcf_missing'] = vcf.get(workflow+'_variant_calling').get('is_normal_bam_used_by_'+workflow+'_missing')
                report_info['used_tumor_bam_gnos_id'] = []
                report_info['used_tumor_bam_gnos_url'] = []
                report_info['is_tumor_bam_used_by_vcf_missing'] = vcf.get(workflow+'_variant_calling').get('is_tumor_bam_used_by_'+workflow+'_missing')    
                for vcf_input in vcf.get(workflow+'_variant_calling').get('workflow_details').get('variant_pipeline_input_info'):
                    if 'normal' in vcf_input.get('attributes').get('dcc_specimen_type').lower():
                        report_info['used_normal_bam_gnos_id'] = vcf_input.get('attributes').get('analysis_id')
                        report_info['used_normal_bam_gnos_url'] = vcf_input.get('attributes').get('analysis_url')
                    if 'tumour' in vcf_input.get('attributes').get('dcc_specimen_type').lower():
                        report_info['used_tumor_bam_gnos_id'].append(vcf_input.get('attributes').get('analysis_id'))
                        report_info['used_tumor_bam_gnos_url'].append(vcf_input.get('attributes').get('analysis_url'))
                
                report_info_list.append(copy.deepcopy(report_info))

    return report_info_list


def add_report_info_3(report_info, report_info_list, es_json):

    duplicate_bwa_bams = es_json.get('duplicated_bwa_alignment_summary')
    if duplicate_bwa_bams.get('exists_mismatch_bwa_bams_in_normal'):
        aliquot = duplicate_bwa_bams.get('normal')
        add_report_info_3_aliquot(aliquot, report_info, report_info_list)   
    
    if duplicate_bwa_bams.get('exists_mismatch_bwa_bams_in_tumor'):
        for aliquot in duplicate_bwa_bams.get('tumor'):
            add_report_info_3_aliquot(aliquot, report_info, report_info_list)
            

def add_report_info_3_aliquot(aliquot, report_info, report_info_list):
    if aliquot.get('exists_mismatch_bwa_bams'):
        report_info['aliquot_id'] = aliquot.get('aliquot_id')
        report_info['dcc_specimen_type'] = aliquot.get('dcc_specimen_type')
        report_info['exists_gnos_id_mismatch'] = aliquot.get('exists_gnos_id_mismatch')
        report_info['exists_md5sum_mismatch'] = aliquot.get('exists_md5sum_mismatch')
        report_info['exists_version_mismatch'] = aliquot.get('exists_version_mismatch')
        aligned_bam = aliquot.get('aligned_bam')
        # find gnos_id which are train2_bams
        report_info['train2_bams_gnos_id'] = set()
        report_info['gnos_id_to_be_reassigned_as_train2_bam'] = set()
        report_info['gnos_id_to_keep'] = set()
        report_info['gnos_id_to_be_removed'] = set()
        gnos_id_all = set()
        gnos_id_without_sanger_call = []
        for bam in aligned_bam:
            gnos_id_all.add(bam.get('gnos_id'))
            if bam.get('is_train2_bam'): 
                report_info['train2_bams_gnos_id'].add(bam.get('gnos_id'))
                report_info['gnos_id_to_keep'].add(bam.get('gnos_id'))
            else:
                if bam.get('is_used_in_sanger_variant_call'):
                    report_info['gnos_id_to_be_reassigned_as_train2_bam'].add(bam.get('gnos_id'))
                    report_info['gnos_id_to_keep'].add(bam.get('gnos_id'))
                else:
                    gnos_id_without_sanger_call.append(bam.get('gnos_id'))

        if not report_info['gnos_id_to_keep']:
            max_num = max(map(gnos_id_without_sanger_call.count, gnos_id_without_sanger_call))
            gnos_id_tmp = [x for x in gnos_id_without_sanger_call if gnos_id_without_sanger_call.count(x) == max_num]
            gnos_id_tmp = gnos_id_tmp[0]
            report_info['gnos_id_to_be_reassigned_as_train2_bam'].add(gnos_id_tmp)
            report_info['gnos_id_to_keep'].add(gnos_id_tmp)

        report_info['gnos_id_to_be_removed'] = gnos_id_all - report_info['gnos_id_to_keep']
        report_info_list.append(copy.deepcopy(report_info))

    return report_info_list


def init_report_dir(metadata_dir, report_name, repo):
    report_dir = metadata_dir + '/reports/' + report_name if not repo else metadata_dir + '/reports/' + report_name + '/' + repo
    if not os.path.exists(report_dir): os.makedirs(report_dir)  # make the folder if not exists
    

    return report_dir


def read_annotations(annotations, type, file_name):
    if not os.path.isfile(file_name):
        return
    with open(file_name, 'r') as r:
        if annotations.get(type): # reset annotation if exists
            del annotations[type]

        if type == 'esad-uk_reheader_uuid':
            annotations[type] = set()
            for line in r:
                if line.startswith('#'): continue
                if len(line.rstrip()) == 0: continue
                annotations[type].add(line.rstrip())
        elif type == 'gender':
            annotations[type] = {}
            reader = csv.DictReader(r, delimiter='\t')
            for row in reader:
                # if not row.get('study_donor_involved_in') == 'PCAWG': continue
                if not row.get('project_code') or not row.get('submitted_donor_id'): continue
                annotations[type][row.get('project_code')+'::'+row.get('submitted_donor_id')] = row.get('donor_sex') if row.get('donor_sex') else None
        elif type == 'gender_update':
            reader = csv.DictReader(r, delimiter='\t')
            for row in reader:
                if not row.get('donor_unique_id'): continue
                annotations['gender'][row.get('donor_unique_id')] = row.get('DCC ICGC21 submitted donor_sex') if row.get('DCC ICGC21 submitted donor_sex') else None

        else:
            print('unknown annotation type: {}'.format(type))
    return annotations

def main(argv=None):

    parser = ArgumentParser(description="Get Donor Info For Specific Query",
             formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--metadata_dir", dest="metadata_dir",
             help="Directory containing metadata manifest files", required=True)
    parser.add_argument("-r", "--gnos_repo", dest="repo",
             help="Specify which GNOS repo to process, process all repos if none specified", required=False)
    parser.add_argument("-q", "--ES_query", dest="q_index",
             help="Specify which ES_query to be used", required=False)


    args = parser.parse_args()
    metadata_dir = args.metadata_dir  # this dir contains gnos manifest files, will also host all reports
    repo = args.repo
    q_index = args.q_index

    if not os.path.isdir(metadata_dir):  # TODO: should add more directory name check to make sure it's right
        sys.exit('Error: specified metadata directory does not exist!')

    q_index = range(len(es_queries)) if not q_index else [int(q_index)] 

    timestamp = str.split(metadata_dir, '/')[-1]
    es_index = 'p_' + ('' if not repo else repo+'_') + re.sub(r'\D', '', timestamp).replace('20','',1)
    es_type = "donor"
    es_host = 'localhost:9200'

    es = Elasticsearch([es_host])
  
    # output result
    report_name = re.sub(r'^generate_', '', os.path.basename(__file__))
    report_name = re.sub(r'\.py$', '', report_name)
    report_dir = init_report_dir(metadata_dir, report_name, repo)

    annotations = {}
    annotations = read_annotations(annotations, 'esad-uk_reheader_uuid', 'esad-uk_uuids.txt')
    annotations = read_annotations(annotations, 'gender', '../pcawg-operations/lists/donor.all_projects.release20.tsv')
    annotations = read_annotations(annotations, 'gender_update', '../pcawg-operations/lists/donor.gender_update.release21.tsv')

    for q in q_index:
        report_tsv_fh = open(report_dir + '/' + es_queries[q].get('name') + '.txt', 'w')  

        # get the list of donors
        donors_list = get_donors_list(es, es_index, es_queries, q)

        report_info_list_full = []
        for donor_unique_id in donors_list:
            # get json doc for each donor                 
            es_json = get_donor_json(es, es_index, donor_unique_id)
            
            report_info_list_donor = create_report_info(donor_unique_id, es_json, q, annotations)

            report_info_list_full.extend(report_info_list_donor)

        # do diff for santa_cruz missing only
        if q in [4, 10, 13, 16]:
            if q==4:
                release_tsv = '../pcawg-operations/data_releases/santa_cruz/santa_cruz_freeze_entry.tsv' 
            elif q==10:
                release_tsv = '../pcawg-operations/data_releases/aug2015/release_aug2015_entry.tsv'
            elif q==13:
                release_tsv = '../pcawg-operations/data_releases/oct2015/release_oct2015_entry.tsv'
            elif q==16:
                release_tsv = '../pcawg-operations/data_releases/mar2016/release_mar2016_entry.tsv'
            else:
                print('No entry for this query!')
            # generate the set of gnos_id
            gnos_id_set = set([l.get('gnos_id') for l in report_info_list_full])
            report_info_list_full = []
            # read bench mark santa_cruz list, hardcode the location of santa_cruz_freeze_json
            with open(release_tsv, 'r') as s:
                reader = csv.DictReader(s, delimiter='\t')
                for row in reader:
                    if not row.get('gnos_id') in gnos_id_set:
                        row_order = OrderedDict()
                        for fn in reader.fieldnames:
                            row_order[fn.strip('#')] = row.get(fn)
                        report_info_list_full.append(row_order)

        report_info_list_full.sort(key=itemgetter('donor_unique_id'))

        header = True  
        for r in report_info_list_full:
            if header:
                report_tsv_fh.write('\t'.join(r.keys()) + '\n')
                header = False 
            # make the list of output from dict
            line = []
            for p in r.keys():
                if isinstance(r.get(p), list):
                    line.append('|'.join(r.get(p)))
                elif isinstance(r.get(p), set):
                    line.append('|'.join(list(r.get(p))))
                elif r.get(p) is None:
                    line.append('')
                else:
                    line.append(str(r.get(p)))
            report_tsv_fh.write('\t'.join(line) + '\n') 
        
        report_tsv_fh.close()            


    return 0


if __name__ == "__main__":
    sys.exit(main())
