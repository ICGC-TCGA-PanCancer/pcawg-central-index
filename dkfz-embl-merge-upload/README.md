SOP for preparing GNOS submissions for DKFZ data fixes and merging with EMBL
=====

Summary
-----
This SOP describes a process to generate DKFZ/EMBL combined GNOS submissions that will include fixed DKFZ variant
calling data files and the original EMBL results.


Main steps
-----

Each of the following step will need to be finished completely before moving on to the next step.

### 1. Download DKFZ and EMBL variant call results ###

You need to organize the data in a *working directory*, this directory contains all data files downloaded from GNOS
and files locally fixed. It is required to structure the directory as following example shows (this example only shows
data folders for two DKFZ GNOS entries and two EMBL GNOS entries):

```
working_dir
└── downloads
    ├── dkfz
    │   ├── 0eaa9e80-c9fb-4459-895d-00bd897e69ec
    │   │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-copyNumberEstimation_1-0-132-1.20150604.somatic.cnv.tar.gz
    │   │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-copyNumberEstimation_1-0-132-1.20150604.somatic.cnv.vcf.gz
    │   │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-copyNumberEstimation_1-0-132-1.20150604.somatic.cnv.vcf.gz.tbi
    │   │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-indelCalling_1-0-132-1.20150604.germline.indel.vcf.gz
    │   │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-indelCalling_1-0-132-1.20150604.somatic.indel.tar.gz
    │   │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-indelCalling_1-0-132-1.20150604.somatic.indel.vcf.gz
    │   │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-indelCalling_1-0-132-1.20150604.somatic.indel.vcf.gz.tbi
    │   │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-snvCalling_1-0-132-1.20150604.germline.snv_mnv.vcf.gz
    │   │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-snvCalling_1-0-132-1.20150604.germline.snv_mnv.vcf.gz.tbi
    │   │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-snvCalling_1-0-132-1.20150604.somatic.snv_mnv.tar.gz
    │   │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-snvCalling_1-0-132-1.20150604.somatic.snv_mnv.vcf.gz
    │   │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-snvCalling_1-0-132-1.20150604.somatic.snv_mnv.vcf.gz.tbi
    │   │   └── fixed_files
    │   │       ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-indelCalling_1-0-132-1.20150604.germline.indel.vcf.gz.tbi
    │   │       ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-copyNumberEstimation_1-0-132-1.20150814.somatic.cnv.tar.gz
    │   │       ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-copyNumberEstimation_1-0-132-1.20150814.somatic.cnv.vcf.gz
    │   │       └── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-copyNumberEstimation_1-0-132-1.20150814.somatic.cnv.vcf.gz.tbi
    │   └── f1894948-a49b-410f-9009-8e960c6b5565 [this is a folder named using GNOS ID]
    │       ├── all_dkfz_gnos_original_files_go_here
    │       └── fixed_files [this is the folder containing fixed files]
    │           └── locally_fixed_files_go_here
    └── embl
        ├── c3d41bd1-0194-4341-a5dd-a15878392416
        │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.germline.sv.bedpe.txt.tar.gz
        │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.germline.sv.readname.txt.tar.gz
        │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.germline.sv.vcf.gz
        │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.germline.sv.vcf.gz.tbi
        │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.somatic.sv.bedpe.txt.tar.gz
        │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.somatic.sv.readname.txt.tar.gz
        │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.somatic.sv.vcf.gz
        │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.somatic.sv.vcf.gz.tbi
        │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.sv.cov.plots.tar.gz
        │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.sv.cov.tar.gz
        │   ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.sv.vcf.gz
        │   └── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.sv.vcf.gz.tbi
        └── ca1a99a7-8875-4262-8c10-c1a75b60c039 [this is a folder named using GNOS ID]
            └── all_embl_gnos_original_files_go_here
```

The list of first batch of donors whose DKFZ call result need to be fixed is packaged with this
tool in file: `to_be_fixed_donor_list.txt`. There are totally 742 donors, it is expected
to have all DKFZ variant call result files downloaded from GNOS to under `working_dir/downloads/dkfz/`,
similarly EMBL results will be under `working_dir/downloads/embl/`. 


### 2. Add fixed DKFZ files to `fixed_files` directory ###

Fixed files are expected to be saved in corresponding `fixed_files` directory under where the original
GNOS files are located. Such as `downloads/dkfz/0eaa9e80-c9fb-4459-895d-00bd897e69ec/fixed_files` in the
above example. Files in this directory will be used for new GNOS submission, they will replace the original
files that have the same naming pattern.

IMPORTANT: Please ensure all fixed files for all donors are added to the expected folders before moving on to the next step!

Note: It was discovered that the `tbi` file for germline indel VCF file was missing. It would be appreciated if DKFZ can
generate the index file and add it in `fixed_files` folder, it will then be included to new submission.


### 3. Run the Python script to generate new GNOS submissions ###

This tool includes one single Python script: `metadata_fix_and_merge.py`, it runs with working directory as argument:
```
python metadata_fix_and_merge.py <working_dir>
```

You can try the script out using *test* as the 'working_dir'. This should be useful to get yourself familiar how it works.


Here is how the script works:

- examine the working directory to ensure all DKFZ and EMBL entries exist, and `fixed_files` subfolder exists for each DKFZ GNOS entry
- download metadata XML for all the above GNOS entries and save them in the perspective folder
- scan the `fixed_files` directory to pickup information for the new DKFZ data files and update GNOS metadata XML
- create `uploads` directory (same place as `downloads` located) and folders with new UUID for storing new GNOS submissions
- create `analysis.xml` file for new submission by merging updated DKFZ metadata with original EMBL metadata
- build symlinks in the new submission linking to files from DKFZ and EMBL download folders

The `uploads` directory is structured as the following example:

```
working_dir
├── downloads
└── uploads
    └── ce914762-9f61-4402-b5fc-196e7cfe1ab8 [new UUID for new GNOS submission]
        ├── analysis.xml
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-copyNumberEstimation_1-0-132-1.20150814.somatic.cnv.tar.gz
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-copyNumberEstimation_1-0-132-1.20150814.somatic.cnv.vcf.gz
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-copyNumberEstimation_1-0-132-1.20150814.somatic.cnv.vcf.gz.tbi
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-indelCalling_1-0-132-1.20150604.germline.indel.vcf.gz
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-indelCalling_1-0-132-1.20150814.germline.indel.vcf.gz.tbi
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-indelCalling_1-0-132-1.20150604.somatic.indel.tar.gz
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-indelCalling_1-0-132-1.20150604.somatic.indel.vcf.gz
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-indelCalling_1-0-132-1.20150604.somatic.indel.vcf.gz.tbi
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-snvCalling_1-0-132-1.20150604.germline.snv_mnv.vcf.gz
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-snvCalling_1-0-132-1.20150604.germline.snv_mnv.vcf.gz.tbi
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-snvCalling_1-0-132-1.20150604.somatic.snv_mnv.tar.gz
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-snvCalling_1-0-132-1.20150604.somatic.snv_mnv.vcf.gz
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.dkfz-snvCalling_1-0-132-1.20150604.somatic.snv_mnv.vcf.gz.tbi
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.germline.sv.bedpe.txt.tar.gz
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.germline.sv.readname.txt.tar.gz
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.germline.sv.vcf.gz
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.germline.sv.vcf.gz.tbi
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.somatic.sv.bedpe.txt.tar.gz
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.somatic.sv.readname.txt.tar.gz
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.somatic.sv.vcf.gz
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.somatic.sv.vcf.gz.tbi
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.sv.cov.plots.tar.gz
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.sv.cov.tar.gz
        ├── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.sv.vcf.gz
        └── d5423a93-0a72-43cf-a7ae-9011f47553c7.embl-delly_1-0-0-preFilter.20150604.sv.vcf.gz.tbi
```


IMPORTANT: To prevent possible duplicated uploads, it is expected to run the script ONLY ONCE to have ALL new submissions
created. In case the script failed all previously generated submissions MUST be deleted before running the script again.

Note: A log file will be produced for each run, the file name starts with a timestamp, e.g., 2015-08-15_16-53-46.process.log.
The log file provides useful information for debugging when failure occurs.


### 4. Submit / Upload new submissions to GNOS ###

With all submissions created, please follow usual GNOS submission steps to submit/upload metadata/data to GNOS. This consists
of two major steps, metadata submission using `cgsubmit` and data uploading using `gtupload`.

