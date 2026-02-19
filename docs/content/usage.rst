Usage
=====

After :doc:`installation <installation>`, the workflow can be ran via the `ATAC` command. To get an overview of the available options, run:

.. code:: bash

   ATAC -h

Required files
--------------

The workflow requires at least the following files/paths to run:

  - A directory containing (assumed to be deduplicated) BAM or CRAM files.
  - A GTF file containing the gene models.
  - a FASTA file containing the reference genome sequence.
  - a bed file containing read attracting regions. This should at least contain the mitochondrial chromosome, but can also contain other regions that you want to be filtered out.

All other command line options are depicted in the `section below <all-command-line-options_>`.

Differential accessibility
--------------------------

Samplesheet
^^^^^^^^^^^

Differential accessibility analysis can be performed by specifying a samplesheet and a comparison file. 
The samplesheet should be tab-separated and should (aside from the `sample` column), contain at least one column with a covariate of interest.
Note that per unique combination of covariates, at least two instances (i.e. replicates) need to be present.


======  =======  =======
sample  factor1  time
======  =======  =======
s1      WT       0
s2      WT       0
s3      WT       1
s4      WT       1      
s5      KO       0
s6      KO       0
s7      KO       1
s8      KO       1
======  =======  =======

Keep in mind that when a samplesheet is provided, only samples that are present in that samplesheet are considered for the analysis.
If no samplesheet is provided, all BAM/CRAM files in the specified directory are included in the analysis.

Comparison
^^^^^^^^^^

Twogroup comparisons
""""""""""""""""""""

The actual comparison to be made need to be specified in a separate YAML file. It is possible to specify multiple comparisons in the same file.
Three different types of analyses can be requested, the most basic one is a simple two-group comparison:

.. code:: yaml

   comparison_name:
     type: twogroup
     design: ~factor1*time
     groups:
       my_control_group:
         factor1: WT
       my_treatment_group:
         factor1: KO

design specification in comparison entries is optional. If not specified, the default design will include all factors in the samplesheet, without any interaction terms.
In the above example, the design includes an interaction term between `factor1` and `time` (and expands to ~factor1 + time + factor1:time).
Allowed characters in the design specification are (aside from the tilde `~`) `+`, `*` and `:`.
It's possible to have multiple values for a factor, and it's also possible to have multiple factors specified per group. 

For example this is valid:

.. code:: yaml

   comparison_name:
     type: twogroup
     design: ~celltype*drug
     groups:
       celltypes12:
         celltype:
           - celltype1
           - celltype2
       celltypes34:
         celltype:
           - celltype3
           - celltype4

This will detect the difference between between celltype1 and celltype2 on one hand, and celltype3 and celltype4 on the other hand.
Keep in mind that in the above example, drug is included in the design but not specified in the groups. This means that the effect will be averaged over the different drug conditions. 
In case you want to restrict the comparison to, for example, DMSO treated samples only, you can specify the following:

.. code:: yaml

   comparison_name:
     type: twogroup
     design: ~celltype*drug
     groups:
       celltypes12:
         celltype:
           - celltype1
           - celltype2
         drug: DMSO
       celltypes34:
         celltype:
           - celltype3
           - celltype4
         drug: DMSO

The output for the differential testing will be stored in a folder of their respective type (in this case, `twogroup`), with subfolder named after the comparison_name (first key in an entry).

LRT
"""

A second comparison type that can be requested in the comparison YAML file is an LRT test. 
This allows a simple way to test multiple factors at the same time. Requesting an LRT can be done as such:

.. code:: yaml

   comparison_name:
     type: lrt
     design: ~celltype*drug
     reduced: ~drug

Above design will test for peaks that are differentially accessible between celltypes. 
If celltype itself has only two levels, this is equivalent to a two-group comparison. However, if there are more than two celltypes, this will test for any difference between the celltypes, without specifying which celltypes are different from each other.
Keep in mind that the reduced design should be nested within the full design, meaning that all factors and interaction terms in the reduced design should also be present in the full design.

Timecourse
""""""""""



.. _all-command-line-options:


All command line options
------------------------

.. click:: aos.atac:main
   :prog: ATAC
   :nested: full