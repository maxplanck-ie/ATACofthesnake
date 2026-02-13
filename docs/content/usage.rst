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

design specification in comparison entries is always optional. If not specified, the default design will include all factors in the samplesheet, without any interaction terms.
In the above example, the design includes an interaction term between `factor1` and `time`. 
Allowed characters in the design specification are (aside from the tilde `~`) `+`, `*` and `:`.




.. _all-command-line-options:


All command line options
------------------------

.. click:: aos.atac:main
   :prog: ATAC
   :nested: full