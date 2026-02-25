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

The default output (which will always be generated) is as such:

 - Directory `bw` containing RPKM and scalefactor (calculated on signal only) bigwig files for each sample.
 - Directory `figures`, containing the read statistics before and after filtering (fragmentsizes.png and alignmentsieve.png, respectively). The FRIP scores (fripscores.png), fraction of mitochondrial reads (mitofraction.png) and a PCA plot of the samples.
 - Directory `peaks`, containing the called peaks for each sample.
 - Directory `peakset`, containing the unionized peak set, a countmatrix, the calculated scalefactors and the peak-to-gene annotations. 
 - Directory `qc`, contains the numeric data underlying the figures
 - Directory `sieve`, contains the filtered (blacklist / NFR) BAM files. 

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
Included are a table with the results of the differential testing (standard edgeR output), which also includes a column with the assigned group label (based on the sign of the estimated log fold change) and the actual peak ID.
Keep in mind that the group label is set irrespective of the FDR, so filter first before you interpret the group labels. Additionally, an MAplot is generated, and potentially a heatmap too (depending on the number of differential sites).


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

The output here is similar to the two-group comparison output (though under the `lrt` subfolder), with the difference that there is no group label in the edgeR output. 
Instead, a rudimentary k-means clustering is ran on the significant peaks (calculating k based on the intertia metrics). This is applied to the heatmap (which again is only generated if there are sufficient significant peaks to begin with).


Timecourse
""""""""""

A third comparison type is a timecourse analysis. Two different timecourse modes are possible, with both of them based on Gaussian process regression.
Note that timecourse analyses could also be performed by using the `LRT` comparison type, and dropping the time factor from the reduced design. 
The assumptions made there are however different: if time is encoded as a continuous variable, the effect is presumed to be monotonic with how time is encoded in the design (lienar, quadratic, etc.).
If time is encoded as a categorical variable, then the effect is not presumed to be monotonic, but the inherent ordering is not taking into account.
To tackle this, the timecourse comparison allows a setup with time being treated as a continuous variable, but without presuming monotonicity. 
This can be specified as such:

.. code:: yaml

    comparison_name:
      type: timecourse
      time: 'time'
      time_type: 'continuous'

Note that here you cannot specify a design, and all covariates in the samplesheet will be included in the analysis by default. 
Additionally, note that testing differential time trajectories across other covariates (i.e. time - covariate interaction term) is not yet implemented.
Running the above comparison will run a Gaussian process regression per peak, and test whether the result is substantially different from a flat kernel by taking the marginal likelihood ratio between them.
Since this cannot be directly tested with a 'regular' chi-squared test, a permutation-based approach is used to calculate a p-value, that is afterwards correct with a Benjamini-Hochberg correction for multiple testing.
The output is written to the `gp` folder, with the main results being written to a table with the likelihood ratio, p-value and FDR for each peak. Additionally, an array with predicted standardized accessibility scores (y_pred) and one with the actual standardized values (y_std) are included in the table.
A similar table, but now only with significant peaks (by default FDR < 0.05) also exists, with the addition that here the peaks are clustered (k-means setup, similar to the one performed in the LRT comparison type). A visualization of these clusters, and their individual peak trajectories (standardized) is included in the folder too.

A second timecourse mode can be used when time is not encoded as a continuous variable, but as an ordinal variable. Typical use cases for this mode are differentiation trajectories.
Here it's important that order needs to be specified, with the 'earliest' time point being the first. Note that not all possible values of the time factor listed in the samplesheet need to be included here (samples with missing levels are simply ignored).

.. code:: yaml

    comparison_name:
      type: timecourse
      time: 'celltype'
      time_type: 'ordinal'
      order:
        - 'preT_DN1'
        - 'preT_DN2b'
        - 'preT_DN3'
        - 'TDP'
        - 'T4'

A similar Gaussian process regression is ran here as in the continuous case, with the exception that the timepoints are re-encoded between 0 and 1, and the distance between subsequent levels is included in the kernel.
Similar to the output in the continuous case, with the addition of the estimated relative distance between time points per peak in a separate table.

.. _all-command-line-options:


All command line options
------------------------

.. click:: aos.atac:main
   :prog: ATAC
   :nested: full