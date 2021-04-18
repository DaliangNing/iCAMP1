# iCAMP
Infer Community Assembly Mechanisms by Phylogenetic bin-based null model analysis (Latest version 1.3.4)

Daliang Ning
- Downloaded **4050** times from 2020.9.9 to 2021.4.17.
- Recommendation: [NST (stochasticity assessment tool)](https://github.com/DaliangNing/NST) 
## News
- 2021.4.18 iCAMP v1.4.3 updated on [github](https://github.com/DaliangNing/iCAMP1/tree/master/RPackage/AllVersions), to allow relative abundances in community matrix and community data transformation.
- 2021.4.1 Frontiers in Microbilogy opens a research topic [**Community Assembly Mechanisms Shaping Microbiome Spatial or Temporal Dynamics**](https://www.frontiersin.org/research-topics/20916/).
- 2021.1.9 iCAMP v1.3.4 updated on CRAN, improved function qpen and added function icamp.cate to summary for different categories (e.g. core vs rare taxa). New example is updated in the subfolder [**Examples**](https://github.com/DaliangNing/iCAMP1/tree/master/Examples/SimpleOTU).
- 2020.9.22 Media reports: [OU-VPR news](https://ou.edu/research-norman/news-events/2020/study-expands-understanding-of-ecological-factors-driving-microbial-community-assembly-in-response-to-warming), [EurekAlert!](https://www.eurekalert.org/pub_releases/2020-09/uoo-efd092220.php), [Phy.Org](https://phys.org/news/2020-09-ecological-factors-microbial-response.html)
- 2020.9.18 iCAMP paper is published on Nature Communications. https://doi.org/10.1038/s41467-020-18560-z
- 2020.9.9  iCAMP v1.2.8 has been published on CRAN. https://cran.r-project.org/web/packages/iCAMP/
- 2020.8.25 iCAMP v1.2.5 fixed some typo and memory.limit issue, and is submitted back to CRAN.
- 2020.8.24 iCAMP v1.2.4 has been submitted to CRAN.
- 2020.8.23 upload iCAMP package (v1.2.4) and the code/data for the first iCAMP manuscript to GitHub.

## Key functions in iCAMP package
- [**iCAMP**](https://doi.org/10.1101/2020.02.22.960872): Quantify relative importance of basic community assembly processes at both community and phylogenetic group ('bin') levels.  
  - Based on phylogenetic marker gene sequencing results, e.g. OTU or ASV table and phylogenetic tree from 16S sequencing data.
  - The processes including homogeneous and heterogeneous selection, homoginizing and limited dispersal, and 'drift' (drift and other processes)
  - Quantitative for each turnover (between two samples) at community level, and for each phylogenetic bin in a group of samples.
  - Each phylogenetic bin is usually a group of taxa (a few dozens to a few hundreds of OTUs or ASVs) from a family or order.
  - key function: icamp.big ([Ning et al 2020 Nat Commun in revision](https://doi.org/10.1101/2020.02.22.960872))
- To implement some other published methods
  - [**NP**](https://doi.org/10.1038/ismej.2015.142): Neutral taxa percentage, i.e. number or relative abundance of taxa following neutral theory model.
    - developed by [Burns et al (2016 ISME J)](https://doi.org/10.1038/ismej.2015.142), based on a neutral theory model [(Sloan et al 2006 EM)](https://doi.org/10.1111/j.1462-2920.2005.00956.x).
    - I add options to perform bootstraping test and re-define taxa abundance profile in one or multiple metacommunities (regional pools).
    - function: snm.comm
  - [**QPEN**](https://doi.org/10.1038/ismej.2013.93): quantifying community assembly processes based on entire-community null model analysis. 
    - developed by Stegen et al ([2013 ISME J](https://doi.org/10.1038/ismej.2013.93), [2015 Front Microbiol](https://doi.org/10.3389/fmicb.2015.00370)).
    - I add options to handle big datasets and re-define taxa abundance profile in the metacommunity.
    - function: qpen
  - [**tNST** and **pNST**](https://doi.org/10.1073/pnas.1904623116): taxonomic and phylogenetic normalized stochasticity ratio.
    - Not in iCAMP, but in our another R package [NST](https://cran.r-project.org/web/packages/NST/index.html)
    - We developed NST ([Ning et al 2019 PNAS](https://doi.org/10.1073/pnas.1904623116)) from previous stochasticity ratio ([Zhou et al 2014 PNAS](https://doi.org/10.1073/pnas.1324044111)).
- Some handy functions for big datasets
  - phylogenetic and taxonomic **null model analysis** at both community and bin levels
    - functions: bNTIn.p, bNTI.bin.big, bNRIn.p, bNRI.bin.big, RC.pc, RC.bin.bigc
  - between-taxa **niche difference** and **phylogenetic distance** of big communities
    - functions: dniche, pdist.big
  - **phylogenetic signal** test within phylogenetic groups
    - function: ps.bin
  - **midpoint root** of big trees
    - function: midpoint.root.big

## Publication
Ning D, Yuan M, Wu L, Zhang Y, Guo X, Zhou X, Yang Y, Arkin AP, Firestone MK, and Zhou J. A quantitative framework reveals ecological drivers of grassland microbial community assembly in response to warming. Nature Communications 11, 4717 (2020) [doi:10.1038/s41467-020-18560-z](https://doi.org/10.1038/s41467-020-18560-z).

## How to use
### System requirements

- Operating systems: Windows, or Mac, or Linux, any versions which can run R (version >= 3.2).

- Dependencis: R (version >=3.2; https://www.r-project.org/), R packages: vegan, parallel, permute, ape, bigmemory, nortest, minpack.lm, Hmisc, stats4.
  - R package NST is necessary to run the funciton tNST and pNST in the example, but not required for running package iCAMP.

- iCAMP current version 1.2.4 has been tested on the current development version of R (4.1.0, 2020-8-18 r79041), R 4.0.2, and R 3.5.3. 

- Any required non-standard hardware: No. However, if you are dealing with a large dataset (e.g. >20,000 taxa), a server with enough CPU threads (e.g. >=20) is preferred to finish the calculation in reasonable time.

### Installation guide

- Downlaod and install R (https://www.r-project.org/).

- Install iCAMP.

  - Install published iCAMP: Open R, use function "install.packages" as below.
  ```
  install.packages("iCAMP")
  ```
  - Install from source file:
    - Download an iCAMP version from this repository iCAMP1/RPackage/AllVersions.
    - Open R, install or update following packages: vegan, parallel, permute, ape, bigmemory, nortest.
    ```
    install.packages(c("vegan", "permute", "ape", "bigmemory", "nortest", "minpack.lm", "Hmisc"))
    ```   
    - In R, click Packages/install package from local file, then select the file. For windows, select the .zip file. For Mac/Linux, select the .gz file. Alternatively, in Linux sytem, if you open R in a terminal, use following command to install from the .gz file (revise "/Path/to/the/folder" to the real path of the .gz file on your computer, revise "xxx" to the version number of iCAMP):
    ```
    install.packages(pkgs="/Path/to/the/folder/iCAMP_xxx.tar.gz", repos = NULL, type="source")
    ```

- The whole installation typically takes several minutes. Usually, <5 min for R installation, <1 min for the iCAMP package, <5 min for installation of other packages.

### Instructions for use
- Before analyze your own data with iCAMP, you may go through a simple example dataset in the folder /Examples/SimpleOTU.

- When analyzing your own data, check the format of the example data files (otu.txt, tree.nwk, treat2col.txt, and environment.txt) in the folder "SimpleOTU". Revise your data files to the same format. It is fine if you do not have environment factor information, just pay attention to the notes specific to no-env.file situation in the file "icamp.test.r".

- Change the folder paths and file names in the "icamp.test.r" to your data as indicated. 

- Change the thread number for parallel computing, memory limitation, and other parameter setting according to your need. You may check the help document of each function for detailed explanation.

- Run the codes and check the output files in the output folder you've specified. You may check the ReadMe.md in /Examples/SimpleOTU for the meaning of each output file, as well as the help documents in the R package for details.

## End
