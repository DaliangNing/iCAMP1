# iCAMP1
Infer Community Assembly Mechanisms by Phylogenetic bin-based null model analysis (Version 1)

Daliang Ning
2020.8.22
## News

## Key functions in iCAMP package
- iCAMP: Quantify relative importance of basic community assembly processes at both community and phylogenetic group ('bin') levels.  
  - Based on phylogenetic marker gene sequencing results, e.g. OTU or ASV table and phylogenetic tree from 16S sequencing data.
  - The processes including homogeneous and heterogeneous selection, homoginizing and limited dispersal, and 'drift' (drift and other processes)
  - Quantitative for each turnover (between two samples) at community level, and for each phylogenetic bin in a group of samples.
  - Each phylogenetic bin is usually a group of taxa (a few dozens to a few hundreds of OTUs or ASVs) from a family or order.
  - key function: icamp.big
- To implement some other published methods
  - NP: Neutral taxa percentage, i.e. number or relative abundance of taxa following neutral theory model.
    - add options to perform bootstraping test and set multiple metacommunities (different regional pools).
    - key function: snm.comm
  - QPEN: quantifying community assembly processes based on entire-community null model analysis. 
    - add options to handle big datasets and set multiple metacommunities.
    - key function: qpen
  - Our another R package NST has functions to calculate taxonomic and phylogenetic normalized stochasticity ratio (tNST and pNST).
- Some handy functions for big datasets
  - phylogenetic and taxonomic null model analysis at both community and bin levels
    - functions: bNTIn.p, bNTI.bin.big, bNRIn.p, bNRI.bin.big, RC.pc, RC.bin.bigc
  - between-taxa niche difference and phylogenetic distance calculation
    - functions: dniche, pdist.big
  - phylogenetic signal test within phylogenetic groups
    - function: ps.bin
  - midpoint root of big trees
    - function: midpoint.root.big

## Publication
Ning, D., Yuan, M., Wu, L., Zhang, Y., Guo, X., Zhou, X. et al. (2020). A quantitative framework reveals the ecological drivers of grassland soil microbial community assembly in response to warming. bioRxiv. doi:10.1101/2020.02.22.960872.

## How to use
### System requirements

- Operating systems: Windows, or Mac, or Linux, any versions which can run R (version >= 3.2).

- Dependencis: R (version >=3.2; https://www.r-project.org/), R packages: vegan, parallel, permute, ape, bigmemory, nortest.

- iCAMP current version 1.1.3 has been tested on R 3.5.3, R3.6.2, and R4.0.0. 

- Any required non-standard hardware: No.

### Installation guide

- Downlaod and install R (https://www.r-project.org/).

- Install iCAMP.

  - Install published iCAMP: Open R, use function "install.packages" as below. (i am still polishing the package, not submitted yet. 2020.8.15)
  ```
  install.packages("iCAMP")
  ```
  - Install another version:
    - Download an iCAMP version from this repository iCAMP1/RPackage/AllVersions.
    - Open R, install or update following packages: vegan, parallel, permute, ape, bigmemory, nortest.
    ```
    install.packages(c("vegan", "parallel", "permute", "ape", "bigmemory", "nortest"))
    ```   
    - In R, click Packages/install package from local file, then select the file. For windows, select the .zip file. For Mac/Linux, select the .gz file. Alternatively, in Linux sytem, if you open R in a terminal, use following command to install from the .gz file (revise "/Path/to/the/folder" to the real path of the .gz file on your computer, revise "xxx" to the version number of iCAMP):
    ```
    install.packages(pkgs="/Path/to/the/folder/iCAMP_xxx.tar.gz", repos = NULL, type="source")
    ```

- The whole installation typically takes several minutes. Usually, <5 min for R installation, <1 min for the iCAMP package, <2 min for installation of other packages.

### Instructions for use
- Before analyze your own data with iCAMP, you may go through a simple example dataset in the folder /Examples/SimpleOTU.

- When analyzing your own data, check the format of the example data files (otu.txt, tree.nwk, treat2col.txt, and environment.txt) in the folder "SimpleOTU". Revise your data files to the same format. It is fine if you do not have environment factor information, just pay attention to the notes specific to no-env.file situation in the file "icamp.test.r".

- Change the folder paths and file names in the "icamp.test.r" to your data as indicated. 

- Change the thread number for parallel computing, memory limitation, and other parameter setting according to your need. You may check the help document of each function for detailed explanation.

- Run the codes and check the output files in the output folder you've specified. You may check the ReadMe.md in /Examples/SimpleOTU for the meaning of each output file, as well as the help documents in the R package for details.

## End
