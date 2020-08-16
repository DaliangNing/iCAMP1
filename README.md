# iCAMP1
Infer Community Assembly Mechanisms by Phylogenetic bin-based null model analysis (Version 1)

## System requirements

- Operating systems: Windows, or Mac, or Linux, any versions which can run R (version >= 3.2).

- Dependencis: R (version >=3.2; https://www.r-project.org/), R packages: vegan, parallel, permute, ape, bigmemory, nortest.

- iCAMP current version 1.1.3 has been tested on R 3.5.3, R3.6.2, and R4.0.0. 

- Any required non-standard hardware: No.

## Installation guide

- Downlaod and install R (https://www.r-project.org/).

- Install iCAMP.

  - Install published iCAMP: Open R, use function "install.packages" as below.
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

