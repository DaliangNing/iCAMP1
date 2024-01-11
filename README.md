# iCAMP
Infer Community Assembly Mechanisms by Phylogenetic bin-based null model analysis ([Latest version 1.6.4, 2023-9-4](https://github.com/DaliangNing/iCAMP1/tree/master/RPackage/AllVersions))

Daliang Ning

- Downloaded **18251** times from 2020.9.9 to 2023.7.3.
- Key publications:
-- Ning, D., Yuan, M., Wu, L. et al. A quantitative framework reveals ecological drivers of grassland microbial community assembly in response to warming. Nature Communications 11, 4717 (2020). [https://doi.org/10.1038/s41467-020-18560-z](https://doi.org/10.1038/s41467-020-18560-z) iCAMP was developed in this study
-- Ning D, Wang Y, Fan Y, et al. Environmental stress mediates groundwater microbial community assembly. Nature Microbiology (2024). [https://doi.org/10.1038/s41564-023-01573-x](https://doi.org/10.1038/s41564-023-01573-x) Genearl framework about stress-assembly relationship was proposed in this study
- Recommendation: [NST (stochasticity assessment tool)](https://github.com/DaliangNing/NST) 
## News
- 2024.1.11 Our paper about stress-assembly relationship is published on Nature Microbiology today. [https://doi.org/10.1038/s41564-023-01573-x](https://doi.org/10.1038/s41564-023-01573-x)
- 2023.9.4 iCAMP v1.6.4 is uploaded. More functions for big data; added options for main function to save intermediate results for resuming after unexpected break.
- 2022.6.1 iCAMP v1.5.12 is published on CRAN.
- 2022.4.10 iCAMP v1.5.7 is uploaded to CRAN.
- 2021.10.28 Studies using iCAMP are published on Ecology Letters (https://doi.org/10.1111/ele.13904) and Water Research (https://doi.org/10.1016/j.watres.2021.117744; https://doi.org/10.1016/j.watres.2021.117295)
- 2021.4.18 iCAMP v1.4.3 updated on [github](https://github.com/DaliangNing/iCAMP1/tree/master/RPackage/AllVersions), to allow relative abundances in community matrix and community data transformation.
- 2021.4.15 iCAMP is highlighted on [DOE website](https://www.energy.gov/science/ber/articles/new-approach-helps-determine-how-much-microbial-community-composition-driven).
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

## How to use
### System requirements

- Operating systems: Windows, or Mac, or Linux, any versions which can run R (version >= 3.2).

- Dependencis: R (version >=3.5; https://www.r-project.org/), R packages:vegan,parallel,permute,ape,bigmemory,nortest,minpack.lm,Hmisc,stats4,DirichletReg,data.table.
  - R package NST is necessary to run the funciton tNST and pNST in the example, but not required for running package iCAMP.

- iCAMP current version 1.5.12 has been tested on the current development version of R (4.3.0 pre-release) and R 4.2.1. 

- Any required non-standard hardware: No. However, if you are dealing with a large dataset (e.g. >20,000 taxa), a server with enough CPU threads (e.g. >=20) is preferred to finish the calculation in reasonable time.

### Installation guide

- Downlaod and install R (https://www.r-project.org/).

- Install iCAMP.

  - Install published iCAMP (version<=1.5.12): Open R, use function "install.packages" as below.
  ```
  install.packages("iCAMP")
  ```
  - Install from source file (version>=1.4.1):
    - Download an iCAMP version from this repository iCAMP1/RPackage/AllVersions.
    - Open R, install or update following packages: vegan, parallel, permute, ape, bigmemory, nortest, minpack.lm, Hmisc, stats4, DirichletReg, data.table.
    ```
    install.packages(c("vegan", "permute", "ape", "bigmemory", "nortest", "minpack.lm", "Hmisc", "stats4", "DirichletReg", "data.table"))
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

## Publications
### Our studies
- Ning D, Yuan M, Wu L, Zhang Y, Guo X, Zhou X, Yang Y, Arkin AP, Firestone MK, and Zhou J. 2020. A quantitative framework reveals ecological drivers of grassland microbial community assembly in response to warming. Nature Communications 11, 4717. https://doi.org/10.1038/s41467-020-18560-z.
- Aslani F, Geisen S, Ning D, Tedersoo L, and Bahram M. 2021. Towards revealing the global diversity and community assembly of soil eukaryotes. Ecology Letters https://doi.org/10.1111/ele.13904.
- Wang A, Shi K, Ning D, Cheng H, Wang H, Liu W, Gao S, Li Z, Han J, Liang B, and Zhou J. 2021. Electrical selection for planktonic sludge microbial community function and assembly. Water Research 206, 117744. https://doi.org/10.1016/j.watres.2021.117744.
- Ceja-Navarro JA, Wang Y, Ning D, Arellano A, Ramanculova L, Yuan MM, Byer A, Craven KD, Saha MC, Brodie EL, Pett-Ridge J, and Firestone MK. 2021. Protist diversity and community complexity in the rhizosphere of switchgrass are dynamic as plants develop. Microbiome 9, 96. https://doi.org/10.1186/s40168-021-01042-9.
- Sun C, Zhang B, Ning D, Zhang Y, Dai T, Wu L, Li T, Liu W, Zhou J, and Wen X. 2021. Seasonal dynamics of the microbial community in two full-scale wastewater treatment plants: Diversity, composition, phylogenetic group based assembly and co-occurrence pattern. Water Research 200, 117295. https://doi.org/10.1016/j.watres.2021.117295.

### Other examples
#### 2022
- Wang Y, Li S, Lang X, Huang X, and Su J. 2022. Effects of microtopography on soil fungal community diversity, composition, and assembly in a subtropical monsoon evergreen broadleaf forest of Southwest China. CATENA 211, 106025. https://doi.org/10.1016/j.catena.2022.106025.
- Song T, Liang Q, Du Z, Wang X, Chen G, Du Z, and Mu D. 2022. Salinity Gradient Controls Microbial Community Structure and Assembly in Coastal Solar Salterns. Genes 13, 385. https://doi.org/10.3390/genes13020385.
- Lv B, Shi J, Li T, Ren L, Tian W, Lu X, Han Y, Cui Y, and Jiang T. 2022. Deciphering the characterization, ecological function and assembly processes of bacterial communities in ship ballast water and sediments. Science of the Total Environment 816, 152721. https://doi.org/10.1016/j.scitotenv.2021.152721.
- Ju Z, Du X, Feng K, Li S, Gu S, Jin D, and Deng Y. 2021. The Succession of Bacterial Community Attached on Biodegradable Plastic Mulches During the Degradation in Soil. Frontiers in Microbiology 12 https://doi.org/10.3389/fmicb.2021.785737.
-	Chen S, Tao J, Chen Y, Wang W, Fan L, and Zhang C. 2022. Interactions Between Marine Group II Archaea and Phytoplankton Revealed by Population Correlations in the Northern Coast of South China Sea. Frontiers in Microbiology 12 https://doi.org/10.3389/fmicb.2021.785532.
-	Zhang S, Li K, Hu J, Wang F, Chen D, Zhang Z, Li T, Li L, Tao J, Liu D, and Che R. 2022. Distinct assembly mechanisms of microbial sub-communities with different rarity along the Nu River. Journal of Soils and Sediments https://doi.org/10.1007/s11368-022-03149-4.
#### 2021
- Stopnisek N and Shade A. 2021. Persistent microbiome members in the common bean rhizosphere: an integrated analysis of space, time, and plant genotype. The Isme Journal 15, 2708-2722. https://doi.org/10.1038/s41396-021-00955-5.
- Dong Y, Sanford RA, Connor L, Chee-Sanford J, Wimmer BT, Iranmanesh A, Shi L, Krapac IG, Locke RA, and Shao H. 2021. Differential structure and functional gene response to geochemistry associated with the suspended and attached shallow aquifer microbiomes from the Illinois Basin, IL. Water Research 202, 117431. https://doi.org/10.1016/j.watres.2021.117431.
- Zhu D, Delgado-Baquerizo M, Ding J, Gillings MR, and Zhu Y-G. 2021. Trophic level drives the host microbiome of soil invertebrates at a continental scale. Microbiome 9, 189. https://doi.org/10.1186/s40168-021-01144-4.
- Sun Y, Zhang M, Duan C, Cao N, Jia W, Zhao Z, Ding C, Huang Y, and Wang J. 2021. Contribution of stochastic processes to the microbial community assembly on field-collected microplastics. Environmental Microbiology 23, 6707-6720. https://doi.org/10.1111/1462-2920.15713.
- Macia-Vicente JG and Popa F. 2022. Local endemism and ecological generalism in the assembly of root-colonizing fungi. Ecological Monographs 92, e1489. https://doi.org/10.1002/ecm.1489.
- Yi M, Fang Y, Hu G, Liu S, Ni J, and Liu T. 2021. Distinct community assembly processes underlie significant spatiotemporal dynamics of abundant and rare bacterioplankton in the Yangtze River. Frontiers of Environmental Science & Engineering 16, 79. https://doi.org/10.1007/s11783-021-1513-4.
- Zheng L, Wang X, Ding A, Yuan D, Tan Q, Xing Y, and Xie E. 2021. Ecological Insights Into Community Interactions, Assembly Processes and Function in the Denitrifying Phosphorus Removal Activated Sludge Driven by Phosphorus Sources. Frontiers in Microbiology 12 https://doi.org/10.3389/fmicb.2021.779369.
- Matar GK, Ali M, Bagchi S, Nunes S, Liu W-T, and Saikaly PE. 2021. Relative Importance of Stochastic Assembly Process of Membrane Biofilm Increased as Biofilm Aged. Frontiers in Microbiology 12 https://doi.org/10.3389/fmicb.2021.708531.
- Yuan H, Li T, Li H, Wang C, Li L, Lin X, and Lin S. 2021. Diversity Distribution, Driving Factors and Assembly Mechanisms of Free-Living and Particle-Associated Bacterial Communities at a Subtropical Marginal Sea. Microorganisms 9, 2445. https://doi.org/10.3390/microorganisms9122445.
- Wang Y, Lu G, Yu H, Du X, He Q, Yao S, Zhao L, Huang C, Wen X, and Deng Y. 2021. Meadow degradation increases spatial turnover rates of the fungal community through both niche selection and dispersal limitation. Science of the Total Environment 798, 149362. https://doi.org/10.1016/j.scitotenv.2021.149362.
- Xie J, Wang X, Xu J, Xie H, Cai Y, Liu Y, and Ding X. 2021. Strategies and Structure Feature of the Aboveground and Belowground Microbial Community Respond to Drought in Wild Rice (Oryza longistaminata). Rice 14, 79. https://doi.org/10.1186/s12284-021-00522-8.
- Zhou R, Wang H, Wei D, Zeng S, Hou D, Weng S, He J, and Huang Z. 2021. Bacterial and eukaryotic community interactions might contribute to shrimp culture pond soil ecosystem at different culture stages. Soil Ecology Letters https://doi.org/10.1007/s42832-021-0082-6.
 
## End
