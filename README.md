# Single-cell RNA-seq pseudotime estimation algorithms
[![Test links](https://github.com/agitter/single-cell-pseudotime/actions/workflows/links.yml/badge.svg)](https://github.com/agitter/single-cell-pseudotime/actions/workflows/links.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1297422.svg)](https://doi.org/10.5281/zenodo.1297422)

Single cells, many algorithms.
The goal of this page is to catalog the many algorithms that estimate pseudotimes for cells based on their gene expression levels.
This problem is also referred to as single-cell trajectory inference or ordering.
It contains method names, software links, and manuscript links and simply seeks to list as many methods as possible without commentary.
Some related methods not specifically designed for RNA-seq (e.g. mass cytometry) are included as well, as are some methods for estimating RNA velocity.
The list also includes methods that are specifically designed to take pseudotemporal data as input.

The initial list was created by Anthony Gitter, but pull requests are very welcome.
Thank you to the other [contributors](https://github.com/agitter/single-cell-pseudotime/graphs/contributors).


## Citation

Anthony Gitter. Single-cell RNA-seq pseudotime estimation algorithms. 2018.
doi:10.5281/zenodo.1297422
https://github.com/agitter/single-cell-pseudotime


## Problem overview

Informally, the pseudotime estimation problem can be stated as:
- **Given:** single-cell gene expression measurements for a heterogeneous collection of cells that is transitioning from biological state **A** to state **B**
- **Return:** a quantitative value for each cell that represents its progress in the **A** to **B** transition

There are many ways to approach this problem, and major algorithmic steps that are common to most (but not all) methods are:
- Reduce the dimension of the dataset
- Find a smooth progression through the low dimensional data, assuming that cells that are nearby one another in the low dimensional space have similar expression levels because they are at similar points in to **A** to **B** process

Dimension reduction sometimes relies on knowledge of important marker genes and sometimes uses the full gene expression matrix.
The trajectory through the low dimensional space can be identified using graph algorithms (e.g., minimum spanning tree or shortest path), principal curves, or probabilistic models (e.g., Gaussian process).

[Bacher and Kendziorski 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0927-y), [Trapnell 2015](http://genome.cshlp.org/content/25/10/1491.full), [Tanay and Regev 2017](http://www.nature.com/nature/journal/v541/n7637/full/nature21350.html), [Moon et al. 2017](https://doi.org/10.1016/j.coisb.2017.12.008), [Tritschler et al. 2019](https://doi.org/10.1242/dev.170506), [Weiler et al. 2021](https://doi.org/10.1101/2021.12.22.473434), and [Ding et al. 2022](https://doi.org/10.1038/s41576-021-00444-7) provide a more comprehensive overview of single-cell RNA-seq and the pseudotime estimation problem.
[Cannoodt et al. 2016](http://onlinelibrary.wiley.com/wol1/doi/10.1002/eji.201646347/abstract) reviews pseudotime inference algorithms.
[Pablo Cordero's blog post](http://hyperparameter.space/blog/a-single-cell-journey-from-mechanistic-to-descriptive-modeling-and-back-again/) discusses why it is hard to recover the correct dynamics of a system from single-cell data.
For more general lists of methods for single-cell RNA-seq see [seandavi/awesome-single-cell](https://github.com/seandavi/awesome-single-cell) and [scRNA-tools](https://www.scrna-tools.org/).
The Hemberg lab [single-cell RNA-seq course](https://scrnaseq-course.cog.sanger.ac.uk/website/index.html) has an [overview of five pseudotime algorithms](https://scrnaseq-course.cog.sanger.ac.uk/website/biological-analysis.html#pseudotime-analysis) with usage examples.
Many modern ideas for pseudotime estimation are descended from [Magwene et al. 2003](https://doi.org/10.1093/bioinformatics/btg081) on reconstructing the order of microarray expression samples.

Single-cell expression data have also inspired new methods for gene regulatory network reconstruction, as reviewed by [Fiers et al. 2018](https://doi.org/10.1093/bfgp/elx046) and [Todorov et al. 2018](https://doi.org/10.1007/978-1-4939-8882-2_10).
Several of these, such as [SINGE](http://doi.org/10.1016/j.celrep.2022.110333), treat pseudotime annotations as time points and extend traditional time series network inference algorithms for single-cell data.
[BEELINE](http://doi.org/10.1038/s41592-019-0690-6), [SERGIO](http://doi.org/10.1016/j.cels.2020.08.003), and [McCalla et al. 2023](https://doi.org/10.1093/g3journal/jkad004) benchmark many of these specialized network inference methods.

## Choosing a method

Some of the distinguishing factors among algorithms include:
- Use of prior knowledge such as capture times (DeLorean) or switch-like marker genes (Ouija)
- Modeling specific types of biological processes such as branching processes in differentiation (multiple methods) or cyclic processes (Oscope)
- Return a single pseudotime or a posterior distribution over pseudotimes for each cell
- Perform additional analyses after inferring pseudotimes such as regulatory network inference or identifying differentially expressed genes over pseudotime

[Saelens et al. 2019](https://doi.org/10.1038/s41587-019-0071-9) performed a comprehensive evaluation of 29 different single-cell trajectory inference methods and discuss the different types of algorithms in more detail.
They benchmark both quantitative performance and assess software quality.
See their [website](https://dynverse.org/) and [GitHub repository](https://github.com/dynverse/dynverse) as well.
[Tian et al. 2018](https://doi.org/10.1101/433102) also include trajectory inference algorithms in their single-cell RNA-seq benchmarking study.
[Escort](https://doi.org/10.1101/2023.12.18.572214) is a framework to help guide the selection of a suitable trajectory inference algorithm for a dataset.

## Algorithms

### Temporal Reconstruction Algorithm
Manuscript: [Reconstructing the temporal ordering of biological samples using microarray data](https://doi.org/10.1093/bioinformatics/btg081)

### Monocle / Monocle 2 / Monocle 3 / Census
Software: https://bioconductor.org/packages/release/bioc/html/monocle.html

Monocle manuscript: [The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells](https://doi.org/10.1038/nbt.2859)

Census manuscript: [Single-cell mRNA quantification and differential analysis with Census](https://doi.org/10.1038/nmeth.4150)

Monocle 2 manuscript: [Reversed graph embedding resolves complex single-cell trajectories](https://doi.org/10.1038/nmeth.4402)

Monocle 3 manuscript: [The single-cell transcriptional landscape of mammalian organogenesis](https://doi.org/10.1038/s41586-019-0969-x)

### Wanderlust / Cycler / Wishbone
Wanderlust software: http://www.c2b2.columbia.edu/danapeerlab/html/wanderlust.html

Wanderlust manuscript: [Single-Cell Trajectory Detection Uncovers Progression and Regulatory Coordination in Human B Cell Development](https://doi.org/10.1016/j.cell.2014.04.005)

Cycler manuscript: [Trajectories of cell-cycle progression from fixed cell populations](https://doi.org/10.1038/nmeth.3545)

Wishbone software: http://www.c2b2.columbia.edu/danapeerlab/html/wishbone.html

Wishbone manuscript: [Wishbone identifies bifurcating developmental trajectories from single-cell data](https://doi.org/10.1038/nbt.3569)

### SCUBA
Software: https://github.com/gcyuan/SCUBA

Manuscript: [Bifurcation analysis of single-cell gene expression data reveals epigenetic landscape](https://doi.org/10.1073/pnas.1408993111)

### Oscope
Software: https://www.biostat.wisc.edu/~kendzior/OSCOPE/

Manuscript: [Oscope identifies oscillatory genes in unsynchronized single-cell RNA-seq experiments](https://doi.org/10.1038/nmeth.3549)

### Diffusion maps / destiny
destiny software: http://bioconductor.org/packages/release/bioc/html/destiny.html

Diffusion maps manuscript (a): [Decoding the regulatory network of early blood development from single-cell gene expression measurements](https://doi.org/10.1038/nbt.3154)

Diffusion maps manuscript (b): [Diffusion maps for high-dimensional single-cell analysis of differentiation data](https://doi.org/10.1093/bioinformatics/btv325)

destiny manuscript: [destiny: diffusion maps for large-scale single-cell data in R](https://doi.org/10.1093/bioinformatics/btv715)

### DeLorean
Software: https://github.com/JohnReid/DeLorean

Manuscript: [Pseudotime estimation: deconfounding single cell time series](https://doi.org/10.1093/bioinformatics/btw372)

### Waterfall
Manuscript: [Single-Cell RNA-Seq with Waterfall Reveals Molecular Cascades underlying Adult Neurogenesis](https://doi.org/10.1016/j.stem.2015.07.013)

### Embeddr
Software: https://github.com/kieranrcampbell/embeddr

Manuscript: [Laplacian eigenmaps and principal curves for high resolution pseudotemporal ordering of single-cell RNA-seq profiles](https://doi.org/10.1101/027219)

### GP-LVM and pseudogp
GP-LVM software: https://github.com/kieranrcampbell/gpseudotime

GP-LVM manuscript: [Bayesian Gaussian Process Latent Variable Models for pseudotime inference in single-cell RNA-seq data](https://doi.org/10.1101/026872)

pseudogp software: https://github.com/kieranrcampbell/pseudogp

pseudogp manuscript: [Order Under Uncertainty: Robust Differential Expression Analysis Using Probabilistic Models for Pseudotime Inference](https://doi.org/10.1371/journal.pcbi.1005212)

### GP-LVM
Analysis code: https://github.com/Teichlab/spectrum-of-differentiation-supplements

Manuscript: [Single-Cell RNA-Sequencing Reveals a Continuous Spectrum of Differentiation in Hematopoietic Cells](https://doi.org/10.1016/j.celrep.2015.12.082)

### SLICER
Software: https://github.com/jw156605/SLICER

Manuscript: [SLICER: inferring branched, nonlinear cellular trajectories from single cell RNA-seq data](https://doi.org/10.1186/s13059-016-0975-3)

### TSCAN
Software: https://github.com/zji90/TSCAN

Manuscript: [TSCAN: Pseudo-time reconstruction and evaluation in single-cell RNA-seq analysis](https://doi.org/10.1093/nar/gkw430)

### SCOUP
Software: https://github.com/hmatsu1226/SCOUP

Manuscript: [SCOUP: a probabilistic model based on the Ornstein–Uhlenbeck process to analyze single-cell expression data during differentiation](https://doi.org/10.1186/s12859-016-1109-3)

### Topslam
Software: https://github.com/mzwiessele/topslam

Manuscript: [Topslam: Waddington Landscape Recovery for Single Cell Experiments
](https://doi.org/10.1101/057778)

### Ouija
Software: https://github.com/kieranrcampbell/ouija and http://www.github.com/kieranrcampbell/ouijaflow

Manuscript: [A descriptive marker gene approach to single-cell pseudotime inference](https://doi.org/10.1093/bioinformatics/bty498)

### CellTrails

Sofware: https://bioconductor.org/packages/release/bioc/html/CellTrails.html

Manuscript: [Transcriptional dynamics of hair-bundle morphogenesis revealed with CellTrails](https://doi.org/10.1016/j.celrep.2018.05.002)

### Slingshot
Software: https://github.com/kstreet13/slingshot

Extended vignette: https://github.com/drisso/bioc2016singlecell/tree/master/vignettes

Manuscript: [Slingshot: Cell lineage and pseudotime inference for single-cell transcriptomics](https://doi.org/10.1101/128843)

Workflow manuscript: [Bioconductor workflow for single-cell RNA sequencing: Normalization, dimensionality reduction, clustering, and lineage inference](https://doi.org/10.12688/f1000research.12122.1)

### GPfates

Software: https://github.com/Teichlab/GPfates

Manuscript: [Temporal mixture modelling of single-cell RNA-seq data resolves a CD4+ T cell fate bifurcation](https://doi.org/10.1101/074971)

### SCIMITAR

Software: https://github.com/dimenwarper/scimitar

Manuscript: [Tracing co-regulatory network dynamics in noisy, single-cell transcriptome trajectories](https://doi.org/10.1101/070151)

### WaveCrest

Software: https://github.com/lengning/WaveCrest

Manuscript: [Single-cell RNA-seq reveals novel regulators of human embryonic stem cell differentiation to definitive endoderm](https://doi.org/10.1186/s13059-016-1033-x)


### CellTree

Software: http://bioconductor.org/packages/release/bioc/html/cellTree.html

Manuscript: [CellTree: an R/bioconductor package to infer the hierarchical structure of cell populations from single-cell RNA-seq data](https://doi.org/10.1186/s12859-016-1175-6)

### Bayesian hierarchical mixture of factor analysers (MFA)

Software: http://www.github.com/kieranrcampbell/mfa

Manuscript: [Probabilistic modeling of bifurcations in single-cell gene expression data using a Bayesian mixture of factor analyzers](http://doi.org/10.12688/wellcomeopenres.11087.1)

### Mpath

Software: https://github.com/JinmiaoChenLab/Mpath

Manuscript: [Mpath maps multi-branching single-cell trajectories revealing progenitor cell progression during development](https://doi.org/10.1038/ncomms11988)

### SCORPIUS

Software: https://github.com/rcannood/SCORPIUS

Manuscript: [SCORPIUS improves trajectory inference and identifies novel modules in dendritic cell development](https://doi.org/10.1101/079509)

### SCODE

Software: https://github.com/hmatsu1226/SCODE

Manuscript: [SCODE: an efficient regulatory network inference algorithm from single-cell RNA-Seq during differentiation](https://doi.org/10.1093/bioinformatics/btx194)

### switchde

Software: https://bioconductor.org/packages/release/bioc/html/switchde.html

Manuscript: [switchde: inference of switch-like differential expression along single-cell trajectories](https://doi.org/10.1093/bioinformatics/btw798)

### MAGIC

Software: https://github.com/pkathail/magic/

Manuscript: [MAGIC: A diffusion-based imputation method reveals gene-gene interactions in single-cell RNA-sequencing data](https://doi.org/10.1101/111591)

### PHATE

Software: https://github.com/KrishnaswamyLab/PHATE

Manuscript: [Visualizing Transitions and Structure for High Dimensional Data Exploration
](https://doi.org/10.1101/120378)

### SOMSC

Manuscript: [SOMSC: Self-Organization-Map for High-Dimensional Single-Cell Data of Cellular States and Their Transitions](https://doi.org/10.1101/124693)

### TASIC

Software: https://www.andrew.cmu.edu/user/sabrinar/TASIC

Manuscript: [TASIC: determining branching models from time series single cell data](https://doi.org/10.1093/bioinformatics/btx173)

### FORKS

Software: https://github.com/macsharma/FORKS

Manuscript: [FORKS: Finding Orderings Robustly using K-means and Steiner trees](https://doi.org/10.1101/132811)

### UNCURL

Software: https://github.com/yjzhang/uncurl_python and https://github.com/mukhes3/UNCURL_release

Manuscript: [Scalable preprocessing for sparse scRNA-seq data exploiting prior knowledge
](https://doi.org/10.1101/142398)

### reCAT

Software: https://github.com/tinglab/reCAT

Manuscript: [Reconstructing cell cycle pseudo time-series via single-cell transcriptome data](https://doi.org/10.1038/s41467-017-00039-z)

### PhenoPath

Software: [Bioconductor package](https://doi.org/10.18129/B9.bioc.phenopath) and https://github.com/kieranrcampbell/phenopath

Manuscript: [Uncovering pseudotemporal trajectories with covariates from single cell and bulk expression data](https://doi.org/10.1038/s41467-018-04696-6)

### Branched Gaussian processes

Software: https://github.com/ManchesterBioinference/BranchedGP

Manuscript: [BGP: identifying gene-specific branching dynamics from single-cell data with a branching Gaussian process](https://doi.org/10.1186/s13059-018-1440-2)

### Branch-recombinant Gaussian Processes

Software: https://github.com/cap76/BranchingGPs

Manuscript: [Nonparametric Bayesian inference of transcriptional branching and recombination identifies regulators of early human germ cell development](https://doi.org/10.1101/167684)

### MATCHER

Software: https://github.com/jw156605/MATCHER and https://pypi.python.org/pypi/pymatcher

Manuscript: [MATCHER: manifold alignment reveals correspondence between single cell transcriptome and epigenome dynamics](https://doi.org/10.1186/s13059-017-1269-0)

### SoptSC

Software: https://github.com/WangShuxiong/SoptSC

Manuscript: [Low-rank Similarity Matrix Optimization Identifies Subpopulation Structure and Orders Single Cells in Pseudotime](https://doi.org/10.1101/168922)

### Di-SNE

Manuscript: [Assessment of clonal kinetics reveals multiple trajectories of dendritic cell development](https://doi.org/10.1101/167635)

### Population Balance Analysis

Software: https://github.com/AllonKleinLab/PBA

Manuscript: [Fundamental limits on dynamic inference from single cell snapshots](https://doi.org/10.1101/170118)

### Scanpy

Software: https://github.com/theislab/scanpy and https://pypi.python.org/pypi/scanpy

Manuscript: [Scanpy for analysis of large-scale single-cell gene expression data](https://doi.org/10.1101/174029)

### TIDES

Software: https://github.com/roshan9128/tides

Manuscript: [Learning Edge Rewiring in EMT from Single Cell Data
](https://doi.org/10.1101/155028)

### WADDINGTON-OT

Software: https://pypi.org/project/wot/

Manuscript: [Reconstruction of developmental landscapes by optimal-transport analysis of single-cell gene expression sheds light on cellular reprogramming](https://doi.org/10.1101/191056)

### pseudodynamics

Software: https://github.com/theislab/pseudodynamics

Manuscript: [Inferring population dynamics from single-cell RNA-sequencing time series data](https://doi.org/10.1038/s41587-019-0088-0)

### Partition-based graph abstraction

Software: https://github.com/theislab/paga

Manuscript: [PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells](https://doi.org/10.1186/s13059-019-1663-x)

### GPseudoRank

Software: https://github.com/magStra/GPseudoRank

Manuscript: [GPseudoRank: MCMC for sampling from posterior distributions of pseudo-orderings using Gaussian processes](https://doi.org/10.1093/bioinformatics/bty664)

### FateID

Software: https://github.com/dgrun/FateID

Manuscript: [FateID infers cell fate bias in multipotent progenitors from single-cell RNA-seq data](https://doi.org/10.1101/218115)

### cycleX

Manuscript: [cycleX: multi-dimensional pseudotime reveals cell cycle and differentiation relationship of dendritic cell progenitors](https://doi.org/10.1101/222372)

### GrandPrix

Software: https://github.com/ManchesterBioinference/GrandPrix

Manuscript: [GrandPrix: Scaling up the Bayesian GPLVM for single-cell data](https://doi.org/10.1101/227843)

### Partial differential equations

Manuscript: [Modeling acute myeloid leukemia in a continuum of differentiation states](https://doi.org/10.1101/237438)

### scdiff

Software: https://github.com/phoenixding/scdiff and https://pypi.python.org/pypi/scdiff/

Manuscript: [Reconstructing differentiation networks and their regulation from time series single cell expression data](https://doi.org/10.1101/gr.225979.117)

### Topographer

Manuscript: [Topographer Reveals Stochastic Dynamics of Cell Fate Decisions from Single-Cell RNA-Seq Data](https://doi.org/10.1101/251207)

### Markov-Chain Entropy

Manuscript: [Quantifying Waddington's epigenetic landscape: a comparison of single-cell potency measures](https://doi.org/10.1093/bib/bby093)

### Microstates

Manuscript: [Machine learning methods to reverse engineer dynamic gene regulatory networks governing cell state transitions](https://doi.org/10.1101/264671)

### DensityPath

Software: https://github.com/ucasdp/DensityPath

Manuscript: [DensityPath: a level-set algorithm to visualize and reconstruct cell developmental trajectories for large-scale single-cell RNAseq data](https://doi.org/10.1101/276311)

### STREAM

Software: https://github.com/pinellolab/stream

Manuscript: [STREAM: Single-cell Trajectories Reconstruction, Exploration And Mapping of omics data](https://doi.org/10.1101/302554)

Website: http://stream.pinellolab.partners.org/

### HopLand

Software: https://github.com/NetLand-NTU/HopLand

Manuscript: [HopLand: single-cell pseudotime recovery using continuous Hopfield network-based modeling of Waddington's epigenetic landscape](https://doi.org/10.1093/bioinformatics/btx232)

### Dynamic Distribution Decomposition

Manuscript: [Dynamic Distribution Decomposition for Single-Cell Snapshot Time Series Identifies Subpopulations and Trajectories during iPSC Reprogramming](https://doi.org/10.1101/367789)

### Continuous state HMMs 

Software: https://github.com/jessica1338/CSHMM-for-time-series-scRNA-Seq

Manuscript: [Continuous State HMMs for Modeling Time Series Single Cell RNA-Seq Data](https://doi.org/10.1101/380568)

### Palantir

Software: https://github.com/dpeerlab/Palantir/

Manuscript: [Palantir characterizes cell fate continuities in human hematopoiesis](https://doi.org/10.1038/s41587-019-0068-4)

### Trajectory inference Based on SNP information

Software: https://github.com/phoenixding/tbsp

Manuscript: [Cell lineage inference from SNP and scRNA-Seq data](https://doi.org/10.1101/401943)

### t-Distributed Gaussian Process Latent Variable Model

Software: https://github.com/architverma1/tGPLVM

Manuscript: [A robust nonlinear low-dimensional manifold for single cell RNA-seq data](https://doi.org/10.1101/443044)

### Sinova

Software: https://github.com/bionova/sinova

Manuscript: [Systematic Reconstruction of Molecular Cascades Regulating GP Development Using Single-Cell RNA-Seq](https://doi.org/10.1016/j.celrep.2016.04.043)

### Lineage tracing on transcriptional landscapes

Software: Multiple repositories

Manuscript: [Lineage tracing on transcriptional landscapes links state to fate during differentiation](https://doi.org/10.1101/467886)

### CALISTA

Software: https://github.com/CABSEL/CALISTA

Manuscript: [CALISTA: Clustering And Lineage Inference in Single-Cell Transcriptional Analysis](https://doi.org/10.1101/257550)

### TRACER

Manuscript: [Mapping Lung Cancer Epithelial-Mesenchymal Transition States and Trajectories with Single-Cell Resolution](https://doi.org/10.1101/570341)

### psupertime

Software: https://github.com/wmacnair/psupertime

Manuscript: [psupertime: supervised pseudotime inference for single cell RNA-seq data with sequential labels](https://doi.org/10.1101/622001)

### OscoNet

Software: https://github.com/alexisboukouvalas/OscoNet

Manuscript: [OscoNet: Inferring oscillatory gene networks](https://doi.org/10.1101/600049)

### Cyclum

Software: https://github.com/KChen-lab/cyclum

Manuscript: [Latent periodic process inference from single-cell RNA-seq data](https://doi.org/10.1101/625566)

### MERLoT

Software: https://github.com/soedinglab/merlot

Manuscript: [Reconstructing complex lineage trees from scRNA-seq data using MERLoT](https://doi.org/10.1093/nar/gkz706)

### scVelo

Software: https://scvelo.org

Manuscript: [Generalizing RNA velocity to transient cell states through dynamical modeling](https://doi.org/10.1038/s41587-020-0591-3)

### Tempora

Software: https://github.com/BaderLab/Tempora

Manuscript: [Tempora: cell trajectory inference using time-series single-cell RNA sequencing data](https://doi.org/10.1101/846907)

### CellCycleTRACER 

Software: https://ibm.biz/cellcycletracer-aas

Manuscript: [CellCycleTRACER accounts for cell cycle and volume in mass cytometry data](https://doi.org/10.1038/s41467-018-03005-5)

### TrajectoryNet

Software: https://github.com/krishnaswamylab/TrajectoryNet and https://github.com/KrishnaswamyLab/Cell-Dynamics-Pipeline

Manuscript: [TrajectoryNet: A Dynamic Optimal Transport Network for Modeling Cellular Dynamics](https://arxiv.org/abs/2002.04461)

Additional manuscript: [Learning transcriptional and regulatory dynamics driving cancer cell plasticity using neural ODE-based optimal transport](https://doi.org/10.1101/2023.03.28.534644)

### redPATH

Software: https://github.com/tinglab/redPATH

Manuscript: [redPATH: Reconstructing the Pseudo Development Time of Cell Lineages in Single-Cell RNA-Seq Data and Applications in Cancer](https://doi.org/10.1101/2020.03.05.977686)

### GraphDR and StructDR

Software: https://github.com/jzthree/quasildr

Manuscript: [An analytical framework for interpretable and generalizable 'quasilinear' single-cell data analysis](https://doi.org/10.1101/2020.04.12.022806)

### Pseudocell Tracer

Software: https://github.com/akds/pseudocell

Manuscript: [Inferring cellular trajectories from scRNA-seq using Pseudocell Tracer](https://doi.org/10.1101/2020.06.26.173179)

### TinGa

Software: https://github.com/Helena-todd/TinGa

Manuscript: [TinGa: fast and flexible trajectory inference with Growing Neural Gas ](https://doi.org/10.1093/bioinformatics/btaa463)

### scDEC

Software: https://github.com/kimmo1019/scDEC

Manuscript: [Simultaneous deep generative modeling and clustering of single cell genomic data](https://doi.org/10.1101/2020.08.17.254730)

### VeTra

Software: https://github.com/wgzgithub/VeTra

Manuscript: [VeTra: a new trajectory inference tool based on RNA velocity](https://doi.org/10.1101/2020.09.01.277095)

### DTFLOW

Software: https://github.com/statway/DTFLOW

Manuscript: [DTFLOW: Inference and Visualization of Single-cell Pseudo-temporal Trajectories Using Diffusion Propagation](https://doi.org/10.1101/2020.09.10.290973)

### CellPath

Software: https://github.com/PeterZZQ/CellPath

Manuscript: [Inference of multiple trajectories in single cell RNA-seq data from RNA velocity](https://doi.org/10.1101/2020.09.30.321125)

### CellRank / CellRank 2

Software: https://cellrank.org

CellRank manuscript: [CellRank for directed single-cell fate mapping](https://doi.org/10.1038/s41592-021-01346-6)

CellRank 2 manuscript: [Unified fate mapping in multiview single-cell data](https://doi.org/10.1101/2023.07.19.549685)

### Revelio

Software: https://github.com/danielschw188/Revelio

Manuscript: [The transcriptome dynamics of single cells during the cell cycle](https://doi.org/10.15252/msb.20209946)

### Cytopath

Software: https://github.com/aron0093/cytopath

Manuscript: [Cytopath: Simulation based inference of differentiation trajectories from RNA velocity fields](https://doi.org/10.1101/2020.12.21.423801)

### VIA

Software: https://github.com/ShobiStassen/VIA

Manuscript: [VIA: Generalized and scalable trajectory inference in single-cell omics data](https://doi.org/10.1101/2021.02.10.430705)

### Global Waddington-OT

Software: https://github.com/zsteve/gWOT

Manuscript: [Towards a mathematical theory of trajectory inference](https://arxiv.org/abs/2102.09204)

### StationaryOT

Software: https://github.com/zsteve/StationaryOT

Manuscript: [Optimal transport analysis reveals trajectories in steady-state systems](https://doi.org/10.1101/2021.03.02.433630)

### Condiments

Software: https://github.com/HectorRDB/condiments

Vignettes: https://hectorrdb.github.io/condimentsPaper/

Manuscript: [Trajectory inference across multiple conditions with condiments: differential topology, progression, differentiation, and expression](https://doi.org/10.1101/2021.03.09.433671)

### DeepCycle

Software: https://github.com/andreariba/DeepCycle

Manuscript: [Cell cycle gene regulation dynamics revealed by RNA velocity and deep-learning](https://doi.org/10.1101/2021.03.17.435887)

### Tricycle

Software: https://github.com/hansenlab/tricycle

Manuscript: [Universal prediction of cell cycle position using transfer learning](https://doi.org/10.1101/2021.04.06.438463)

### scShaper

Software: https://github.com/elolab/scshaper

Manuscript: [scShaper: ensemble method for fast and accurate linear trajectory inference from single-cell RNA-seq data](https://doi.org/10.1101/2021.05.03.442435)

### CCPE

Software: https://github.com/LiuJJ0327/CCPE

Manuscript: [CCPE: Cell Cycle Pseudotime Estimation for Single Cell RNA-seq Data](https://doi.org/10.1101/2021.06.13.448263)

### scDVF

Software: https://github.com/gersteinlab/scDVF

Manuscript: [scDVF: Data-driven Single-cell Transcriptomic Deep Velocity Field Learning with Neural Ordinary Differential Equations](https://doi.org/10.1101/2022.02.15.480564)

### Tempo

Software: https://github.com/bauerbach95/tempo

Manuscript: [Tempo: an unsupervised Bayesian algorithm for circadian phase inference in single-cell transcriptomics](https://doi.org/10.1101/2022.03.15.484454)

### DeepVelo

Software: https://github.com/bowang-lab/DeepVelo

Manuscript: [DeepVelo: Deep Learning extends RNA velocity to multi-lineage systems with cell-specific kinetics](https://doi.org/10.1101/2022.04.03.486877)

### LRT

Software: https://github.com/JuanXie19/LRT

Manuscript: [LRT: T Cell Trajectory Inference by Integrative Analysis of Single Cell TCR-seq and RNA-seq data](https://doi.org/10.1101/2022.04.14.488320)

### scTour

Software: https://github.com/LiQian-XC/sctour

Manuscript: [scTour: a deep learning architecture for robust inference and accurate prediction of cellular dynamics](https://doi.org/10.1101/2022.04.17.488600)

### UniTVelo

Software: https://github.com/StatBiomed/UniTVelo

Manuscript: [UniTVelo: temporally unified RNA velocity reinforces single-cell trajectory inference](https://doi.org/10.1101/2022.04.27.489808)

### VITAE

Software: https://github.com/jaydu1/VITAE

Manuscript: [Model-based Trajectory Inference for Single-Cell RNA Sequencing Using Deep Learning with a Mixture Prior](https://doi.org/10.1101/2020.12.26.424452)

### Real-time axis for T cells

Manuscript: [From pseudotime to true dynamics: reconstructing a real-time axis for T cells differentiation](https://doi.org/10.1101/2022.06.09.495431)

### GeneTrajectory

Software: https://github.com/KlugerLab/GeneTrajectory

Manuscript [Gene Trajectory Inference for Single-cell Data by Optimal Transport Metrics](https://doi.org/10.1101/2022.07.08.499404)

### scFates

Software: https://scfates.readthedocs.io/en/latest/

Manuscript: [scFates: a scalable python package for advanced pseudotime and bifurcation analysis from single cell data](https://doi.org/10.1101/2022.07.09.498657)

### veloVI

Software: https://github.com/YosefLab/velovi

Manuscript: [Deep generative modeling of transcriptional dynamics for RNA velocity analysis in single cells](https://doi.org/10.1101/2022.08.12.503709)

### MIRA

Software: https://github.com/cistrome/MIRA

Manuscript: [MIRA: joint regulatory modeling of multimodal expression and chromatin accessibility in single cells](https://doi.org/10.1038/s41592-022-01595-z)

### Pyro-Velocity

Software: https://github.com/pinellolab/pyrovelocity

Manuscript: [Pyro-Velocity: Probabilistic RNA Velocity inference from single-cell data](https://doi.org/10.1101/2022.09.12.507691)

### Totem

Software: https://github.com/elolab/Totem

Manuscript: [Cell-connectivity-guided trajectory inference from single-cell data](https://doi.org/10.1093/bioinformatics/btad515)

### scLTNN

Software: https://github.com/Starlitnightly/scltnn

Manuscript: [Identify the origin and end cells and infer the trajectory of cellular fate automatically](https://doi.org/10.1101/2022.09.28.510020)

### UPMM

Manuscript: [Modeling Single-Cell Dynamics Using Unbalanced Parameterized Monge Maps](https://doi.org/10.1101/2022.10.04.510766)

### PhyloVelo

Software: https://github.com/kunwang34/PhyloVelo

Manuscript: [Cell division history encodes directional information of fate transitions](https://doi.org/10.1101/2022.10.06.511094)

### MultiVelo

Software: https://github.com/welch-lab/MultiVelo

Manuscript: [Multi-omic single-cell velocity models epigenome–transcriptome interactions and improves cell fate prediction](https://doi.org/10.1038/s41587-022-01476-y)

### SCTC

Software: https://github.com/hailinphysics/sctc

Manuscript: [SCTC: inference of developmental potential from single-cell transcriptional complexity](https://doi.org/10.1101/2022.10.14.512265)

### MIOFlow

Software: https://github.com/KrishnaswamyLab/MIOFlow

Manuscript: [Manifold Interpolating Optimal-Transport Flows for Trajectory Inference](https://openreview.net/forum?id=ahAEhOtVif)

### DEAPLOG

Software: https://github.com/ZhangHongbo-Lab/DEAPLOG

Manuscript: [A method for differential expression analysis and pseudo-temporal locating and ordering of genes in single-cell transcriptomic data](https://doi.org/10.1101/2022.12.21.521359)

### GCSTI

Software: https://github.com/xznhy/rna-seq

Manuscript: [GCSTI: A Single-Cell Pseudotemporal Trajectory Inference Method Based on Graph Compression](https://doi.org/10.1109/TCBB.2023.3266109)

### moscot

Software: https://moscot-tools.org/

Manuscript: [Mapping cells through time and space with moscot](https://doi.org/10.1101/2023.05.11.540374)

### DELVE

Software: https://github.com/jranek/delve

Manuscript: [Feature selection for preserving biological trajectories in single-cell data](https://doi.org/10.1101/2023.05.09.540043)

### Velvet

Software: https://github.com/rorymaizels/velvet

Manuscript: [Deep dynamical modelling of developmental trajectories with temporal transcriptomics](https://doi.org/10.1101/2023.07.06.547989)

### ENTRAIN

Software: https://github.com/theimagelab/entrain

Manuscript: [ENTRAIN: integrating trajectory inference and gene regulatory networks with spatial data to co-localize the receptor-ligand interactions that specify cell fate](https://doi.org/10.1101/2023.07.09.548284)

### TFvelo

Software: https://github.com/xiaoyeye/TFvelo

Manuscript: [TFvelo: gene regulation inspired RNA velocity estimation](https://doi.org/10.1101/2023.07.12.548785)

### TopGen

Manuscript: [Unraveling cell differentiation mechanisms through topological exploration of single-cell developmental trajectories](https://doi.org/10.1101/2023.07.28.551057)

### cell2fate

Software: https://github.com/BayraktarLab/cell2fate

Manuscript: [Model-based inference of RNA velocity modules improves cell fate prediction](https://doi.org/10.1101/2023.08.03.551650)

### FBSDE

Software: https://github.com/Diebrate/population_model

Manuscript: [Modeling Single Cell Trajectory Using Forward-Backward Stochastic Differential Equations](https://doi.org/10.1101/2023.08.10.552373)

### cy2path

Software: https://github.com/aron0093/cy2path

Manuscript: [Factorial state-space modelling for kinetic clustering and lineage inference](https://doi.org/10.1101/2023.08.21.554135)

### scEGOT

Software: https://github.com/yachimura-lab/scEGOT

Manuscript: [scEGOT: Single-cell trajectory inference framework based on entropic Gaussian mixture optimal transport](https://doi.org/10.1101/2023.09.11.557102)

### Ricci flow

Manuscript: [Charting cellular differentiation trajectories with Ricci flow](https://doi.org/10.1101/2023.07.20.549833)

### Flat NB-VAE

Manuscript: [Modelling single-cell RNA-seq trajectories on a flat statistical manifold](https://openreview.net/pdf?id=sXRpvW3cRR)

### Sceptic

Software: https://github.com/Noble-Lab/Sceptic

Manuscript: [Pseudotime analysis for time-series single-cell sequencing and imaging data](https://doi.org/10.1101/2023.11.03.565575)

### VeloCycle

Software: https://github.com/lamanno-epfl/velocycle/

Manuscript: [Statistical inference with a manifold-constrained RNA velocity model uncovers cell cycle speed modulations](https://doi.org/10.1101/2024.01.18.576093)

### Chronocell

Software: https://github.com/pachterlab/FGP_2024

Manuscript: [Trajectory inference from single-cell genomics data with a process time model](https://doi.org/10.1101/2024.01.26.577510)

### PSEUDOTIMEABC

Software: https://github.com/keita-iida/PSEUDOTIMEABC

Manuscript: [Identifying key regulatory genes in drug resistance acquisition: Modeling pseudotime trajectories of single-cell transcriptome](https://doi.org/10.1101/2024.04.25.591115)

### scTEP

Software: https://cran.r-project.org/package=scTEP

Manuscript: [A robust and accurate single-cell data trajectory inference method using ensemble pseudotime](https://doi.org/10.1186/s12859-023-05179-2)

### TrajAtlas

Software: https://github.com/GilbertHan1011/TrajAtlas

Manuscript: [Trajectory-centric Framework TrajAtlas reveals multi-scale differentiation heterogeneity among cells, genes, and gene module in osteogenesis](https://doi.org/10.1101/2024.05.28.596174)













## Related topics

### Cicero

Manuscript: [Chromatin accessibility dynamics of myogenesis at single cell resolution](https://doi.org/10.1101/155473)

### Effects of imputation on cell ordering

Manuscript: [Comparison of computational methods for imputing single-cell RNA-sequencing data](https://doi.org/10.1101/241190)

### PROSSTT

Software: https://github.com/soedinglab/prosstt

Manuscript: [PROSSTT: probabilistic simulation of single-cell RNA-seq data for complex differentiation processes](https://doi.org/10.1101/256941)

### CellAlign

Software: https://github.com/shenorrLab/cellAlign

Manuscript: [Alignment of single-cell trajectories to compare cellular expression dynamics](https://doi.org/10.1038/nmeth.4628)

### CONFESS

Software: http://bioconductor.org/packages/release/bioc/html/CONFESS.html

Manuscript: [CONFESS: Fluorescence-based single-cell ordering in R](https://doi.org/10.1101/407932)

### Trajectory alignment

Software: https://www.cell.com/cms/10.1016/j.cels.2018.07.006/attachment/2b57ebff-a502-4819-b8ed-8a87d17a4ae7/mmc4.zip

Manuscript: [Aligning Single-Cell Developmental and Reprogramming Trajectories Identifies Molecular Determinants of Myogenic Reprogramming Outcome](https://doi.org/10.1016/j.cels.2018.07.006)

### ImageAEOT

Manuscript: [Autoencoder and Optimal Transport to Infer Single-Cell Trajectories of Biological Processes](https://doi.org/10.1101/455469)

### RNA velocity

Software: http://velocyto.org/

Manuscript: [RNA velocity of single cells](https://doi.org/10.1038/s41586-018-0414-6)

### devMap

Manuscript: [High Resolution Comparison of Cancer-Related Developmental Processes Using Trajectory Alignment](https://doi.org/10.1101/469601)

### Trajan

Software: https://github.com/canzarlab/Trajan

Manuscript: [Dynamic pseudo-time warping of complex single-cell trajectories](https://doi.org/10.1101/522672)

### SINGE

Software: https://github.com/gitter-lab/SINGE

Manuscript: [Network inference with Granger causality ensembles on single-cell transcriptomics](https://doi.org/10.1016/j.celrep.2022.110333)

### GPseudoClust

Software: https://github.com/magStra/GPseudoClust

Manuscript: [GPseudoClust: deconvolution of shared pseudo-trajectories at single-cell resolution](https://doi.org/10.1101/567115)

### tradeSeq

Software: http://www.bioconductor.org/packages/release/bioc/html/tradeSeq.html

Manuscript: [Trajectory-based differential expression analysis for single-cell sequencing data](https://doi.org/10.1038/s41467-020-14766-3)

### CORGI

Software: https://github.com/YutongWangUMich/corgi

Manuscript: [A gene filter for comparative analysis of single-cell RNA-sequencing trajectory datasets](https://doi.org/10.1101/637488)

### BEELINE

Software: https://github.com/murali-group/Beeline

Manuscript: [Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data](https://doi.org/10.1101/642926)

### Dynamo

Software: https://github.com/aristoteleo/dynamo-release

Manuscript: [Mapping Vector Field of Single Cells](https://doi.org/10.1101/696724)

### SERGIO

Software: https://github.com/PayamDiba/SERGIO

Manuscript: [A single-cell expression simulator guided by gene regulatory networks](https://doi.org/10.1101/716811)

### GeneSwitches

Software: https://geneswitches.ddnetbio.com/

Manuscript: [GeneSwitches : Ordering gene-expression and functional events in single-cell experiments](https://doi.org/10.1101/832626)

### CAPITAL

Software: https://github.com/ykat0/capital

Manuscript: [Alignment of time-course single-cell RNA-seq data with CAPITAL](https://doi.org/10.1101/859751)

### Transcriptional uncertainty landscapes

Manuscript: [Universality of cell differentiation trajectories revealed by a reconstruction of transcriptional uncertainty landscapes from single-cell transcriptomic data](https://doi.org/10.1101/2020.04.23.056069)

### Pseudo-Location

Manuscript: [Pseudo-Location: A novel predictor for predicting pseudo-temporal gene expression patterns using spatial functional regression](https://doi.org/10.1101/2020.06.11.145565)

### fishpond

Software: http://bioconductor.org/packages/devel/bioc/html/fishpond.html and https://github.com/skvanburen/scUncertaintyPaperCode

Manuscript: [Compression of quantification uncertainty for scRNA-seq counts](https://doi.org/10.1101/2020.07.06.189639)

### scHOT

Software: https://bioconductor.org/packages/scHOT

Manuscript: [Investigating higher-order interactions in single-cell data with scHOT](https://doi.org/10.1038/s41592-020-0885-x)

### PRESCIENT

Software: https://github.com/gifford-lab/prescient

Manuscript: [Generative modeling of single-cell population time series for inferring cell differentiation landscapes](https://doi.org/10.1101/2020.08.26.269332)

### Mathematics of RNA Velocity

Manuscript: [On the Mathematics of RNA Velocity I: Theoretical Analysis](https://doi.org/10.1101/2020.09.19.304584)

### Mathematics of RNA Velocity II

Manuscript: [On the Mathematics of RNA Velocity II: Algorithmic Aspects](https://doi.org/10.1101/2023.06.09.544270)

### PseudotimeDE

Software: https://github.com/SONGDONGYUAN1994/PseudotimeDE

Manuscript: [PseudotimeDE: inference of differential gene expression along cell pseudotime with well-calibrated p-values from single-cell RNA sequencing data](https://doi.org/10.1101/2020.11.17.387779)

### TIPS

Manuscript: [TIPS: Trajectory Inference of Pathway Significance through Pseudotime Comparison for Functional Assessment of single-cell RNAseq Data](https://doi.org/10.1101/2020.12.17.423360)

### VeloSim

Software: https://github.com/PeterZZQ/VeloSim

Manuscript: [VeloSim: Simulating single cell gene-expression and RNA velocity](https://doi.org/10.1101/2021.01.11.426277)

<details>
<summary>Abstract</summary>
The availability of high throughput single-cell RNA-Sequencing data allows researchers to study the molecular mechanisms that drive the temporal dynamics of cells during differentiation or development. Recent computational methods that build upon single-cell sequencing technology, such as trajectory inference or RNA-velocity estimation, provide a way for researchers to analyze the state of each cell during a continuous dynamic process. However, with the surge of such computational methods, there is still a lack of simulators that can model the cell temporal dynamics, and provide ground truth data to benchmark the computational methods.

Hereby we present VeloSim, a simulation software that can simulate the gene-expression kinetics in cells along continuous trajectories. VeloSim is able to take any trajectory structure composed of basic elements including “linear” and “cycle” as input, and outputs unspliced mRNA count matrix, spliced mRNA count matrix, cell pseudo-time and true RNA velocity of the cells. We demonstrate how VeloSim can be used to benchmark trajectory inference and RNA-velocity estimation methods with different amounts of biological and technical variation within the datasets. VeloSim is implemented into an R package available at https://github.com/PeterZZQ/VeloSim.
</details>

### SnapATAC

Software: https://github.com/r3fang/SnapATAC

Manuscript: [Comprehensive analysis of single cell ATAC-seq data with SnapATAC](https://doi.org/10.1038/s41467-021-21583-9)

### Spectral single cell

Software: https://github.com/mornitzan/spectral_sc

Manuscript: [Revealing lineage-related signals in single-cell gene expression using random matrix theory](https://doi.org/10.1073/pnas.1913931118)

### schubness

Software: https://github.com/EliseAld/schubness

Manuscript: [Hubness reduction improves clustering and trajectory inference in single-cell transcriptomic data](https://doi.org/10.1101/2021.03.18.435808)

### TreeVAE

Software: https://github.com/khalilouardini/treeVAE-reproducibility

Manuscript: [Reconstructing unobserved cellular states from paired single-cell lineage tracing and transcriptomics data](https://doi.org/10.1101/2021.05.28.446021)

### CoSpar

Software: https://cospar.readthedocs.io/

Manuscript: [Learning dynamics by computational integration of single cell genomic and lineage information](https://doi.org/10.1101/2021.05.06.443026)

### Lamian

Software: https://github.com/Winnie09/Lamian and https://github.com/Winnie09/trajectory_variability

Manuscript: [A statistical framework for differential pseudotime analysis with multiple single-cell RNA-seq samples](https://doi.org/10.1101/2021.07.10.451910)

### TedSim

Software: https://github.com/Galaxeee/TedSim

Manuscript: [TedSim: temporal dynamics simulation of single cell RNA-sequencing data and cell division history](https://doi.org/10.1101/2021.06.21.449283)

### Single-cell generalized trend model

Software: https://github.com/ElvisCuiHan/scGTM

Manuscript: [Single-cell generalized trend model (scGTM): a flexible and interpretable model of gene expression trend along cell pseudotime](https://doi.org/10.1093/bioinformatics/btac423)

### Expression and Velocity Integration

Software: https://github.com/jranek/EVI

Manuscript: [Integrating temporal single-cell gene expression modalities for trajectory inference and disease prediction](https://doi.org/10.1101/2022.03.01.482381)

### scSTEM

Software: https://github.com/alexQiSong/scSTEM

Manuscript: [scSTEM: clustering pseudotime ordered single-cell data](https://doi.org/10.1186/s13059-022-02716-9)

### SlowMoMan

Software: https://yunwilliamyu.github.io/SlowMoMan/

Manuscript: [SlowMoMan: A web app for discovery of important features along user-drawn trajectories in 2D embeddings](https://doi.org/10.1101/2022.08.23.505019)

### Dictys

Software: https://github.com/pinellolab/dictys

Manuscript: [Dictys: dynamic gene regulatory network dissects developmental continuum with single-cell multi-omics](https://doi.org/10.1101/2022.09.14.508036)

### LEAP

Software: https://cran.r-project.org/web/packages/LEAP/index.html

Manuscript: [LEAP: constructing gene co-expression networks for single-cell RNA-sequencing data using pseudotime ordering](https://doi.org/10.1093/bioinformatics/btw729)

### TrAGEDy

Software: https://github.com/No2Ross/TrAGEDy

Manuscript: [TrAGEDy: Trajectory Alignment of Gene Expression Dynamics](https://doi.org/10.1101/2022.12.21.521424)

### Genes2Genes

Software: https://github.com/Teichlab/Genes2Genes

Manuscript: [Gene-level alignment of single cell trajectories informs the progression of in vitro T cell differentiation](https://doi.org/10.1101/2023.03.08.531713)

### popInfer

Software: https://github.com/maclean-lab/popInfer

Manuscript: [Gene regulatory network inference with popInfer reveals dynamic regulation of hematopoietic stem cell quiescence upon diet restriction and aging](https://doi.org/10.1101/2023.04.18.537360)

### NeuroVelo

Manuscript: [NeuroVelo: interpretable learning of cellular dynamics from single-cell transcriptomic data](https://doi.org/10.1101/2023.11.17.567500)

### scLANE

Software: https://github.com/jr-leary7/scLANE

Manuscript: [Interpretable trajectory inference with single-cell Linear Adaptive Negative-binomial Expression (scLANE) testing](https://doi.org/10.1101/2023.12.19.572477)

### ExDyn

Software: https://github.com/kojikoji/exdyn

Manuscript: [Inferring extrinsic factor-dependent single-cell transcriptome dynamics using a deep generative model](https://doi.org/10.1101/2024.04.01.587302)

### Hodge decomposition

Software: https://github.com/WeilabMSU/HHD

Manuscript: [Hodge Decomposition of Single-Cell RNA Velocity](https://doi.org/10.1021/acs.jcim.4c00132)

### ConsensusVelo

Manuscript: [Quantifying uncertainty in RNA velocity](https://doi.org/10.1101/2024.05.14.594102)

### Mellon

Software: https://github.com/settylab/Mellon

Manuscript: [Quantifying cell-state densities in single-cell phenotypic landscapes using Mellon](https://doi.org/10.1038/s41592-024-02302-w)

### RNA velocity benchmark

Software: https://github.com/czbiohub-sf/comparison-RNAVelo (currently private or broken)

Manuscript: [Challenges and Progress in RNA Velocity: Comparative Analysis Across Multiple Biological Contexts](https://doi.org/10.1101/2024.06.25.600667)

### noSpliceVelo

Software: https://github.com/Tarun-Mahajan/noSpliceVelo

Manuscript: [noSpliceVelo infers gene expression dynamics without separating unspliced and spliced transcripts](https://doi.org/10.1101/2024.08.08.607261)
