# Single-cell RNA-seq pseudotime estimation algorithms
[![Build Status](https://travis-ci.org/agitter/single-cell-pseudotime.svg?branch=master)](https://travis-ci.org/agitter/single-cell-pseudotime)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1297422.svg)](https://doi.org/10.5281/zenodo.1297422)

Single cells, many algorithms.
The goal of this page is to catalog the many algorithms that estimate pseudotimes for cells based on their gene expression levels.
This problem is also referred to as single-cell trajectory inference
or ordering.
Ultimately, it will contain method names, software links, manuscript links, and brief comments on the particular strengths of each method.
Initially, it seeks simply to list as many methods as possible.
Some related methods not specifically designed for RNA-seq (e.g. mass cytometry) are included as well.
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

Dimension reduction sometimes relies on knowledge of important marker genes and sometimes uses the full gene expression matrix.  The trajectory through the low dimensional space can be identified using graph algorithms (e.g., minimum spanning tree or shortest path), principal curves, or probabilistic models (e.g., Gaussian process).

[Bacher and Kendziorski 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0927-y), [Trapnell 2015](http://genome.cshlp.org/content/25/10/1491.full), [Tanay and Regev 2017](http://www.nature.com/nature/journal/v541/n7637/full/nature21350.html) and [Moon et al. 2017](https://doi.org/10.1016/j.coisb.2017.12.008) provide a more comprehensive overview of single-cell RNA-seq and the pseudotime estimation problem.
[Cannoodt et al. 2016](http://onlinelibrary.wiley.com/wol1/doi/10.1002/eji.201646347/abstract) reviews pseudotime inference algorithms.
[Pablo Cordero's blog post](http://hyperparameter.space/blog/a-single-cell-journey-from-mechanistic-to-descriptive-modeling-and-back-again/) discusses why it is hard to recover the correct dynamics of a system from single-cell data.
For more general lists of methods for single-cell RNA-seq see [OMICtools](https://omictools.com/single-cell-rna-seq-category), [seandavi/awesome-single-cell](https://github.com/seandavi/awesome-single-cell), and [scRNA-tools](https://www.scrna-tools.org/).
The Hemberg lab [single-cell RNA-seq course](https://hemberg-lab.github.io/scRNA.seq.course/index.html) has an [overview of five pseudotime algorithms](https://hemberg-lab.github.io/scRNA.seq.course/biological-analysis.html#pseudotime-analysis) with usage examples.

Single-cell expression data have also inspired new methods for gene regulatory network reconstruction, as reviewed by [Fiers et al. 2018](https://doi.org/10.1093/bfgp/elx046) and [Todorov et al. 2018](https://doi.org/10.1007/978-1-4939-8882-2_10).
Several of these, such as [SCINGE](https://doi.org/10.1101/534834), treat pseudotime annotations as time points and extend traditional time series network inference algorithms for single-cell data.

## Choosing a method

Some of the distinguishing factors among algorithms include:
- Use of prior knowledge such as capture times (DeLorean) or switch-like marker genes (Ouija)
- Modeling specific types of biological processes such as branching processes in differentiation (multiple methods) or cyclic processes (Oscope)
- Return a single pseudotime or a posterior distribution over pseudotimes for each cell
- Perform additional analyses after inferring pseudotimes such as regulatory network inference or identifying differentially expressed genes over pseudotime

[Saelens et al. 2018](https://doi.org/10.1101/276907) performed a comprehensive evaluation of 29 different single-cell trajectory inference methods and discuss the different types of algorithms in more detail.
They benchmark both quantitative performance and assess software quality.
See their [GitHub repository](https://github.com/dynverse/dynverse) as well.
[Tian et al. 2018](https://doi.org/10.1101/433102) also include trajectory inference algorithms in their single-cell RNA-seq benchmarking study.

## Algorithms

### Monocle / Monocle 2 / Census
Software: https://bioconductor.org/packages/release/bioc/html/monocle.html

Monocle manuscript: [The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells](https://doi.org/10.1038/nbt.2859)

Census manuscript: [Single-cell mRNA quantification and differential analysis with Census](https://doi.org/10.1038/nmeth.4150)

Monocle 2 manuscript: [Reversed graph embedding resolves complex single-cell trajectories](https://doi.org/10.1038/nmeth.4402)

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

### LEAP

Software: https://cran.r-project.org/web/packages/LEAP/index.html

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

Manuscript: [Beyond pseudotime: Following T-cell maturation in single-cell RNAseq time series](https://doi.org/10.1101/219188)

### Approximate graph abstraction

Software: https://github.com/theislab/graph_abstraction

Manuscript: [Graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells](https://doi.org/10.1101/208819)

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

Manuscript: [Quantifying Waddington's epigenetic landscape: a comparison of single-cell potency measures](https://doi.org/10.1101/257220)

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

Manuscript: [Palantir characterizes cell fate continuities in human hematopoiesis](https://doi.org/10.1101/385328)

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

### SCINGE

Software: https://github.com/gitter-lab/SCINGE

Manuscript: [Network Inference with Granger Causality Ensembles on Single-Cell Transcriptomic Data](https://doi.org/10.1101/534834)
