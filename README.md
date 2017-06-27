# Single-cell RNA-seq pseudotime estimation algorithms

Single cells, many algorithms.  The goal of this page is to catalog the many
algorithms that estimate pseudotimes for cells based on their gene expression
levels.  This problem is also referred to as single-cell trajectory inference
or ordering.  Ultimately, it will contain method names, software links, manuscript
links, and brief comments on the particular strengths of each method.
Initially, it seeks simply to list as many methods as possible.  Some related
methods not specifically designed for RNA-seq (e.g. mass cytometry) are
included as well.  The list also includes methods
that are specifically designed to take pseudotemporal data as input.

The initial list was created by Anthony Gitter, but pull requests are very welcome.


## Problem overview

Informally, the pseudotime estimation problem can be stated as:
- **Given:** single-cell gene expression measurements for a heterogeneous
collection of cells that is transitioning from biological state **A** to
state **B**
- **Return:** a quantitative value for each cell that represents its
progress in the **A** to **B** transition

There are many ways to approach this problem, and major algorithmic steps that
are common to most (but not all) methods are:
- Reduce the dimension of the dataset
- Find a smooth progression through the low dimensional data, assuming
that cells that are nearby one another in the low dimensional space have similar
expression levels because they are at similar points in to **A** to **B**
process

Dimension reduction sometimes relies on knowledge of important marker genes and
sometimes uses the full gene expression matrix.  The trajectory through
the low dimensional space can be identified using graph algorithms (e.g.,
minimum spanning tree or shortest path), principal curves, or probabilistic
models (e.g., Gaussian process).

https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0927-y
http://genome.cshlp.org/content/25/10/1491.full and
http://www.nature.com/nature/journal/v541/n7637/full/nature21350.html
provide a more comprehensive
overview of single-cell RNA-seq and the pseudotime estimation problem.
http://onlinelibrary.wiley.com/wol1/doi/10.1002/eji.201646347/abstract
reviews pseudotime inference algorithms.
OMICtools has compiled a more general list of methods for single-cell RNA-seq
https://omictools.com/single-cell-rna-seq-category.


## Choosing a method

Many algorithms compare themselves to Monocle, a pioneering method in this area,
but do not evaluate Monocle 2 or other recent approaches.  There is not yet a
systematic benchmark.  Qualitatively, some of the distinguishing factors among
algorithms include:
- Use of prior knowledge such as capture times (DeLorean) or switch-like marker
genes (Ouija)
- Modeling specific types of biological processes such as branching processes
in differentiation (multiple methods) or cyclic processes (Oscope)
- Return a single pseudotime or a posterior distribution over pseudotimes
for each cell


## Algorithms

### Monocle / Monocle 2 / Census
Software: https://bioconductor.org/packages/release/bioc/html/monocle.html

Monocle manuscript: http://www.nature.com/nbt/journal/v32/n4/full/nbt.2859.html

Census manuscript: http://doi.org/10.1038/nmeth.4150

Monocle 2 manuscript: https://doi.org/10.1101/110668

### Wanderlust / Cycler / Wishbone
Wanderlust software: http://www.c2b2.columbia.edu/danapeerlab/html/wanderlust.html

Wanderlust manuscript: http://www.cell.com/cell/abstract/S0092-8674(14)00471-1

Cycler software: http://www.cellcycler.org/

Cycler manuscript: http://www.nature.com/nmeth/journal/v12/n10/full/nmeth.3545.html

Wishbone software: http://www.c2b2.columbia.edu/danapeerlab/html/wishbone.html

Wishbone manuscript: http://www.nature.com/nbt/journal/v34/n6/full/nbt.3569.html

### SCUBA
Software: https://github.com/gcyuan/SCUBA

Manuscript: http://www.pnas.org/content/111/52/E5643.abstract

### Oscope
Software: https://www.biostat.wisc.edu/~kendzior/OSCOPE/

Manuscript: http://www.nature.com/nmeth/journal/v12/n10/full/nmeth.3549.html

### Diffusion maps / destiny
Diffusion maps software: https://www.helmholtz-muenchen.de/icb/research/groups/machine-learning/projects/single-cell-diffusion-map/index.html

Diffusion maps manuscripts: http://www.nature.com/nbt/journal/v33/n3/full/nbt.3154.html and http://bioinformatics.oxfordjournals.org/content/31/18/2989.abstract

destiny software: http://www.helmholtz-muenchen.de/icb/research/groups/quantitative-single-cell-dynamics/software/destiny/index.html

destiny manuscript: http://bioinformatics.oxfordjournals.org/content/32/8/1241.full

### DeLorean
Software: https://github.com/JohnReid/DeLorean

Manuscript: http://bioinformatics.oxfordjournals.org/content/early/2016/06/16/bioinformatics.btw372.abstract

### Waterfall
Manuscript: http://www.cell.com/cell-stem-cell/abstract/S1934-5909(15)00312-4

### Embeddr
Software: https://github.com/kieranrcampbell/embeddr

Manuscript: http://biorxiv.org/content/early/2015/09/18/027219

### GP-LVM and pseudogp
GP-LVM software: https://github.com/kieranrcampbell/gpseudotime

GP-LVM manuscript: http://biorxiv.org/content/early/2015/09/15/026872

pseudogp software: https://github.com/kieranrcampbell/pseudogp

pseudogp manuscript: http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005212

### GP-LVM
Analysis code: https://github.com/Teichlab/spectrum-of-differentiation-supplements

Manuscript: http://www.cell.com/cell-reports/abstract/S2211-1247(15)01538-7

### SLICER
Software: https://github.com/jw156605/SLICER

Manuscript: http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0975-3

### TSCAN
Software: https://github.com/zji90/TSCAN

Manuscript: http://nar.oxfordjournals.org/content/early/2016/05/13/nar.gkw430.abstract

### SCOUP
Software: https://github.com/hmatsu1226/SCOUP

Manuscript: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1109-3

### Topslam
Software: https://github.com/mzwiessele/topslam

Manuscript: http://biorxiv.org/content/early/2016/06/08/057778

### Ouija
Software: https://github.com/kieranrcampbell/ouija

Manuscript: http://biorxiv.org/content/early/2016/06/23/060442

### Slingshot
Software: https://github.com/kstreet13/slingshot

Extended vignette: https://github.com/drisso/bioc2016singlecell/tree/master/vignettes

Manuscript: http://biorxiv.org/content/early/2017/04/19/128843

### GPfates

Software: https://github.com/Teichlab/GPfates

Manuscript: http://biorxiv.org/content/early/2016/09/13/074971

### SCIMITAR

Software: https://github.com/dimenwarper/scimitar

Manuscript: http://biorxiv.org/content/early/2016/08/18/070151

### WaveCrest

Software: https://github.com/lengning/WaveCrest

Manuscript: http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1033-x

### LEAP

Software: https://cran.r-project.org/web/packages/LEAP/index.html

### CellTree

Software: http://bioconductor.org/packages/release/bioc/html/cellTree.html

Manuscript: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1175-6

### Bayesian hierarchical mixture of factor analysers (MFA)

Software: http://www.github.com/kieranrcampbell/mfa

Manuscript: http://doi.org/10.12688/wellcomeopenres.11087.1

Preprint: http://biorxiv.org/content/early/2016/09/21/076547

### Mpath

Software: https://github.com/JinmiaoChenLab/Mpath

Manuscript: http://www.nature.com/ncomms/2016/160630/ncomms11988/full/ncomms11988.html

### SCORPIUS

Software: https://github.com/rcannood/SCORPIUS

Manuscript: http://biorxiv.org/content/early/2016/10/06/079509

### SCODE

Software: https://github.com/hmatsu1226/SCODE

Manuscript: http://biorxiv.org/content/early/2016/11/22/088856

### switchde

Software: https://bioconductor.org/packages/release/bioc/html/switchde.html

Manuscript: https://doi.org/10.1093/bioinformatics/btw798

### MAGIC

Software: https://github.com/pkathail/magic/

Manuscript: https://doi.org/10.1101/111591

### PHATE

Software: https://github.com/SmitaKrishnaswamy/PHATE

Manuscript: https://doi.org/10.1101/120378

### SOMSC

Manuscript: https://doi.org/10.1101/124693

### TASIC

Software: https://www.andrew.cmu.edu/user/sabrinar/TASIC

Manuscript: https://doi.org/10.1093/bioinformatics/btx173

### FORKS

Software: https://github.com/macsharma/FORKS

Manuscript: https://doi.org/10.1101/132811

### UNCURL

Software: https://github.com/mukhes3/UNCURL_release

Manuscript: https://doi.org/10.1101/142398

### reCAT

Software: https://github.com/tinglab/reCAT

Manuscript: https://doi.org/10.1038/s41467-017-00039-z

## Related topics

### Cicero

Manuscript: https://doi.org/10.1101/155473
