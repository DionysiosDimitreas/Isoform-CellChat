# IsoCellChat: Isoform-Resolved Cell-Cell Communication

**IsoCellChat** is an R package designed to perform isoform-aware cell-cell communication analysis. Unlike standard tools that aggregate expression to the gene level, IsoCellChat retains transcript-level resolution, applies protein localization filtering (DeepLoc), and uses a domain-aware database to infer specific isoform-isoform interactions.

## Installation

You can install the development version of IsoCellChat from GitHub:

```r
# install.packages("devtools")
devtools::install_github("DionysiosDimitreas/IsoCellChat")
