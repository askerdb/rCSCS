### Install

Install with:

 `devtools::install_github("askerdb/rCSCS")`

### Example

Example data with fecal samples of Crohns diseased children (CD) and Healthy children (HC), run in both positive and negative ionisation mode

(requires the packages rCSCS, ape, and ggplot2)
```R
#Download the GNPS job
cicra <- rCSCS::prepare_GNPS("b0524246804a4b50a8a4ec6244a8be2e", "cicra_example") 
#Run CSCS on the relevant samples (The ones not labelled full scan or standard)
cicra_cscs <- rCSCS::cscs(cicra$features[-grep("FS|Standard", colnames(cicra$features))], cicra$css)
#Run PCoA
cicra_pcoa <- ape::pcoa(cicra_cscs)
#Extract the disease status
cicra_meta = unlist(lapply(strsplit(rownames(cicra_pcoa$vectors), "_" ), function(x) x[[1]]))
#Plot the data
ggplot2::ggplot(data.frame(cicra_pcoa$vectors[,1:2], Meta = cicra_meta, Name = row.names(cicra_pcoa$vectors)) ,
               ggplot2::aes(x = Axis.1, y = Axis.2, color = Meta)) + ggplot2::geom_point()  +
  ggplot2::geom_label(ggplot2::aes(label = Name))

```
