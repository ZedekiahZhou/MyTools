#' @include myplot.R

# ==========================================================================
# Author:   Aaron
# Aim: GO and KEGG analysis, with plot and tsv output
# Version:  1.0
# Date:     Nov 4, 2019
# ==========================================================================

require(clusterProfiler)
require(DOSE)
require(GO.db)
require(ggplot2)


#' GO and KEGG analysis, with plot and tsv output
#'
#' @param genelist Gene list to perform enrichment analysis on
#' @param IDtype GeneID type: c("ENSEMBL", "ENTREZID", "SYMBOL", "UNIPROT"), default is "SYMBOL"
#' @param species c("hs" or "mm"), default is "hs"
#' @param ont ontology of GO: c("BP", "MF", "CC", "ALL"), default is "BP"
#' @param outdir directory to place the enrichment *.tsv file and plot file, default is the work directory
#' @param data.name Gene list name, serving as prefix of output file, default is NULL
#' @param color plot color pass to myGOhist()
#' @param n number of term to plot, pass to myGOhist()
#'
#' @details Gene list will first converted using bitr(), then perform GO and KEGG pathway enrichment, output the enrichment terms
#' in .tsv file and hist plot in .pdf
#'

GOKEGG <- function(genelist, IDtype = "SYMBOL", species = "hs",
                   ont = "BP",
                   outdir = "./", data.name = NULL, color = NULL,
                   n = 10) {
    # check species
    if (species == "hs") {
        require(org.Hs.eg.db)
        species1 <- "hsa"
        db <- "org.Hs.eg.db"
    } else if (species == "mm") {
        require(org.Mm.eg.db)
        species1 <- "mmu"
        db <- "org.Mm.eg.db"
    } else stop("Unrecognized species!")

    # check IDtype
    if (IDtype == "ENSEMBL") {
        genelist <- sub("\\.[0-9]+", "", genelist)
        readable <-  TRUE
    }

    # convert gene IDs
    gene <- bitr(genelist, fromType=IDtype, OrgDb=db,
                 toType=c("UNIPROT", "SYMBOL", "ENSEMBL", "ENTREZID"))


    # GO BP
    egoBP <- enrichGO(gene         = unique(gene$ENTREZID),
                      OrgDb         = db,
                      keyType       = 'ENTREZID',
                      ont           = ont,
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      pool = T, readable = T)

    egoBP <- data.frame(egoBP)
    if (nrow(egoBP)==0) {
        warning("No GO BP terms enriched!")
    } else {
        write.table(egoBP, file = paste0(outdir, "/", data.name, ".GO.tsv"),
                    quote = F, sep = "\t", row.names = F)
        p <- myGOhist(egoBP, n = n, color = color)
        pdf(file = paste0(outdir, "/", data.name, ".GO.pdf"), width = 7,
            height = 0.2*min(n, nrow(egoBP)) + 1)
        print(p)
        dev.off()
    }

    # KEGG
    kk <- enrichKEGG(gene         = unique(gene$UNIPROT),
                     keyType = 'uniprot',
                     organism     = species1,
                     qvalueCutoff = 0.05)
    kk <- data.frame(kk)
    if (nrow(kk)==0) {
        warning("No KEGG terms enriched!")
    } else {
        write.table(kk, file = paste0(outdir, "/", data.name, ".KEGG.tsv"),
                    quote = F, sep = "\t", row.names = F)
        p <- myGOhist(kk, n = n, color = color)
        pdf(file = paste0(outdir, "/", data.name, ".KEGG.pdf"), width = 7,
            height = 0.2*min(n, nrow(kk)) + 1)
        print(p)
        dev.off()
    }
}
