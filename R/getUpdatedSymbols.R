#' This is some description of this function
#'
#' @title Get rid of your outdated gene symbols
#'
#' @description Update previous gene symbols to up-to-date approved gene symbols.
#'
#' @details Gene symbols are updated accoding to "hgnc_20210601.txt" file download from HGNC database (custom downloads).
#' Previous symbols are updated, approved symbols are preserved.
#'
#' @note Alias symbols will not be updated to reduce multi-mapping which is confused.
#'
#' 如果Previous symbols还是匹配到多个Approved symbols, 或者多个Previous symbols 匹配到Approved symbols, 则按照Approved symbols排序后优先保留前面的
#'
#' @param geneSymbols input gene symbols to update
#'
#' @return a data.frame with two variables: the original and updated symbols


getUpdatedSymbols <- function(geneSymbols) {
    require(tidyr)
    require(dplyr)
    hgnc <- read.table(system.file("extdata", "hgnc_20210601.txt", package = "MyTools"),
                       sep = "\t", header = T, comment.char = "", quote = "")
    hgnc.previous <- hgnc %>% separate_rows(Previous.symbols, sep = ", ")
    hgnc.previous <- hgnc.previous %>% filter(Previous.symbols != "")

    approved <- geneSymbols[geneSymbols %in% hgnc$Approved.symbol] # approved symbol

    not.approved <- setdiff(geneSymbols, approved) # not approved symbol
    previous <- hgnc.previous %>%
        filter(Previous.symbols %in% not.approved) %>% # update previous symbol in not approved symbols
        filter(!(Approved.symbol %in% approved)) # remove symbols which mismatch to already approved ones
    previous <- previous[order(previous$Approved.symbol, previous$Previous.symbols), ]
    previous <- previous %>%
        filter(!duplicated(Approved.symbol)) %>%
        filter(!duplicated(Previous.symbols)) # filter multi-map

    # combine approved and previous symbol
    approved <- data.frame(original = approved, updated = approved)
    previous <- data.frame(original = previous$Previous.symbols, updated = previous$Approved.symbol)
    trans.table <- rbind(approved, previous)
    # add deleted symbols
    deleted <- setdiff(geneSymbols, trans.table$original)
    deleted <- data.frame(original = deleted, updated = rep("", length(deleted)))
    trans.table <- rbind(trans.table, deleted)
    # return data frame
    rownames(trans.table) <- trans.table$original
    trans.table <- trans.table[geneSymbols, ]
}



