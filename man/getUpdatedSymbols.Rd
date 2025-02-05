% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getUpdatedSymbols.R
\name{getUpdatedSymbols}
\alias{getUpdatedSymbols}
\title{Get rid of your outdated gene symbols}
\usage{
getUpdatedSymbols(geneSymbols)
}
\arguments{
\item{geneSymbols}{input gene symbols to update}
}
\value{
a data.frame with two variables: the original and updated symbols
}
\description{
Update previous gene symbols to up-to-date approved gene symbols.
}
\details{
This is some description of this function

Gene symbols are updated accoding to "hgnc_20210601.txt" file download from HGNC database (custom downloads).
Previous symbols are updated, approved symbols are preserved.
}
\note{
Alias symbols will not be updated to reduce multi-mapping which is confused.

如果Previous symbols还是匹配到多个Approved symbols, 或者多个Previous symbols 匹配到Approved symbols, 则按照Approved symbols排序后优先保留前面的
}
