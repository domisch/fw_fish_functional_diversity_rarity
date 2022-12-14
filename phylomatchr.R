# Original function from the PDcalc package, available at https://github.com/davidnipperess/PDcalc/tree/master/R

#' Matches taxa in a set of sites with the tips of a phylogenetic tree
#'
#' @param x is the community data given as a \code{data.frame} or \code{matrix} 
#'   with species/OTUs as columns and samples/sites as rows (like in the 
#'   \code{vegan} package). Columns are labelled with the names of the 
#'   species/OTUs. Rows are labelled with the names of the samples/sites. Data 
#'   can be either abundance or incidence (0/1). Column labels must match tip 
#'   labels in the phylogenetic tree exactly!
#' @param phy is a rooted phylogenetic tree with branch lengths stored as a 
#'   phylo object (as in the \code{ape} package) with terminal nodes labelled
#'   with names matching those of the community data table. Note that the
#'   function trims away any terminal taxa not present in the community data
#'   table, so it is not necessary to do this beforehand.
#'
#' @return A list of two vectors. The first vector is the taxa found in the tree
#'   but not in the sites. The second vector is the taxa found in the sites but
#'   not in the tree.
#' @export
#' @examples
phylomatchr <- function(x,phy) {
  
  missing_in_x <- setdiff(phy$tip.label, colnames(x))
  missing_in_phy <- setdiff(colnames(x), phy$tip.label)
  
  return(list("taxa.missing.from.x"=missing_in_x,
              "taxa.missing.from.phy"=missing_in_phy))
  
}