require(phangorn)

message("USAGE: Rscript generate_ultrametric_tree.R")

input.args <- "./spe8_ortho_concate_ml_gamma_rescaled.newick"

tree.rooted <- ape::root(read.tree(input.args), outgroup = "A.arabicum", resolve.root = TRUE)

#' We know the divergence times between the following clades based on fossil
#' data and other estimates [see references at bottom of this script]:
aly.ath.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("A.thaliana", 
    "A.lyrata")))
chi.aly.ath.cru.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% 
    c("A.lyrata", "A.thaliana", "C.rubella", "C.hirsuta")))
bra.esa.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("B.rapa", 
    "E.salsugineum")))
ingroup.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("A.thaliana", 
    "E.salsugineum")))
aet.split <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("A.arabicum", 
    "A.thaliana")))
chr.df <- data.frame(node = c(aly.ath.mrca, chi.aly.ath.cru.mrca, bra.esa.mrca, ingroup.mrca, 
    aet.split), age.min = c(13, 35.6, 38.4, 43.2, 45), stringsAsFactors = FALSE)
chr.df$age.max <- chr.df$age.min
chr.df[[which(chr.df$node == aet.split), "age.max"]] <- 60
chr.df$soft.bounds <- FALSE

#' Scale the rooted tree
ctrl <- chronos.control(nb.rate.cat = 1)
tree.rooted.chrono <- chronos(tree.rooted, calibration = chr.df, lambda = 3.2)  # Lambda value obtained from references
#' Round the edge lengths:
tree.rooted.chrono$edge.length <- round(tree.rooted.chrono$edge.length, digits = 1)
#' Save the results
write.tree(tree.rooted.chrono, "./spe8_ortho_concate_ml_gamma_rescaled_ultrametric.newick")

#' Plot the tree:
pdf( "./spe8_ortho_concate_ml_gamma_rescaled_ultrametric.pdf" )
plot.phylo(tree.rooted.chrono, show.node.label=TRUE)
edgelabels(tree.rooted.chrono$edge.length)
axisPhylo()
dev.off()

png( "./spe8_ortho_concate_ml_gamma_rescaled_ultrametric.png" )
plot.phylo(tree.rooted.chrono, show.node.label=TRUE)
edgelabels(tree.rooted.chrono$edge.length)
axisPhylo()
dev.off()

#' REFERENCES:
#' [1] http://www.mobot.org/mobot/research/apweb/orders/brassicalesweb.htm
#' [2] Yang, Ruolin, David J. Jarvis, Hao Chen, Mark Beilstein, Jane Grimwood,
#'     Jerry Jenkins, ShengQiang Shu, et al. “The Reference Genome of the
#'     Halophytic Plant Eutrema Salsugineum.” Plant Genetics and Genomics 4
#'     (2013): 46. doi:10.3389/fpls.2013.00046.                                                                                                            
