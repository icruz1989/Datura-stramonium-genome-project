require(phangorn)
require(ape)
#### Set species
input.args <- "./SpeciesTree_rooted_node_labels.txt"
tree.rooted <- ape::root(read.tree(input.args), outgroup = "Petuniainflata", resolve.root = TRUE)
nto.na.ns.nt.dte.dti.cag.cam.st.sp.sly.spi.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("Nicotianatomentosiformis","Nicotianaattenuata","Nicotianasylvestris","Nicotianatabacum","DaturastramoniumTeo1","DaturastramoniumTic23","Capsicumannuumglabriusculum","Capsicumannuummorelia","Solanumtuberosum","Solanumpennelli","Solanumlycopersicum","Solanumpimpinellifolium")))
nto.na.ns.nt.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("Nicotianatomentosiformis","Nicotianaattenuata","Nicotianasylvestris","Nicotianatabacum")))
dte.dti.cag.cam.st.sp.sly.spi.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("DaturastramoniumTeo1","DaturastramoniumTic23","Capsicumannuumglabriusculum","Capsicumannuummorelia","Solanumtuberosum","Solanumpennelli","Solanumlycopersicum","Solanumpimpinellifolium")))
cag.cam.st.sp.sly.spi.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("Capsicumannuumglabriusculum","Capsicumannuummorelia","Solanumtuberosum","Solanumpennelli","Solanumlycopersicum","Solanumpimpinellifolium")))
cag.cam.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("Capsicumannunmglabriusculum","Capsicumannummorelia")))
dte.dti.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("DaturastramoniumTeo1","DaturastramoniumTic23")))
st.sp.sly.spi.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("Solanumtuberosum","Solanumpennelli","Solanumlycopersicum","Solanumpimpinellifolium")))
sp.sly.spi.mrca <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("Solanumpennelli","Solanumlycopersicum","Solanumpimpinellifolium")))
aet.split <- mrca.phylo(tree.rooted, which(tree.rooted$tip.label %in% c("Petuniainflata","Solanumlycopersicum")))
chr.df <- data.frame(node = c(nto.na.ns.nt.dte.dti.cag.cam.st.sp.sly.spi.mrca, nto.na.ns.nt.mrca, dte.dti.cag.cam.st.sp.sly.spi.mrca, cag.cam.st.sp.sly.spi.mrca, cag.cam.mrca, dte.dti.mrca, st.sp.sly.spi.mrca, sp.sly.spi.mrca, aet.split), age.min = c(31, 10, 30, 19, 0.01, 7.9, 1.5, 35), stringsAsFactors = FALSE)

### Time Tree
chr.df$age.max <- chr.df$age.min
chr.df[[which(chr.df$node == aet.split), "age.max"]] <- 35
chr.df$soft.bounds <- FALSE
ctrl <- chronos.control(nb.rate.cat = 1)
tree.rooted.chrono <- chronos(tree.rooted, calibration = chr.df, lambda = 3.2)

#' Scale the rooted tree
ctrl <- chronos.control(nb.rate.cat = 1)
tree.rooted.chrono <- chronos(tree.rooted, calibration = chr.df, lambda = 3.2)  # Lambda value obtained from references
#' Round the edge lengths:
tree.rooted.chrono$edge.length <- round(tree.rooted.chrono$edge.length, digits = 1)
#' Save the results
write.tree(tree.rooted.chrono, "./allsol_forCAFE_concate_ml_gamma_rescaled_ultrametric.newick")

#' Plot the tree:
pdf( "./all_sol_forCAFE_concate_ml_gamma_rescaled_ultrametric3.pdf" )
plot.phylo(tree.rooted.chrono, show.node.label=TRUE)
edgelabels(tree.rooted.chrono$edge.length)
axisPhylo()
dev.off()

#### paste nodes to tree (arbol sin nodos)
fam.tree <- read.tree(file="allsol__corrected_concate_ml_gamma_rescaled_ultrametric_from_cluster.newick")
fam.tree.no.node.labels <- fam.tree
fam.tree.no.node.labels$node.label <- NULL
fam.tree.no.node.labels.path <- file.path(fam.dir, paste(fam.name, "allsol__corrected_concate_ml_gamma_rescaled_ultrametric_from_cluster.newick", sep = ""))

fam.dir <- "/Users/mijaildelacruz/Desktop/Rorthofinder"











