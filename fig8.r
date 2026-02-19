rm(list = ls())
library(VennDiagram)
library(igraph)
library(grid)
target_gene <- "ZBTB16"

generate_mirna_names <- function(n) {
  numbers <- sample(1:2000, n)
  suffixes <- c("-3p", "-5p", "", "a", "b")
  paste0("hsa-miR-", numbers, sample(suffixes, n, replace = TRUE))
}

#making pool of mirnas
all_mirnas <- generate_mirna_names(500)

set.seed(456)
#making 3 datasets
db_targetscan <- sample(all_mirnas, 60)
db_miranda    <- sample(all_mirnas, 55)
db_diana      <- sample(all_mirnas, 50)

#venn diagram
futile.logger::flog.threshold(futile.logger::ERROR)
venn_list <- list(TargetScan = db_targetscan, Miranda = db_miranda, DIANA = db_diana)

venn_plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  fill = c("#E69F00", "#56B4E9", "#009E73"),
  alpha = 0.6,
  main = paste("miRNA Overlap for", target_gene)
)

grid.newpage()
grid.draw(venn_plot)

#finding intersections
common_nodes <- unique(c(
  intersect(db_targetscan, db_miranda),
  intersect(db_miranda, db_diana),
  intersect(db_targetscan, db_diana)
))

edges <- data.frame(from = target_gene, to = common_nodes)
g <- graph_from_data_frame(edges, directed = FALSE)

# گرافیک شبکه
V(g)$color <- ifelse(V(g)$name == target_gene, "#FF3333", "#66CCFF")
V(g)$size <- ifelse(V(g)$name == target_gene, 35, 12)
V(g)$label.cex <- 0.7
V(g)$label.font <- 2

plot(g, 
     layout = layout_with_kk(g),
     vertex.frame.color = "white",
     edge.color = "gray80",
     main = paste("Biological Network of", target_gene))


png("Venn_Biological.png", width=900, height=900, res=130)
grid.draw(venn_plot)
dev.off()

png("Network_Biological.png", width=900, height=900, res=130)
plot(g, layout = layout_with_kk(g))
dev.off()