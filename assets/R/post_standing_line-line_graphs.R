library(igraph)

red_nodes <- c(1, 4, 6)                       # 1 = boundary, 4 and 6 = interior
cols      <- ifelse(1:8 %in% red_nodes, "red", "purple")

# left: one vertical line of 8 nodes
g1 <- make_graph(~ 1-2-3-4-5-6-7-8)
l1 <- cbind(rep(0, 8), 8:1)

# right: the same line cut in half -> two lines of 4 nodes
g2 <- make_graph(~ 1-2-3-4, 5-6-7-8)
l2 <- cbind(c(rep(-0.5, 4), rep(0.5, 4)), c(8:5, 8:5))

png("line_graphs_igraph.png", width = 1500, height = 1400, res = 200)
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))

plot(g1, layout = l1, rescale = FALSE, xlim = c(-1.5, 1.5), ylim = c(0.3, 8.7),
     asp = 0, vertex.color = cols, vertex.size = 18, vertex.label.color = "white",
     vertex.label.font = 2, vertex.frame.color = "black",
     edge.color = "grey30", edge.width = 3, main = "Initial queue")
segments(-0.7, 4.5, 0.7, 4.5, lty = 2); text(0.85, 4.5, "cut", cex = 0.8)

plot(g2, layout = l2, rescale = FALSE, xlim = c(-1.5, 1.5), ylim = c(0.3, 8.7),
     asp = 0, vertex.color = cols, vertex.size = 18, vertex.label.color = "white",
     vertex.label.font = 2, vertex.frame.color = "black",
     edge.color = "grey30", edge.width = 3, main = "Queue after splitting")

legend("bottom", legend = c("susceptible", "infectious"), pch = 21,
       pt.bg = c("purple", "red"), pt.cex = 2, horiz = TRUE, bty = "n")
dev.off()
