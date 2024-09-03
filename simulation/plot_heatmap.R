mat <- readRDS("simulation/v3/multi_fitted.rds")
pseudotime <- readRDS("Covid_Su/data/tActivate_pseudotime.rds")
cellanno = sub(':.*', '', names(pseudotime))
names(cellanno) <- names(pseudotime)
samp_group <- readRDS("simulation/v3/metasamp.rds")
cell_group <- samp_group[cellanno]
names(cell_group) <- names(cellanno)
head(cell_group)

noise = 0.1
set.seed(2024)
mat <- mat + matrix(rnorm(length(mat), mean = 0, sd = noise), nrow = nrow(mat), ncol = ncol(mat))
pl <- (mat - Matrix::rowMeans(mat)) / matrixStats::rowSds(mat)

suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))

# Without group label----
Gene.set <- factor(rep(c("Set 1", "Set 2", "Set 3"), times = c(10, 10, 5)))
Gene.clv <- brewer.pal(8, "Set2")[1:3]
names(Gene.clv) <- levels(Gene.set)

png(file = "simulation/v3/noise_0.1/simulation_heatmap_time.png", res = 300, units = "in", width = 4, height = 3.5)
ht <- Heatmap(matrix = pl[, names(pseudotime)], name = "Expression",
              col = colorRamp2(c(quantile(pl, 0.05), 0, quantile(pl, 0.95)), brewer.pal(7, "RdYlBu")[c(7, 5, 1)]),
              cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F,
              top_annotation = columnAnnotation(Pseudotime = pseudotime,
                                                col = list(Pseudotime = colorRamp2(range(pseudotime), c("white", "hotpink")))),
              left_annotation = rowAnnotation(Gene = Gene.set,
                                              col = list(Gene = Gene.clv)),
              row_split = Gene.set, column_split = NULL,
              use_raster = F, heatmap_legend_param = list(legend_direction = "horizontal"))
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "right", legend_grouping = "original")
dev.off()

# With group label----
Gene.set <- factor(rep(c("Set 1", "Set 2", "Set 3"), times = c(10, 10, 5)))
Gene.clv <- brewer.pal(8, "Set2")[1:3]
names(Gene.clv) <- levels(Gene.set)
Sample.group <- factor(rep(c("Group 1", "Group 2", "Group 3"), times = as.numeric(table(cell_group))))
Sample.clv <- c(brewer.pal(9, "Set1")[3], brewer.pal(11, "Spectral")[c(11, 10)])
names(Sample.clv) <- levels(Sample.group)

png(file = "simulation/v3/noise_0.1/simulation_heatmap_time_group.png", res = 300, units = "in", width = 5, height = 4.5)
ht <- Heatmap(matrix = pl, name = "Expression",
              col = colorRamp2(c(quantile(pl, 0.05), 0, quantile(pl, 0.95)), brewer.pal(7, "RdYlBu")[c(7, 5, 1)]),
              cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F,
              top_annotation = columnAnnotation(Sample = Sample.group,
                                                Pseudotime = c(pseudotime[cell_group == "Group 1"],
                                                               pseudotime[cell_group == "Group 2"],
                                                               pseudotime[cell_group == "Group 3"]),
                                                col = list(Pseudotime = colorRamp2(range(pseudotime), c("white", "hotpink")),
                                                           Sample = Sample.clv)),
              left_annotation = rowAnnotation(Gene = Gene.set,
                                              col = list(Gene = Gene.clv)),
              row_split = Gene.set, column_split = Sample.group,
              use_raster = F, heatmap_legend_param = list(legend_direction = "horizontal"))
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "right", legend_grouping = "original")
dev.off()
