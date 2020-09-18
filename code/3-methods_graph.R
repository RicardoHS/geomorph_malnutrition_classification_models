# script to compute the methodological graph combination


library(igraph)

methods_graph = graph_from_literal(DATA-+DISTANCES, 
                                   DATA-+PROCRUSTES, 
                                   PROCRUSTES-+DISTANCES,
                                   PROCRUSTES-+ALLOMETRIC_SHAPE_REDUCTION,
                                   PROCRUSTES-+CLASSIFICATION_MODEL,
                                   PROCRUSTES-+COLINEARITY_TRIM,
                                   PROCRUSTES-+PCA,
                                   DISTANCES-+COLINEARITY_TRIM,
                                   DISTANCES-+PCA,
                                   DISTANCES-+CLASSIFICATION_MODEL,
                                   PCA-+CLASSIFICATION_MODEL,
                                   ALLOMETRIC_SHAPE_REDUCTION-+COLINEARITY_TRIM,
                                   ALLOMETRIC_SHAPE_REDUCTION-+CLASSIFICATION_MODEL,
                                   COLINEARITY_TRIM-+CLASSIFICATION_MODEL)


plot(methods_graph)

all_paths = all_simple_paths(methods_graph, from=1, to=5)
