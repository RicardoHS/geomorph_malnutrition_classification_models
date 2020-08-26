library('CCA')

library(MASS) # para lda
library(caret) # models and classification utils

source('functions.r')

file_path = './data/Muestra_563_UC3M_fixed.TPS'

tps = read_tps(file_path, remove_point_5 = TRUE, n_cases = 563)
x_coord_names = tps$x_colnames
y_coord_names = tps$y_colnames
x_procrustes_coord_names = sprintf("x_%s_proc",1:length(x_coord_names))
y_procrustes_coord_names = sprintf("y_%s_proc",1:length(y_coord_names))
x_allom_coord_names = sprintf("x_%s_allom",1:length(x_coord_names))
y_allom_coord_names = sprintf("y_%s_allom",1:length(y_coord_names))

data_filter = function(x) x %>% filter(diagnosis=='MAM' | diagnosis=='SAM')

results_cols = c("Accuracy", "Sensitivity","Specificity", "Precision", "Recall", "F1")
results = data.frame(matrix(ncol = length(results_cols), nrow = 0))
colnames(results) = results_cols
  
df = tps$df
df = perform_procrustes(df, x_coord_names, y_coord_names, plot_procrustes=FALSE)
  
# Canonical correlation of the landmarks

X = df %>% select(x_procrustes_coord_names)
Y = df %>% select(y_procrustes_coord_names)

canonical_correlation = cc(X, Y)

plt.cc(canonical_correlation, var.label = TRUE, ind.names = c(x_coord_names, y_coord_names))
