library('igraph')
library('umap')
library('tidyverse')

source('functions.r')

file_path = './data/Muestra_563_UC3M_fixed.TPS'
tps = read_tps(file_path, remove_point_5 = TRUE, n_cases = 563)
x_coord_names = tps$x_colnames
y_coord_names = tps$y_colnames
x_procrustes_coord_names = sprintf("x_%s_proc",1:length(x_coord_names))
y_procrustes_coord_names = sprintf("y_%s_proc",1:length(y_coord_names))
x_allom_coord_names = sprintf("x_%s_allom",1:length(x_coord_names))
y_allom_coord_names = sprintf("y_%s_allom",1:length(y_coord_names))

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

all_paths = all_simple_paths(methods_graph, from=1, to=5)


#################
DATA = function(data){
  tps$df
}

PROCRUSTES = function(data){
  procrustes_df = perform_procrustes(data, x_coord_names, y_coord_names, plot_procrustes=FALSE)
  procrustes_df %>% select(x_procrustes_coord_names, y_procrustes_coord_names)
}

DISTANCES = function(data){
  x_coords = colnames(data)[colnames(data) %>% str_detect('^x_')]
  y_coords = colnames(data)[colnames(data) %>% str_detect('^y_')]
  distances_raw = get_all_distances(data, x_coords, y_coords)
  distances_raw
}

PCA = function(data){
  pca_result = perform_pca_and_filter(data, columns=NULL, max_explained_var=0.95)
  pca_result
}

COLINEARITY_TRIM = function(data){
  x_columns = colnames(data)
  non_colinear_cols = get_noncolinear_columns(data, columns = x_columns, max_percentage=0.96, plot_correlation=FALSE)
  data %>% select(non_colinear_cols)
}

ALLOMETRIC_SHAPE_REDUCTION = function(data){
  gpa = perform_procrustes(DATA(), x_coord_names, y_coord_names, plot_procrustes=FALSE, return_object=TRUE)
  data = remove_allometric_shape_effect(data, gpa)
  data %>% select(c(x_allom_coord_names, y_allom_coord_names))
}
#################3

umap.config = umap::umap.defaults
umap.config$n_neighbors = 25

set.seed(42)
for(path in all_paths){
  t_name = path$name[-c(1,length(path$name))] %>% paste(collapse='->')
  method_data = DATA()
  for(stage in path$name[-c(1, length(path$name))]){
    method_data = match.fun(stage)(method_data)
  }
  
  df.umap = umap(scale(method_data), config = umap.config)
  gg.data = data.frame(x=df.umap$layout[,1], y=df.umap$layout[,2], diagnosis=tps$df$diagnosis)
  gg1 = ggplot(gg.data, aes(x=x,y=y, color=diagnosis)) + geom_jitter()  + geom_point(alpha = 1/10, size = 0.1, stroke = 0, shape = 16) + 
    ggtitle(t_name, subtitle = '') + theme(plot.title = element_text(size = 8))
  
  df.umap = cbind(method_data, DATA() %>% select(sex))
  df.umap['sex'] = lapply(df.umap['sex'], factor)
  df.umap['sex'] = lapply(df.umap['sex'], as.numeric)
  df.umap = umap(scale(df.umap), config = umap.config)
  gg.data = data.frame(x=df.umap$layout[,1], y=df.umap$layout[,2], diagnosis=tps$df$diagnosis)
  gg2 = ggplot(gg.data, aes(x=x,y=y, color=diagnosis)) + geom_jitter()  + geom_point(alpha = 1/10, size = 0.1, stroke = 0, shape = 16) + 
    ggtitle(t_name, subtitle = '+ sex') + theme(plot.title = element_text(size = 8))
  
  df.umap = cbind(method_data, DATA() %>% select(age_group))
  df.umap['age_group'] = lapply(df.umap['age_group'], factor)
  df.umap['age_group'] = lapply(df.umap['age_group'], as.numeric)
  df.umap = umap(scale(df.umap), config = umap.config)
  gg.data = data.frame(x=df.umap$layout[,1], y=df.umap$layout[,2], diagnosis=tps$df$diagnosis)
  gg3 = ggplot(gg.data, aes(x=x,y=y, color=diagnosis)) + geom_jitter()  + geom_point(alpha = 1/10, size = 0.1, stroke = 0, shape = 16) + 
    ggtitle(t_name, subtitle = '+ age_group') + theme(plot.title = element_text(size = 8))
  
  df.umap = cbind(method_data, DATA() %>% select(sex, age_group))
  df.umap['sex'] = lapply(df.umap['sex'], factor)
  df.umap['sex'] = lapply(df.umap['sex'], as.numeric)
  df.umap['age_group'] = lapply(df.umap['age_group'], factor)
  df.umap['age_group'] = lapply(df.umap['age_group'], as.numeric)
  df.umap = umap(scale(df.umap), config = umap.config)
  gg.data = data.frame(x=df.umap$layout[,1], y=df.umap$layout[,2], diagnosis=tps$df$diagnosis)
  gg4 = ggplot(gg.data, aes(x=x,y=y, color=diagnosis)) + geom_jitter()  + geom_point(alpha = 1/10, size = 0.1, stroke = 0, shape = 16) + 
    ggtitle(t_name, subtitle = '+ sex + age_group')+ theme(plot.title = element_text(size = 8))
  
  gg_c = ggpubr::ggarrange(gg1, gg2, gg3, gg4, nrow = 2, ncol=2, common.legend = TRUE, legend="right")
  file_name = paste('../cluster_images/', t_name %>% str_replace('->','_'), '.png', sep='')
  ggsave(plot=gg_c, filename=file_name, width=24, height = 24, units='cm', device = 'png')
}

