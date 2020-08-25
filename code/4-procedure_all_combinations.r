library(MASS) # para lda
library(caret) # models and classification utils
library(igraph)

source('functions.r')

file_path = './data/Muestra_563_UC3M_fixed.TPS'

tps = read_tps(file_path, remove_point_5 = TRUE, n_cases = 563)
x_coord_names = tps$x_colnames
y_coord_names = tps$y_colnames
x_procrustes_coord_names = sprintf("x_%s_proc",1:length(x_coord_names))
y_procrustes_coord_names = sprintf("y_%s_proc",1:length(y_coord_names))
x_allom_coord_names = sprintf("x_%s_allom",1:length(x_coord_names))
y_allom_coord_names = sprintf("y_%s_allom",1:length(y_coord_names))
x_extra_columns = c('sex', 'age_group')

diagnosis_combinations = combn(levels(tps$df$diagnosis), 2, simplify = FALSE)

##### DEFINE GRAPH FUNCTIONS

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

CLASSIFICATION_MODEL = function(data){
  x_columns = c(colnames(data), x_extra_columns)
  x_data = cbind(data, DATA() %>% select(x_extra_columns, 'diagnosis'))
  
  res_svm = test_LOOCV_svmRadial(x_data, x_columns = x_columns, data_filter = data_filter, test_name=paste(t_name,'->SVM',sep=''))
  res_rf = test_LOOCV_rf(x_data, x_columns = x_columns, data_filter = data_filter, test_name=paste(t_name,'->RF',sep=''))
  res_nnet = test_LOOCV_nnet(x_data, x_columns = x_columns, data_filter = data_filter, test_name=paste(t_name,'->NNET',sep=''))
  rbind(res_svm, res_rf, res_nnet)
}

CLASSIFICATION_MODEL_FULL = function(data){
  x_columns = c(colnames(data), x_extra_columns)
  x_data = cbind(data, DATA() %>% select(x_extra_columns, 'diagnosis'))
  
  res_svm = test_LOOCV_svmRadial(x_data, x_columns = x_columns, data_filter = data_filter, test_name=paste(t_name,'->SVM',sep=''), results_cols = c("Accuracy"))
  res_rf = test_LOOCV_rf(x_data, x_columns = x_columns, data_filter = data_filter, test_name=paste(t_name,'->RF',sep=''), results_cols = c("Accuracy"))
  res_nnet = test_LOOCV_nnet(x_data, x_columns = x_columns, data_filter = data_filter, test_name=paste(t_name,'->NNET',sep=''), results_cols = c("Accuracy"))
  rbind(res_svm, res_rf, res_nnet)
}

#####

# Binary-fashion way of modeling
for(comb in diagnosis_combinations){
  cat('Diagnosis:', comb, '\n')
  data_filter = function(x) x %>% filter(diagnosis==comb[1] | diagnosis==comb[2])
  
  results_cols = c("Accuracy", "Specificity", "Precision", "Recall", "F1")
  results = data.frame(matrix(ncol = length(results_cols), nrow = 0))
  colnames(results) = results_cols
  
  ##################################################################################
  ##################################################################################
  
  for(path in all_paths){
    t_name = path$name[-c(1,length(path$name))] %>% paste(collapse='->')
    cat(t_name,'\n')
    
    method_data = DATA()
    for(stage in path$name[-1]){
      method_data = match.fun(stage)(method_data)
    }
    results = results %>% rbind(method_data)
  }
  
  print(results)
  write.csv(results,file=paste("results_",paste(comb,collapse='_'),'.csv', sep=''))
}

#################################################################################
# Multiclass models
data_filter = function(x) x


##################################################################################
##################################################################################

for(path in all_paths){
  
  results_cols = c("Accuracy")
  results = data.frame(matrix(ncol = length(results_cols), nrow = 0))
  colnames(results) = results_cols
  
  t_name = path$name[-c(1,length(path$name))] %>% paste(collapse='->')
  cat(t_name,'\n')
  
  method_data = DATA()
  for(stage in path$name[-c(1,length(path$name))]){
    method_data = match.fun(stage)(method_data)
  }
  results = results %>% rbind(CLASSIFICATION_MODEL_FULL(method_data))
}

print(results)
write.csv(results,file="results_ALL_CLASES.csv")
