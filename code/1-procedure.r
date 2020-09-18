library(MASS) # for lda
library(caret) # models and classification utils

# custom utils
source('functions.r')

#here you have to change the path to use your own data. Probably you will need to tune the read_tps function inside functions.r
file_path = './data/Muestra_563_UC3M_fixed.TPS'

tps = read_tps(file_path, remove_point_5 = TRUE, n_cases = 563)

# some columns names to use it latter
x_coord_names = tps$x_colnames
y_coord_names = tps$y_colnames
x_procrustes_coord_names = sprintf("x_%s_proc",1:length(x_coord_names))
y_procrustes_coord_names = sprintf("y_%s_proc",1:length(y_coord_names))
x_allom_coord_names = sprintf("x_%s_allom",1:length(x_coord_names))
y_allom_coord_names = sprintf("y_%s_allom",1:length(y_coord_names))

# compute all the diagnosis types combinations
diagnosis_combinations = combn(levels(tps$df$diagnosis), 2, simplify = FALSE)

for(comb in diagnosis_combinations){
  print(comb)
  data_filter = function(x) x %>% filter(diagnosis==comb[1] | diagnosis==comb[2])
  
  results_cols = c("Accuracy", "Sensitivity","Specificity", "Precision", "Recall", "F1")
  results = data.frame(matrix(ncol = length(results_cols), nrow = 0))
  colnames(results) = results_cols
  
  ##################################################################################
  ##################################################################################
  
  #### Procrustes  -> SVM
  # test name
  t_name = 'procrustes_svm'
  df = tps$df
  
  df = perform_procrustes(df, x_coord_names, y_coord_names, plot_procrustes=FALSE)
  
  res = test_LOOCV_svmRadial(df, c(x_procrustes_coord_names, y_procrustes_coord_names), data_filter = data_filter, test_name=t_name)
  results = results %>% rbind(res)
  
  
  #### Procrustes -> allo -> SVM
  t_name = 'procrustes_allo_svm'
  df = tps$df
  
  gpa = perform_procrustes(df, x_coord_names, y_coord_names, plot_procrustes=FALSE, return_object=TRUE)
  df = remove_allometric_shape_effect(df, gpa)
  
  res = test_LOOCV_svmRadial(df, c(x_allom_coord_names, y_allom_coord_names), data_filter = data_filter, test_name=t_name)
  results = results %>% rbind(res)
  
  #### Procrustes -> Remove Colinear -> SVM
  t_name = 'procrustes_colinear_svm'
  df = tps$df
  
  df = perform_procrustes(df, x_coord_names, y_coord_names, plot_procrustes=FALSE)
  (non_colinear_cols = get_noncolinear_columns(df, columns = c(x_procrustes_coord_names, y_procrustes_coord_names), max_percentage=0.96, plot_correlation=FALSE))
  
  res = test_LOOCV_svmRadial(df, non_colinear_cols, data_filter = data_filter, test_name=t_name)
  results = results %>% rbind(res)
  
  #### Procrustes -> allo -> Remove Colinear -> SVM
  t_name = 'procrustes_allo_collinear_svm'
  df = tps$df
  
  gpa = perform_procrustes(df, x_coord_names, y_coord_names, plot_procrustes=FALSE, return_object=TRUE)
  df = remove_allometric_shape_effect(df, gpa)
  (non_colinear_cols = get_noncolinear_columns(df, columns = c(x_allom_coord_names, y_allom_coord_names), max_percentage=0.96, plot_correlation=FALSE))
  
  res = test_LOOCV_svmRadial(df, non_colinear_cols, data_filter = data_filter, test_name=t_name)
  results = results %>% rbind(res)
  
  #### Distances -> PCA95 -> SVM
  t_name = 'distances_pca95_svm'
  df = tps$df
  distances_raw = get_all_distances(df, x_coord_names, y_coord_names)
  distances_pca = perform_pca_and_filter(distances_raw, columns=NULL, max_explained_var=0.95)
  df = cbind(df, distances_pca)
  
  res = test_LOOCV_svmRadial(df, colnames(distances_pca), data_filter = data_filter, test_name=t_name)
  results = results %>% rbind(res)
  
  #### Distances -> PCA99 -> SVM
  t_name = 'distances_pca99_svm'
  df = tps$df
  distances_raw = get_all_distances(df, x_coord_names, y_coord_names)
  distances_pca = perform_pca_and_filter(distances_raw, columns=NULL, max_explained_var=0.99)
  df = cbind(df, distances_pca)
  
  res = test_LOOCV_svmRadial(df, colnames(distances_pca), data_filter = data_filter, test_name=t_name)
  results = results %>% rbind(res)
  
  #### Distances -> PCA999 -> SVM
  t_name = 'distances_pca999_svm'
  df = tps$df
  distances_raw = get_all_distances(df, x_coord_names, y_coord_names)
  distances_pca = perform_pca_and_filter(distances_raw, columns=NULL, max_explained_var=0.999)
  df = cbind(df, distances_pca)
  
  res = test_LOOCV_svmRadial(df, colnames(distances_pca), data_filter = data_filter, test_name=t_name)
  results = results %>% rbind(res)
  
  #### Distances -> Remove Colinear96 -> SVM
  t_name = 'distances_colinear96_svm'
  df = tps$df
  distances_raw = get_all_distances(df, x_coord_names, y_coord_names)
  (non_colinear_cols = get_noncolinear_columns(distances_raw, columns = NULL, max_percentage=0.96, plot_correlation=FALSE))
  df = cbind(df, distances_raw)
  
  res = test_LOOCV_svmRadial(df, non_colinear_cols, data_filter = data_filter, test_name=t_name)
  results = results %>% rbind(res)
  
  #### Distances -> Remove Colinear99 -> SVM
  t_name = 'distances_colinear99_svm'
  df = tps$df
  distances_raw = get_all_distances(df, x_coord_names, y_coord_names)
  (non_colinear_cols = get_noncolinear_columns(distances_raw, columns = NULL, max_percentage=0.99, plot_correlation=FALSE))
  df = cbind(df, distances_raw)
  
  res = test_LOOCV_svmRadial(df, non_colinear_cols, data_filter = data_filter, test_name=t_name)
  results = results %>% rbind(res)
  
  #### Distances -> SVM
  t_name = 'distances_svm'
  df = tps$df
  distances_raw = get_all_distances(df, x_coord_names, y_coord_names)
  df = cbind(df, distances_raw)
  
  res = test_LOOCV_svmRadial(df, colnames(distances_raw), data_filter = data_filter, test_name=t_name)
  results = results %>% rbind(res)
  
  ######## Procrustes -> Distances -> Remove Colinear -> SVM
  t_name = 'procrustes_distances_colinear_svm'
  df = tps$df
  
  df = perform_procrustes(df, x_coord_names, y_coord_names, plot_procrustes=FALSE)
  distances_procrustes = get_all_distances(df, x_procrustes_coord_names, y_procrustes_coord_names)
  (non_colinear_cols = get_noncolinear_columns(distances_procrustes, columns = NULL, max_percentage=0.96, plot_correlation=FALSE))
  df = cbind(df, distances_procrustes)
  
  res = test_LOOCV_svmRadial(df, non_colinear_cols, data_filter = data_filter, test_name=t_name)
  results = results %>% rbind(res)
  
  ######## Procrustes -> allo -> Distances -> Remove Colinear -> SVM
  t_name = 'procrustes_allo_distances_colinear_svm'
  df = tps$df
  
  gpa = perform_procrustes(df, x_coord_names, y_coord_names, plot_procrustes=FALSE, return_object=TRUE)
  df = remove_allometric_shape_effect(df, gpa)
  distances_procrustes = get_all_distances(df, x_allom_coord_names, y_allom_coord_names)
  (non_colinear_cols = get_noncolinear_columns(distances_procrustes, columns = NULL, max_percentage=0.96, plot_correlation=FALSE))
  df = cbind(df, distances_procrustes)
  
  res = test_LOOCV_svmRadial(df, non_colinear_cols, data_filter = data_filter, test_name=t_name)
  results = results %>% rbind(res)
  
  ##################################################################################
  ##################################################################################
  #diagnosis_table = df %>% data_filter %>% dplyr::select(diagnosis) %>% table
  #results = results %>% mutate(n_MAM=diagnosis_table['MAM'], n_NOR=diagnosis_table['NOR'], n_RIS=diagnosis_table['RIS'], n_SAM=diagnosis_table['SAM'])
  print(results)
}

