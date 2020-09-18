# script to repeat the experiments but using the state of the art methods

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

diagnosis_combinations = combn(levels(tps$df$diagnosis), 2, simplify = FALSE)

DATA = function(data){
  tps$df
}

PROCRUSTES = function(data){
  procrustes_df = perform_procrustes(data, x_coord_names, y_coord_names, plot_procrustes=FALSE)
  procrustes_df %>% dplyr::select(x_procrustes_coord_names, y_procrustes_coord_names)
}

COLINEARITY_TRIM = function(data){
  x_columns = colnames(data)
  non_colinear_cols = get_noncolinear_columns(data, columns = x_columns, max_percentage=0.96, plot_correlation=FALSE)
  data %>% dplyr::select(non_colinear_cols)
}

ALLOMETRIC_SHAPE_REDUCTION = function(data){
  gpa = perform_procrustes(DATA(), x_coord_names, y_coord_names, plot_procrustes=FALSE, return_object=TRUE)
  data = remove_allometric_shape_effect(data, gpa)
  data %>% dplyr::select(c(x_allom_coord_names, y_allom_coord_names))
}

LDA = function(data, results_cols = c("Accuracy",  "Specificity", "Precision", "Recall", "F1"), seed=42){
  x_columns = c(colnames(data), x_extra_columns)
  t_data = cbind(data, DATA() %>% dplyr::select(x_extra_columns, 'diagnosis'))
  
  t_data = t_data %>% data_filter
  t_data$diagnosis = droplevels(t_data$diagnosis)
  
  set.seed(seed)
  
  (formula = as.formula(paste('diagnosis ~ ', paste(x_columns, collapse = '+'))))
  train_control<- trainControl(method = 'cv', number = 10, classProbs = TRUE, verboseIter = FALSE, savePredictions = 'final')
  clf <- train(formula, data = t_data, trControl=train_control, method = "lda", trace = FALSE) #,linout = TRUE)
  cm = caret::confusionMatrix(data=clf$pred$pred, reference = clf$pred$obs)
  result = c(cm$overall,cm$byClass) %>% t() %>% as.data.frame() %>% dplyr::select(results_cols)
  
  classes = unique(t_data$diagnosis)
  for(c in classes){
    total_c = clf$pred %>% filter(obs==c) %>% dim
    correct_c = clf$pred %>% filter(obs==c) %>% mutate(correct = if_else(pred == obs, 1, 0)) %>% dplyr::select(correct) %>% sum
    result = result %>% mutate(!!paste(" total_",c,sep=""):=total_c[1], !!paste("correct_",c,sep=""):=correct_c)
  }
  
  rownames(result) = test_name
  result
}

######################################################

for(x_extra_columns in list(NULL, c('sex'), c('age_group'), c('sex', 'age_group'))){
  for(comb in diagnosis_combinations){
    cat('Diagnosis:', comb, '\n')
    data_filter = function(x) x %>% filter(diagnosis==comb[1] | diagnosis==comb[2])
    
    results_cols = c("Accuracy", "Specificity", "Precision", "Recall", "F1")
    results = data.frame(matrix(ncol = length(results_cols), nrow = 0))
    colnames(results) = results_cols
    
    #PROCRUSTES + LDA
    test_name = "PROCRUSTES-LDA"
    method_data = LDA(PROCRUSTES(DATA()))
    results = results %>% rbind(method_data)
    #PROCRUSTES + ALLOM + LDA
    test_name = 'PROCRUSTES-ALLOM_REDUCTION-LDA'
    method_data = LDA(ALLOMETRIC_SHAPE_REDUCTION(PROCRUSTES(DATA())))
    results = results %>% rbind(method_data)
    
    print(results)
    write.csv(results,file=paste("../results/results_stateofart_",paste(x_extra_columns,collapse='-'),"_",paste(comb,collapse='_'),".csv", sep=''))
  }
}
