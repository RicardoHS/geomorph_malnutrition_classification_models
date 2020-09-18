# script used to generate the analysis by sex and age_group using the best models

library(MASS) # para lda
library(caret) # models and classification utils
library(igraph)
library(tidyverse)

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
  data %>% dplyr::select(non_colinear_cols)
}

ALLOMETRIC_SHAPE_REDUCTION = function(data){
  gpa = perform_procrustes(DATA(), x_coord_names, y_coord_names, plot_procrustes=FALSE, return_object=TRUE)
  data = remove_allometric_shape_effect(data, gpa)
  data %>% dplyr::select(c(x_allom_coord_names, y_allom_coord_names))
}

LDA = function(data, results_cols = c("Accuracy",  "Specificity", "Precision", "Recall", "F1"), seed=42, return_model=FALSE){
  x_columns = c(colnames(data), x_extra_columns)
  t_data = cbind(data, DATA() %>% dplyr::select(x_extra_columns, 'diagnosis'))
  
  t_data = t_data %>% data_filter
  t_data$diagnosis = droplevels(t_data$diagnosis)
  
  set.seed(seed)
  
  (formula = as.formula(paste('diagnosis ~ ', paste(x_columns, collapse = '+'))))
  train_control<- trainControl(method = 'cv', number = 10, classProbs = TRUE, verboseIter = FALSE, savePredictions = 'final')
  clf <- train(formula, data = t_data, trControl=train_control, method = "lda", trace = FALSE) #,linout = TRUE)
  if(return_model){
    return(clf)
  }
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

NNET = function(data, return_model=FALSE){
  x_columns = c(colnames(data), x_extra_columns)
  x_data = cbind(data, DATA() %>% dplyr::select(x_extra_columns, 'diagnosis'))
  test_LOOCV_nnet(x_data, x_columns = x_columns, data_filter = data_filter, test_name=paste(t_name,'->NNET',sep=''), return_model=return_model)
}

SVM = function(data, return_model=FALSE){
  x_columns = c(colnames(data), x_extra_columns)
  x_data = cbind(data, DATA() %>% dplyr::select(x_extra_columns, 'diagnosis'))
  test_LOOCV_svmRadial(x_data, x_columns = x_columns, data_filter = data_filter, test_name=paste(t_name,'->SVM',sep=''), return_model=return_model)
}

RF = function(data, return_model=FALSE){
  x_columns = c(colnames(data), x_extra_columns)
  x_data = cbind(data, DATA() %>% dplyr::select(x_extra_columns, 'diagnosis'))
  test_LOOCV_rf(x_data, x_columns = x_columns, data_filter = data_filter, test_name=paste(t_name,'->RF',sep=''), return_model=return_model)
}


###################################################### STATE-OF-THE-ART
plots = list()

## MAM-NOR
diag=c('MAM','NOR')
x_extra_columns = c('sex', 'age_group')
data_filter = function(x) x %>% filter(diagnosis==diag[1] | diagnosis==diag[2])
clf = LDA(PROCRUSTES(DATA()), return_model = TRUE)
is_correct = clf$pred$pred==clf$pred$obs
df = DATA() %>% data_filter
df[clf$pred$rowIndex,'is_correct'] = is_correct
plots[[1]] = ggplot(df, aes(x=sex, fill=factor(is_correct))) + 
                  geom_bar(position = 'fill') + 
                  ggtitle(paste(paste(diag, collapse='-'), ': State of the art results by Age Group')) + 
                  facet_grid(.~age_group) +
                  labs(fill='is_correct') + ylab('ratio')

## MAM-RIS
diag=c('MAM','RIS')
x_extra_columns = c('age_group')
data_filter = function(x) x %>% filter(diagnosis==diag[1] | diagnosis==diag[2])
clf = LDA(PROCRUSTES(DATA()), return_model = TRUE)
is_correct = clf$pred$pred==clf$pred$obs
df = DATA() %>% data_filter
df[clf$pred$rowIndex,'is_correct'] = is_correct
plots[[2]] = ggplot(df, aes(x=sex, fill=factor(is_correct))) + 
            geom_bar(position = 'fill') + 
            ggtitle(paste(paste(diag, collapse='-'), ': State of the art results by Age Group')) + 
            facet_grid(.~age_group) +
            labs(fill='is_correct') + ylab('ratio')

## MAM-SAM
diag=c('MAM','SAM')
x_extra_columns = c('sex', 'age_group')
data_filter = function(x) x %>% filter(diagnosis==diag[1] | diagnosis==diag[2])
clf = LDA(PROCRUSTES(DATA()), return_model = TRUE)
is_correct = clf$pred$pred==clf$pred$obs
df = DATA() %>% data_filter
df[clf$pred$rowIndex,'is_correct'] = is_correct
plots[[3]] = ggplot(df, aes(x=sex, fill=factor(is_correct))) + 
            geom_bar(position = 'fill') + 
            ggtitle(paste(paste(diag, collapse='-'), ': State of the art results by Age Group')) + 
            facet_grid(.~age_group) +
            labs(fill='is_correct') + ylab('ratio')

## NOR-RIS
diag=c('NOR','RIS')
x_extra_columns = c()
data_filter = function(x) x %>% filter(diagnosis==diag[1] | diagnosis==diag[2])
clf = LDA(PROCRUSTES(DATA()), return_model = TRUE)
is_correct = clf$pred$pred==clf$pred$obs
df = DATA() %>% data_filter
df[clf$pred$rowIndex,'is_correct'] = is_correct
plots[[4]] = ggplot(df, aes(x=sex, fill=factor(is_correct))) + 
            geom_bar(position = 'fill') + 
            ggtitle(paste(paste(diag, collapse='-'), ': State of the art results by Age Group')) + 
            facet_grid(.~age_group) +
            labs(fill='is_correct') + ylab('ratio')

## NOR-SAM
diag=c('NOR', 'SAM')
x_extra_columns = c('sex', 'age_group')
data_filter = function(x) x %>% filter(diagnosis==diag[1] | diagnosis==diag[2])
clf = LDA(PROCRUSTES(DATA()), return_model = TRUE)
is_correct = clf$pred$pred==clf$pred$obs
df = DATA() %>% data_filter
df[clf$pred$rowIndex,'is_correct'] = is_correct
plots[[5]] = ggplot(df, aes(x=sex, fill=factor(is_correct))) + 
            geom_bar(position = 'fill') + 
            ggtitle(paste(paste(diag, collapse='-'), ': State of the art results by Age Group')) + 
            facet_grid(.~age_group) +
            labs(fill='is_correct') + ylab('ratio')

## RIS-SAM
diag=c('RIS','SAM')
x_extra_columns = c('sex', 'age_group')
data_filter = function(x) x %>% filter(diagnosis==diag[1] | diagnosis==diag[2])
clf = LDA(PROCRUSTES(DATA()), return_model = TRUE)
is_correct = clf$pred$pred==clf$pred$obs
df = DATA() %>% data_filter
df[clf$pred$rowIndex,'is_correct'] = is_correct
plots[[6]] = ggplot(df, aes(x=sex, fill=factor(is_correct))) + 
            geom_bar(position = 'fill') + 
            ggtitle(paste(paste(diag, collapse='-'), ': State of the art results by Age Group')) + 
            facet_grid(.~age_group) +
            labs(fill='is_correct') + ylab('ratio')

gg_c = ggpubr::ggarrange(plotlist = plots, nrow = 2, ncol=3, common.legend = TRUE, legend="right")
print(gg_c)







###################################################### THESIS RESULTS
plots = list()

## MAM-NOR
diag=c('MAM','NOR')
x_extra_columns = c('age_group')
data_filter = function(x) x %>% filter(diagnosis==diag[1] | diagnosis==diag[2])
clf = NNET(PROCRUSTES(DATA()), return_model = TRUE)
is_correct = clf$pred$pred==clf$pred$obs
df = DATA() %>% data_filter
df[clf$pred$rowIndex,'is_correct'] = is_correct
plots[[1]] = ggplot(df, aes(x=sex, fill=factor(is_correct))) + 
  geom_bar(position = 'fill') + 
  ggtitle(paste(paste(diag, collapse='-'), ': Thesis results by Age Group')) + 
  facet_grid(.~age_group) +
  labs(fill='is_correct') + ylab('ratio')

## MAM-RIS
diag=c('MAM','RIS')
x_extra_columns = c('age_group')
data_filter = function(x) x %>% filter(diagnosis==diag[1] | diagnosis==diag[2])
clf = NNET(COLINEARITY_TRIM(PROCRUSTES(DATA())), return_model = TRUE)
is_correct = clf$pred$pred==clf$pred$obs
df = DATA() %>% data_filter
df[clf$pred$rowIndex,'is_correct'] = is_correct
plots[[2]] = ggplot(df, aes(x=sex, fill=factor(is_correct))) + 
  geom_bar(position = 'fill') + 
  ggtitle(paste(paste(diag, collapse='-'), ': Thesis results by Age Group')) + 
  facet_grid(.~age_group) +
  labs(fill='is_correct') + ylab('ratio')

## MAM-SAM
diag=c('MAM','SAM')
x_extra_columns = c('age_group')
data_filter = function(x) x %>% filter(diagnosis==diag[1] | diagnosis==diag[2])
clf = NNET(PCA(DISTANCES(DATA())), return_model = TRUE)
is_correct = clf$pred$pred==clf$pred$obs
df = DATA() %>% data_filter
df[clf$pred$rowIndex,'is_correct'] = is_correct
plots[[3]] = ggplot(df, aes(x=sex, fill=factor(is_correct))) + 
  geom_bar(position = 'fill') + 
  ggtitle(paste(paste(diag, collapse='-'), ': Thesis results by Age Group')) + 
  facet_grid(.~age_group) +
  labs(fill='is_correct') + ylab('ratio')

## NOR-RIS
diag=c('NOR','RIS')
x_extra_columns = c()
data_filter = function(x) x %>% filter(diagnosis==diag[1] | diagnosis==diag[2])
clf = SVM(COLINEARITY_TRIM(DISTANCES(DATA())), return_model = TRUE)
is_correct = clf$pred$pred==clf$pred$obs
df = DATA() %>% data_filter
df[clf$pred$rowIndex,'is_correct'] = is_correct
plots[[4]] = ggplot(df, aes(x=sex, fill=factor(is_correct))) + 
  geom_bar(position = 'fill') + 
  ggtitle(paste(paste(diag, collapse='-'), ': Thesis results by Age Group')) + 
  facet_grid(.~age_group) +
  labs(fill='is_correct') + ylab('ratio')

## NOR-SAM
diag=c('NOR', 'SAM')
x_extra_columns = c('sex', 'age_group')
data_filter = function(x) x %>% filter(diagnosis==diag[1] | diagnosis==diag[2])
clf = NNET(PCA(DISTANCES(DATA())), return_model = TRUE)
is_correct = clf$pred$pred==clf$pred$obs
df = DATA() %>% data_filter
df[clf$pred$rowIndex,'is_correct'] = is_correct
plots[[5]] = ggplot(df, aes(x=sex, fill=factor(is_correct))) + 
  geom_bar(position = 'fill') + 
  ggtitle(paste(paste(diag, collapse='-'), ': Thesis results by Age Group')) + 
  facet_grid(.~age_group) +
  labs(fill='is_correct') + ylab('ratio')

## RIS-SAM
diag=c('RIS','SAM')
x_extra_columns = c('sex', 'age_group')
data_filter = function(x) x %>% filter(diagnosis==diag[1] | diagnosis==diag[2])
clf = NNET(PCA(DISTANCES(DATA())), return_model = TRUE)
df = DATA() %>% data_filter
is_correct = clf$pred$pred==clf$pred$obs
df[clf$pred$rowIndex,'is_correct'] = is_correct
plots[[6]] = ggplot(df, aes(x=sex, fill=factor(is_correct))) + 
  geom_bar(position = 'fill') + 
  ggtitle(paste(paste(diag, collapse='-'), ': Thesis results by Age Group')) + 
  facet_grid(.~age_group) +
  labs(fill='is_correct') + ylab('ratio')

gg_c = ggpubr::ggarrange(plotlist = plots, nrow = 2, ncol=3, common.legend = TRUE, legend="right")
print(gg_c)

