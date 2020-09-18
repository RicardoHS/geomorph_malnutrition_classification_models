library(geomorph) # procrustes but with geomorph https://rdrr.io/cran/geomorph/man/gpagen.html
library(tidyverse)
library(corrplot)

sed_image = function(filepath){
  # swap spaces for hyphens in lines beginning with "IMAGE="
  sed_file = NULL
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if(substr(line,1,6)=='IMAGE='){
      line=gsub(" ", "-", line)
    }
    sed_file = paste(sed_file, line, sep='\n')
  }
  close(con)
  sed_file
}


read_tps = function(file_path, rows_by_case = 24, n_cases = 10, remove_point_5=TRUE){
# This function reads the data from the .tps file and transforms it into an ordered dataframe
  # See the tps as if it were a .csv, the file is not quite correct but it is a good entry point.
  sed_file = sed_image(file_path)
  data = read.csv(text=sed_file, sep = " ", header = FALSE)
  
  # this part selects the line of the tps that starts with IMAGE = and performs a series of transformations, in order:
  # data [((1: n_cases) * rows_by_case) -1,] - select the rows. This returns a dataframe with two columns when it should be one (because the read.csv doesn't do all the work for us)
  # unite ('image', sep = "") - Join the two columns in a new image call.
  # t () - Transpose the dataframe for the next transformation
  # as.vector () - converts the transposed dataframe to a vector, which is what the following functions need
  # str_remove ('IMAGE =') - Remove IMAGE = from all vector elements
  # str_replace ("", "-") - Change whitespace to hyphens (some town names have spaces and that's not good for data transformations)
  # str_split ('_', simplify = TRUE) - Split each element into columns based on the underscore to split it. This generates an array
  # as.data.frame () - Converts array to dataframe
  # generate the tidy df from the IMAGE str
  df = data[((1:n_cases)*rows_by_case)-1,] %>% unite('image', sep="") %>% t() %>% as.vector() %>% 
    str_remove('IMAGE=') %>% 
    str_replace(" ", "-") %>% 
    str_split('_', simplify = TRUE) %>% 
    as.data.frame()
  # names the columns of the dataframe
  colnames(df) = c('country', 'region','village','body_part','participant', 'sex', 'age_group', 'crit', 'diagnosis', 'ethny')
  
  # add the scales - add a new column called scales to the dataframe
  df$scales = data[(1:n_cases)*rows_by_case, 1] %>% t() %>% as.vector() %>% str_remove('SCALE=') %>% as.numeric()
  
  # add the data points
  # another read.csv problem. All LMs have been treated as factors rather than numerical.
  # here we calculate the indices of the rows of the dataframe that are not LM, that is, the lines that start with LM= IMAGE= SCALE=
  excluded_rows = c(((0:(n_cases-1))*rows_by_case)+1, ((0:(n_cases-1))*rows_by_case)+23, ((0:(n_cases-1))*rows_by_case)+24)
  # here we select with -excluded_rows all the columns that are not the excluded cases of the previous line, that is, the MLs and convert them to numeric
  data_points = data[-excluded_rows,] %>% lapply(function(x) as.numeric(as.character(x)))
  
  # to add all the points to the dataframe more simply.
  # We generate a matrix with the x points and convert it to a dataframe
  x = matrix(data_points$V1, nrow = n_cases, byrow = TRUE) %>% as.data.frame()
  if(remove_point_5){
    x = x[,-5]
  }
  #we name the columns
  colnames(x) = sprintf("x_%s",seq(1:dim(x)[2]))
  # we generate a matrix with the points and convert it to dataframe
  y = matrix(-data_points$V2, nrow = n_cases, byrow = TRUE)%>% as.data.frame()
  if(remove_point_5){
    y = y[,-5]
  }
  colnames(y) = sprintf("y_%s",seq(1:dim(y)[2]))
  
  # we put together the dataframe we had with the dataframe of the LM points
  list(df = cbind(df,x,y), x_colnames=colnames(x), y_colnames=colnames(y))
}


perform_procrustes = function(df, x_coord_names, y_coord_names, plot_procrustes=FALSE, return_object=FALSE){
  # perform generalized procrustes analysis
  A = array(NA, dim=c(length(x_coord_names),2,nrow(df)))
  for(i in 1:nrow(df)){
    A[,1,i] = df[i,] %>% dplyr::select(x_coord_names) %>% as.vector() %>% as.numeric() 
    A[,2,i] = df[i,] %>% dplyr::select(y_coord_names) %>% as.vector() %>% as.numeric() 
  }
  
  # perform gpa
  gpa = gpagen(A)
  p_gpa = gpa$coords
  if(plot_procrustes){
    plot(gpa)
  }
  
  # we add the procurstes coordinates to the dataframe
  for(i in 1:length(x_coord_names)){
    df[, sprintf('x_%s_proc', i)] = as.vector(p_gpa[i,1,])
    df[, sprintf('y_%s_proc', i)] = as.vector(p_gpa[i,2,])
  }
  
  if(return_object){
    gpa
  }else{
    df 
  }
}


get_noncolinear_columns = function(df, columns, max_percentage=0.80, plot_correlation=FALSE){
  # Feautre Colinearity Trim: get rid of the coordinates that have a very high correlation
  if(!is.null(columns)){
    df = df %>% dplyr::select(columns)
  }else{
    df = df
  }
  
  tmp <- cor(df)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  
  data.new <- df[,!apply(tmp,2,function(x) any(abs(x) > max_percentage))]
  non_colinear_cols = colnames(data.new)
  
  if(plot_correlation){
    df %>% dplyr::select(non_colinear_cols) %>% cor() %>% corrplot::corrplot(method='number', type = 'upper')
  }
  
  non_colinear_cols
}

get_distances_colnames = function(n_cols){
  # Compute the distance matrix column names (just the names no the matrix)
  names = NULL
  for(i in 1:(n_cols-1)){
    for(j in (i+1):n_cols){
      names = c(names, paste('LM',i,'toLM',j, sep=''))
    }
  }
  return(names)
}

get_all_distances = function(df, x_coord_names, y_coord_names){
  # Compute the distance matrix
  distances = array(NA, dim=c(nrow(df), sum(1:(length(x_coord_names)-1))))
  for(i in 1:nrow(df)){
    distances[i,] = df[i,]%>% dplyr::select(x_coord_names, y_coord_names) %>% matrix(nrow = 2, byrow = TRUE) %>% t() %>% dist() %>% as.vector()
  }
  colnames(distances) = get_distances_colnames(length(x_coord_names))
  data.frame(distances)
}

get_train_test_split = function(df, split_percentage=1/3){
  # split data on train-test
  train_index <- sample(1:nrow(df), split_percentage * nrow(df))
  t_train <- df[train_index, ]
  t_test <- df[-train_index, ]
  list(train=t_train, test=t_test)
}

test_svmRadial = function(df, x_columns, data_filter, test_name, 
                          results_cols = c("Accuracy","Specificity", "Precision", "Recall", "F1"),
                          split_percetange=0.66,  seed=42){
  # perform SVM test using simple train-test split
  t_data = df %>% data_filter
  t_data$diagnosis = droplevels(t_data$diagnosis)
  
  set.seed(seed)
  split = get_train_test_split(t_data, split_percentage = split_percetange)
  
  (formula = as.formula(paste('diagnosis ~ ', paste(x_columns, collapse = '+'))))
  clf <- train(formula, data = split$train, method = "svmRadial")
  preds = predict(clf, split$test)
  cm = caret::confusionMatrix(data=preds, reference = split$test$diagnosis)
  result = c(cm$overall,cm$byClass) %>% t() %>% as.data.frame() %>% dplyr::select(results_cols)
  rownames(result) = test_name
  result
}

test_LOOCV_svmRadial = function(df, x_columns, data_filter, test_name, 
                          results_cols = c("Accuracy",  "Specificity", "Precision", "Recall", "F1"), seed=42, return_model=FALSE){
  # perform SVM test using leave one out cross validation
  t_data = df %>% data_filter
  t_data$diagnosis = droplevels(t_data$diagnosis)
  
  set.seed(seed)
  
  (formula = as.formula(paste('diagnosis ~ ', paste(x_columns, collapse = '+'))))
  train_control<- trainControl(method="LOOCV", savePredictions = 'final')
  clf <- train(formula, data = t_data, trControl=train_control, method = "svmRadial")
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

test_LOOCV_rf = function(df, x_columns, data_filter, test_name, 
                                results_cols = c("Accuracy",  "Specificity", "Precision", "Recall", "F1"), seed=42, return_model=FALSE){
  # perform random forest test using leave one out cross validation
  t_data = df %>% data_filter
  t_data$diagnosis = droplevels(t_data$diagnosis)
  
  set.seed(seed)
  
  (formula = as.formula(paste('diagnosis ~ ', paste(x_columns, collapse = '+'))))
  train_control<- trainControl(method="LOOCV", savePredictions = 'final')
  clf <- train(formula, data = t_data, trControl=train_control, method = "rf", ntree = 50)
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

test_LOOCV_nnet = function(df, x_columns, data_filter, test_name, 
                         results_cols = c("Accuracy",  "Specificity", "Precision", "Recall", "F1"), seed=42, return_model=FALSE){
  # perform neural network test using leave one out cross validation
  t_data = df %>% data_filter
  t_data$diagnosis = droplevels(t_data$diagnosis)
  
  set.seed(seed)
  
  (formula = as.formula(paste('diagnosis ~ ', paste(x_columns, collapse = '+'))))
  train_control<- trainControl(method = 'cv', number = 10, classProbs = TRUE, verboseIter = FALSE, savePredictions = 'final')
  clf <- train(formula, data = t_data, trControl=train_control, tuneGrid=expand.grid(size=c(2,4,5), decay=c(0.1,0.01)), method = "nnet", trace = FALSE) #,linout = TRUE)
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

perform_pca_and_filter = function(df, columns, max_explained_var=0.95){
  # compute PCA and reduce dimensionality till max_explained_var is reached
  if(!is.null(columns)){
    df = df %>% dplyr::select(columns)
  }else{
    df = df
  }
  pca <- prcomp(df, scale = TRUE)
  prop_var <- pca$sdev^2 / sum(pca$sdev^2)
  prop_var_cum <- cumsum(prop_var)
  num_pcs = which(prop_var_cum > max_explained_var)[1]
  pca$x[,1:num_pcs]
}

remove_allometric_shape_effect = function(df, procrustes_object){
  # remove the allometric shape effect
  po = procrustes_object
  x_df = data.frame(lm(matrix(po$coords[,1,],nrow = length(po$coords[1,1,])) ~ po$Csize)$residuals)
  y_df = data.frame(lm(matrix(po$coords[,2,],nrow = length(po$coords[1,1,])) ~ po$Csize)$residuals)
  colnames(x_df) = sprintf("x_%s_allom",1:dim(po$coords)[1])
  colnames(y_df) = sprintf("y_%s_allom",1:dim(po$coords)[1])
  cbind(df, x_df, y_df)
}
