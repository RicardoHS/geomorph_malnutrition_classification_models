library(geomorph) # procrustes but with geomorph https://rdrr.io/cran/geomorph/man/gpagen.html
library(tidyverse) # para manejar las estructuras de datos mejor
library(corrplot)

sed_image = function(filepath){
  # cambia espacios por guiones en lineas que empiezan por "IMAGE="
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
# Esta funcion lee los datos del archivo .tps y los transforma en un dataframe ordenado 
  # lee el tps como si fuese un .csv, el archivo no es del todo correcto pero es un buen punto de entrada.
  sed_file = sed_image(file_path)
  data = read.csv(text=sed_file, sep = " ", header = FALSE)
  
  # esta parte selecciona la linea del tps que empieza por IMAGE= y le hace una serie de transformaciones, por orden:
  # data[((1:n_cases)*rows_by_case)-1,] - selecciona las rows. Esto devuelve un dataframe con dos columnas cuando deberia ser una (porq el read.csv no nos hace todo el trabajo)
  # unite('image', sep="") - Junta las dos columnas en una nueva llamada image.
  # t() - Transpone el dataframe para la siguiente transformacion
  # as.vector() - convierte el dataframe transpuesto en un vector, que es lo que necesitan las siguientes funciones
  # str_remove('IMAGE=') - Elimina IMAGE= de todos los elementos del vector
  # str_replace(" ", "-") - Cambia los espacios en blanco por guiones (algunos nombres de pueblo tienen espacio y eso no es bueno para las transformaciones de datos)
  # str_split('_', simplify = TRUE) - Divide cada elemento en columnas basandose en el guion bajo para dividirlo. Esto genera una matriz
  # as.data.frame() - Convierte la matriz en dataframe
  # generate the tidy df from the IMAGE str
  df = data[((1:n_cases)*rows_by_case)-1,] %>% unite('image', sep="") %>% t() %>% as.vector() %>% 
    str_remove('IMAGE=') %>% 
    str_replace(" ", "-") %>% 
    str_split('_', simplify = TRUE) %>% 
    as.data.frame()
  # pone nombres a las columnas del dataframe
  colnames(df) = c('country', 'region','village','body_part','participant', 'sex', 'age_group', 'crit', 'diagnosis', 'ethny')
  
  # add the scales - Añade una columna nueva llamada scales al dataframe 
  df$scales = data[(1:n_cases)*rows_by_case, 1] %>% t() %>% as.vector() %>% str_remove('SCALE=') %>% as.numeric()
  
  # add the data points
  # otro problema del read.csv. Todos los LM los ha tratado como factores en vez de numericos.
  # aqui calculamos los indices de las rows del dataframe que no son LM, es decir las lineas que empiezan con LM= IMAGE= SCALE=
  excluded_rows = c(((0:(n_cases-1))*rows_by_case)+1, ((0:(n_cases-1))*rows_by_case)+23, ((0:(n_cases-1))*rows_by_case)+24)
  # aqui seleccionamos con -excluded_rows todas las columnas que no son los casos excluidos de la linea anterior, es decir, los LM y las convertimos a numerico
  data_points = data[-excluded_rows,] %>% lapply(function(x) as.numeric(as.character(x)))
  
  # para añadir todos los puntos al dataframe de forma mas simple. 
  # Generamos una matriz con los puntos x y la convertimos a dataframe
  x = matrix(data_points$V1, nrow = n_cases, byrow = TRUE) %>% as.data.frame()
  if(remove_point_5){
    x = x[,-5]
  }
  # ponemos nombres a las columnas
  colnames(x) = sprintf("x_%s",seq(1:dim(x)[2]))
  # Generamos una matriz con los puntos y y la convertimos a dataframe
  y = matrix(-data_points$V2, nrow = n_cases, byrow = TRUE)%>% as.data.frame()
  if(remove_point_5){
    y = y[,-5]
  }
  colnames(y) = sprintf("y_%s",seq(1:dim(y)[2]))
  
  # juntamos el dataframe que teniamos con los dataframe de los puntos LM
  list(df = cbind(df,x,y), x_colnames=colnames(x), y_colnames=colnames(y))
}


perform_procrustes = function(df, x_coord_names, y_coord_names, plot_procrustes=FALSE, return_object=FALSE){
  A = array(NA, dim=c(length(x_coord_names),2,nrow(df)))
  for(i in 1:nrow(df)){
    A[,1,i] = df[i,] %>% dplyr::select(x_coord_names) %>% as.vector() %>% as.numeric() 
    A[,2,i] = df[i,] %>% dplyr::select(y_coord_names) %>% as.vector() %>% as.numeric() 
  }
  
  # ahcemos gpa
  gpa = gpagen(A)
  p_gpa = gpa$coords
  if(plot_procrustes){
    plot(gpa)
  }
  
  # añadimos las coordenadas procurstes al dataframe
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
  # nos libramos de las coordenadas que tienen una correlacion muy alta
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

get_all_distances = function(df, x_coord_names, y_coord_names){
  distances = array(NA, dim=c(nrow(df), sum(1:(length(x_coord_names)-1))))
  for(i in 1:nrow(df)){
    distances[i,] = df[i,]%>% dplyr::select(x_coord_names, y_coord_names) %>% matrix(nrow = 2, byrow = TRUE) %>% t() %>% dist() %>% as.vector()
  }
  
  data.frame(distances)
}

get_train_test_split = function(df, split_percentage=1/3){
  # seleccionamos los datos que queremos
  train_index <- sample(1:nrow(df), split_percentage * nrow(df))
  t_train <- df[train_index, ]
  t_test <- df[-train_index, ]
  list(train=t_train, test=t_test)
}

test_svmRadial = function(df, x_columns, data_filter, test_name, 
                          results_cols = c("Accuracy","Specificity", "Precision", "Recall", "F1"),
                          split_percetange=0.66,  seed=42){
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
                          results_cols = c("Accuracy",  "Specificity", "Precision", "Recall", "F1"), seed=42){
  t_data = df %>% data_filter
  t_data$diagnosis = droplevels(t_data$diagnosis)
  
  set.seed(seed)
  
  (formula = as.formula(paste('diagnosis ~ ', paste(x_columns, collapse = '+'))))
  train_control<- trainControl(method="LOOCV", savePredictions = 'final')
  clf <- train(formula, data = t_data, trControl=train_control, method = "svmRadial")
  cm = caret::confusionMatrix(data=clf$pred$pred, reference = clf$pred$obs)
  result = c(cm$overall,cm$byClass) %>% t() %>% as.data.frame() %>% dplyr::select(results_cols)
  
  classes = unique(t_data$diagnosis)
  for(c in classes){
    total_c = clf$pred %>% filter(obs==c) %>% dim
    correct_c = clf$pred %>% filter(obs==c) %>% mutate(correct = if_else(pred == obs, 1, 0)) %>% select(correct) %>% sum
    result = result %>% mutate(!!paste(" total_",c,sep=""):=total_c[1], !!paste("correct_",c,sep=""):=correct_c)
  }
  
  rownames(result) = test_name
  result
}

test_LOOCV_rf = function(df, x_columns, data_filter, test_name, 
                                results_cols = c("Accuracy",  "Specificity", "Precision", "Recall", "F1"), seed=42){
  t_data = df %>% data_filter
  t_data$diagnosis = droplevels(t_data$diagnosis)
  
  set.seed(seed)
  
  (formula = as.formula(paste('diagnosis ~ ', paste(x_columns, collapse = '+'))))
  train_control<- trainControl(method="LOOCV", savePredictions = 'final')
  clf <- train(formula, data = t_data, trControl=train_control, method = "rf", ntree = 50)
  cm = caret::confusionMatrix(data=clf$pred$pred, reference = clf$pred$obs)
  result = c(cm$overall,cm$byClass) %>% t() %>% as.data.frame() %>% dplyr::select(results_cols)
  
  classes = unique(t_data$diagnosis)
  for(c in classes){
    total_c = clf$pred %>% filter(obs==c) %>% dim
    correct_c = clf$pred %>% filter(obs==c) %>% mutate(correct = if_else(pred == obs, 1, 0)) %>% select(correct) %>% sum
    result = result %>% mutate(!!paste(" total_",c,sep=""):=total_c[1], !!paste("correct_",c,sep=""):=correct_c)
  }
  
  rownames(result) = test_name
  result
}

test_LOOCV_nnet = function(df, x_columns, data_filter, test_name, 
                         results_cols = c("Accuracy",  "Specificity", "Precision", "Recall", "F1"), seed=42){
  t_data = df %>% data_filter
  t_data$diagnosis = droplevels(t_data$diagnosis)
  
  set.seed(seed)
  
  (formula = as.formula(paste('diagnosis ~ ', paste(x_columns, collapse = '+'))))
  train_control<- trainControl(method = 'cv', number = 10, classProbs = TRUE, verboseIter = FALSE, savePredictions = 'final')
  clf <- train(formula, data = t_data, trControl=train_control, tuneGrid=expand.grid(size=c(2,4,5), decay=c(0.1,0.01)), method = "nnet", trace = FALSE) #,linout = TRUE)
  cm = caret::confusionMatrix(data=clf$pred$pred, reference = clf$pred$obs)
  result = c(cm$overall,cm$byClass) %>% t() %>% as.data.frame() %>% dplyr::select(results_cols)
  
  classes = unique(t_data$diagnosis)
  for(c in classes){
    total_c = clf$pred %>% filter(obs==c) %>% dim
    correct_c = clf$pred %>% filter(obs==c) %>% mutate(correct = if_else(pred == obs, 1, 0)) %>% select(correct) %>% sum
    result = result %>% mutate(!!paste(" total_",c,sep=""):=total_c[1], !!paste("correct_",c,sep=""):=correct_c)
  }
  
  rownames(result) = test_name
  result
}

perform_pca_and_filter = function(df, columns, max_explained_var=0.95){
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
  po = procrustes_object
  x_df = data.frame(lm(matrix(po$coords[,1,],nrow = length(po$coords[1,1,])) ~ po$Csize)$residuals)
  y_df = data.frame(lm(matrix(po$coords[,2,],nrow = length(po$coords[1,1,])) ~ po$Csize)$residuals)
  colnames(x_df) = sprintf("x_%s_allom",1:dim(po$coords)[1])
  colnames(y_df) = sprintf("y_%s_allom",1:dim(po$coords)[1])
  cbind(df, x_df, y_df)
}
