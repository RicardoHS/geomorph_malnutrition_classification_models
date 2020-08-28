############## MODELING
df = cbind(DISTANCES(DATA()), DATA() %>% dplyr::select(diagnosis, sex, age_group))
x_columns = c(colnames(DISTANCES(DATA())), 'sex', 'age_group')
df = df %>% dplyr::mutate(diagnosis=case_when(diagnosis!='SAM' ~ 'NoSAM', TRUE ~ 'SAM'))
df = df %>% data_filter
df$diagnosis = as.factor(df$diagnosis)
df$diagnosis = droplevels(df$diagnosis)



train_control<- trainControl(method="LOOCV", savePredictions = 'final')
clf <- train(formula, data = df, trControl=train_control, method = "rf", n_tree=50)
cm = caret::confusionMatrix(data=clf$pred$pred, reference = clf$pred$obs)
cm

############# PLOTS
plot(varImp(clf), top=30)
importance = varImp(clf)$importance[[1]]

coord_plots_x = DATA()[1,] %>% dplyr::select(x_coord_names) %>% as_vector() 
coord_plots_y = DATA()[1,] %>% dplyr::select(y_coord_names) %>% as_vector() 

plot(coord_plots_x, coord_plots_y, bg="red", pch=21, asp=1)
text(coord_plots_x, coord_plots_y, x_coord_names, cex=0.6, pos=4, col="red") 

most_important_distances = get_distances_colnames(length(x_coord_names))[order(importance, decreasing =TRUE)[1:4]]
for(mid in most_important_distances){
  points = mid %>% str_remove('LM') %>% str_split('toLM') %>% unlist() %>% as.numeric()
  segments(coord_plots_x[points[1]], coord_plots_y[points[1]], 
           coord_plots_x[points[2]], coord_plots_y[points[2]])
}
