install.packages("randomForest")
install.packages("ggraph")
install.packages("igraph")
library(randomForest)
# Create the Random Forest model.
# The randomForest function accepts a "formula" structure as its main parameter. In this case, "Species" will be the variable
# to be predicted, while the others will be the predictors.
myLittleForest <- randomForest(Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, data = iris)

# Print the summary of our model.
print(myLittleForest)
print(importance(myLittleForest, type=2))
plot(myLittleForest)

my_data1 <- read.csv("https://ibm.box.com/shared/static/fzceg5vdj9hxpf7aopgvfgobi1g4vb4v.csv")

head(my_data1)
#Sldprice - House sale price
#rooms - Number of rooms
#beds - No of bedrooms
#d_cbd - Distance to centre of town
#hway_1 - Within 5 km of highway
#sway_1 - Within 1 km of subway
#hh_avinc - average household income
#detach - detached
#brick - brick
#air_con - air condition
#bsmt_fin - finished basement

plot(my_data1$sldprice)
## removing NAs from the data 

new_data <- na.omit(my_data1)
fit1 <- randomForest(sldprice~hh_avinc+rooms+beds+sway_1+hway_1+d_cbd+detach+air_con+brick+bsmt_fin,data=new_data,importance=TRUE)

print(fit1)

round(importance(fit1,type=1),2)
fit2 <- randomForest(sldprice~hh_avinc+rooms+beds+sway_1+hway_1+d_cbd+detach+air_con+brick+bsmt_fin,data=new_data,proximity=TRUE,action=na.omit,importance=TRUE)
print(fit2)
round(importance(fit2,type=1),2)

par(mfrow=c(1,2))
    plot(fit1)
plot(fit2)

library(dplyr)
library(ggraph)
library(igraph)

tree_func <- function(final_model, 
                      tree_num) {
  
  # get tree by index
  tree <- randomForest::getTree(final_model, 
                                k = tree_num, 
                                labelVar = TRUE) %>%
    tibble::rownames_to_column() %>%
    # make leaf split points to NA, so the 0s won't get plotted
    mutate(`split point` = ifelse(is.na(prediction), `split point`, NA))
  
  # prepare data frame for graph
  graph_frame <- data.frame(from = rep(tree$rowname, 2),
                            to = c(tree$`left daughter`, tree$`right daughter`))
  
  # convert to graph and delete the last node that we don't want to plot
  graph <- graph_from_data_frame(graph_frame) %>%
    delete_vertices("0")
  
  # set node labels
  V(graph)$node_label <- gsub("_", " ", as.character(tree$`split var`))
  V(graph)$leaf_label <- as.character(tree$prediction)
  V(graph)$split <- as.character(round(tree$`split point`, digits = 2))
  
  # plot
  plot <- ggraph(graph, 'dendrogram') + 
    theme_bw() +
    geom_edge_link() +
    geom_node_point() +
    geom_node_text(aes(label = node_label), na.rm = TRUE, repel = TRUE) +
    geom_node_label(aes(label = split), vjust = 2.5, na.rm = TRUE, fill = "white") +
    geom_node_label(aes(label = leaf_label, fill = leaf_label), na.rm = TRUE, 
                    repel = TRUE, colour = "white", fontface = "bold", show.legend = FALSE) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 18))
  
  print(plot)
}


tree_func(fit2,10)






