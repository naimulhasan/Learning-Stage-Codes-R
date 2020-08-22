install.packages("class")
install.packages("caret")
install.packages("mlbench")
install.packages("e1071")

library(class)
library(caret)
library(mlbench)
require(mlbench)
library(e1071)
library(base)
require(base)
data(Sonar)
head(Sonar)
View(Sonar)

cat("number of rows and columns are:", nrow(Sonar), ncol(Sonar))

?base::table

base::table(Sonar$Class) 
apply(Sonar, 2, function(x) sum(is.na(x))) 
#for a matrix 1 indicates rows, 2 indicates columns, c(1, 2) indicates rows and columns
?apply

SEED <- 123
set.seed(SEED)
data <- Sonar[base::sample(nrow(Sonar)), ] # shuffle data first
View(data)
bound <- floor(0.7 * nrow(data))
df_train <- data[1:bound, ] #70% at training data set
df_test <- data[(bound + 1):nrow(data), ] #30% at test data set
cat("number of training and test samples are ", nrow(df_train), nrow(df_test))

cat("number of training classes: \n", base::table(df_train$Class)/nrow(df_train))
cat("\n")
cat("number of test classes: \n", base::table(df_test$Class)/nrow(df_test))

X_train <- subset(df_train, select=-Class)
y_train <- df_train$Class
X_test <- subset(df_test, select=-Class) # exclude Class for prediction
y_test <- df_test$Class

?knn

model_knn <- knn(train=X_train,
                 test=X_test,
                 cl=y_train,  # class labels
                 k=3)
model_knn

conf_mat <- base::table(y_test, model_knn)
conf_mat

cat("Test accuracy: ", sum(diag(conf_mat))/sum(conf_mat))

?knn.cv

knn_loocv <- knn.cv(train=X_train, cl=y_train, k=3)
knn_loocv

conf_mat_cv <- base::table(y_train, knn_loocv)
conf_mat_cv
cat("LOOCV accuracy: ", sum(diag(conf_mat_cv)) / sum(conf_mat_cv))

SEED <- 2016
set.seed(SEED)
# create the training data 70% of the overall Sonar data.
in_train <- createDataPartition(Sonar$Class, p=0.7, list=FALSE) # create training indices
ndf_train <- Sonar[in_train, ]
ndf_test <- Sonar[-in_train, ]

?trainControl


ctrl <- trainControl(method="repeatedcv", number=5, repeats=2)

?expand.grid()

nn_grid <- expand.grid(k=c(1,3,5,7))
nn_grid
set.seed(SEED)
?train

best_knn <- train(Class~., data=ndf_train,
                  method="knn",
                  trControl=ctrl, 
                  preProcess = c("center", "scale"),  # standardize
                  tuneGrid=nn_grid)
best_knn

SEED <- 123 
set.seed(SEED) 
ctrl <- trainControl(method="repeatedcv", number=5, repeats=5) 
nn_grid <- expand.grid(k=c(1, 3, 5, 7)) 
best_knn_reduced <- train( Class~., data=ndf_train, method="knn", 
                           trControl=ctrl, preProcess=c("center", "scale","YeoJohnson"))
X_test <- subset(ndf_test, select=-Class) 
pred_reduced <- predict(best_knn_reduced, newdata=X_test, model="best") 
conf_mat_best_reduced <- confusionMatrix(ndf_test$Class, pred_reduced) 
conf_mat_best_reduced

View(Sonar)


# Load the data and view its structure
balance_scale <- read.csv("https://ibm.box.com/shared/static/684jzm7e6fbbssg87yc2v4dy53dgkdew.txt", sep = ",")
str(balance_scale)

# View the first few rows of the data using the head function
# Note: The raw data does not contain any column names
head(balance_scale)

# Add column names
colnames(balance_scale) <- c("Class_Name","Left_Weight", "Left_Distance", "Right_Weight", "Right_Distance")
head(balance_scale)
# Note: We do not need to standardize the data in this instance since the numerical data values lie on the same scale.
# Calculate the products and differences
Right_Product <- balance_scale[,4]*balance_scale[,5]
Left_Product <- balance_scale[,2]*balance_scale[,3]
Differences <- Right_Product-Left_Product
# Add columns for Right_Product, Left_Product and Differences
balance_scale$Right_Product <- Right_Product
balance_scale$Left_Product <- Left_Product
balance_scale$Differences <- Differences

# Use the sample function to create a vector of indices that will be used
set.seed(1234)
ind <- sample(2, nrow(balance_scale), replace=TRUE, prob=c(0.7, 0.3))


# Create the training and test data from the dataset using ind
bscale.train <- balance_scale[ind==1, 6:8]
bscale.test <- balance_scale[ind==2, 6:8]



# Create the target vectors for the training and test data from the dataset using ind
bscale.trainLabels <- balance_scale[ind==1, 1]
bscale.testLabels <- balance_scale[ind==2, 1]

# Use the knn command to make predictions on the Class_Name of the test data
knn_class <- knn(train = bscale.train, test = bscale.test, cl = bscale.trainLabels, k=3)

# Find the number of incorrectly classified points
correct <- which(knn_class == bscale.testLabels, arr.ind = TRUE)
incorrect <- which(knn_class != bscale.testLabels, arr.ind = TRUE)
cat("Number of incorrectly classified points:",length(incorrect),"\n")

# Find the proportion of correctly classified points
proportion_correct <- length(correct)/length(bscale.testLabels)
cat("Proportion of correctly classified points", proportion_correct,"\n")

# Run the knn regression using the kknn command
knn_reg <- kknn(formula = Differences ~ ., train=bscale.train, test=bscale.test, k=3)

# Find the number of incorrectly classified points
incorrect_reg <- which(knn_reg$fitted.values != bscale.test$Differences, arr.ind = TRUE)
cat("Number of incorrectly classified points:", length(incorrect_reg), "\n");

# Find the proportion of correctly classified points
correct_reg <- which(knn_reg$fitted.values == bscale.test$Differences, arr.ind = TRUE)
cat("Proportion of correctly classified points", length(correct_reg)/length(bscale.test$Differences), "\n")


# Display the first few rows of the regression estimates of the differences and their true values
head(cbind(knn_reg$fitted.values,bscale.test$Differences))

best_reg <- train.kknn(formula = Differences ~ ., data=bscale.train, kmax=8)
best_reg$best.parameters

# Run the knn regression again using the kknn command with k=2
knn_reg2 <- kknn(formula = Differences ~ ., train=bscale.train, test=bscale.test, k=2)

# Find the number of incorrectly classified points
incorrect_reg2 <- which(knn_reg2$fitted.values != bscale.test$Differences, arr.ind = TRUE)
cat("Number of incorrectly classified points:", length(incorrect_reg2),"\n")
# Find the proportion of correctly classified points
correct_reg2 <- which(knn_reg2$fitted.values == bscale.test$Differences, arr.ind = TRUE)
cat("Proportion of correctly classified points",length(correct_reg2)/length(bscale.test$Differences),"\n")

# Display the first few rows of the new regression estimates of the differences and their true values
head(cbind(knn_reg2$fitted.values,bscale.test$Differences))

# The kknn function can also be used for classification
knn_class2 <- kknn(formula = Class_Name ~ ., train=subset(balance_scale, select=c(Class_Name,Right_Product,Left_Product,Differences))[ind==1,], test=subset(balance_scale, select=c(Class_Name,Right_Product,Left_Product,Differences))[ind==2,], k=3)
# Find the number of incorrectly classified points
incorrect_class2 <- which(knn_class2$fitted.values != bscale.testLabels, arr.ind = TRUE)
cat("Number of incorrectly classified points:",length(incorrect_class2),"\n")
best_class <- train.kknn(formula = Class_Name ~ ., data=subset(balance_scale, select=c(Class_Name,Right_Product,Left_Product,Differences))[ind==1,], kmax=8)
best_class$best.parameters

# Using k=1
knn_class3 <- kknn(formula = Class_Name ~ ., train=subset(balance_scale, select=c(Class_Name,Right_Product,Left_Product,Differences))[ind==1,], test=subset(balance_scale, select=c(Class_Name,Right_Product,Left_Product,Differences))[ind==2,], k=1)
# Find the number of incorrectly classified points
incorrect_class3 <- which(knn_class3$fitted.values != bscale.testLabels, arr.ind = TRUE)
cat("Number of incorrectly classified points:", length(incorrect_class3),"\n")







