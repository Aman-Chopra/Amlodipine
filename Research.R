#Tips before running :

#This script assumes that the data is present in a file named "Pedal.xlsx" 
#present in a folr named "R_Files" located in "~/Documents/R_Files"
#If this script is to be run on any other location or file :
#1 . Set the parameters to setwd() as the location of the file (dataset in xlsx form)
#2 . Rename the file name parameter in read.xlsx() to the dataset name with file extension (along with quotes to notify it as a string). 

#PRE-REQUISITE LIBRARIES : 

print("Installing Required Libraries")
#install.packages("xlsx")
#install.packages("corrplot")
#install.packages("e1071")
#install.packages("caret")
#install.packages("randomForest")
#install.packages("gmodels")
#install.packages("class")
#install.packages("Hmisc")
#install.packages("Boruta")
#install.packages("mlbench")
#install.packages("gbm")
#install.packages("GA");

#Importing dataset
library("xlsx")
setwd("C:/Users/win-8/Documents/fwdphilipshackaboutteam27")
dataset <- read.xlsx("Pedal.xlsx",1)

#Cleaning the data
#Removing non-numeric data of names because it is not important
#dataset <- subset(dataset, select = -c(Name,Dose Mg) )
dataset <- dataset[, -c(1,31,32,33,34,35)]
dataset <- head(dataset,-4)

#Applying preprocessing steps
#Calculating mean value to impute
#Some values of only three out of 36 attributes were missing
library(Hmisc)
LVEF <- with(dataset, impute(LVEF, mean))
Fractional <- with(dataset, impute(Fractional_Na, mean))
IVF <- with(dataset, impute(IVC, mean))

#Shuffling the dataset for unbiased results
aman<-dataset
set.seed(98)
random<- runif(nrow(dataset))
dataset<-aman[order(random),]

# Summary of the dataset
print(summary(dataset))
str(dataset)

#Correlation of variables amongst themselves and to Groups
library(corrplot)
print(cor(dataset, use="complete.obs", method="spearman"))
COR <- cor(as.matrix(dataset[,32]), as.matrix(dataset[,-1]), use="pairwise", method="spearman")
print(COR)
corrplot(COR, method="number")

#Finding the highly correlated variables
library(mlbench)
library(caret)
correlationMatrix <- cor(dataset[,1:32])
print(correlationMatrix)
# find attributes that are highly corrected (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.5)
print(highlyCorrelated)

# Dividing the dataset based on three different groups to see their variations
group_1 <- subset(dataset , Groups == 1)	
group_2 <- subset(dataset , Groups == 2)
group_3 <- subset(dataset , Groups == 3)

# See the initial pattern of clustering
library(e1071)
plot(dataset$Weight_Kgs,dataset$Plasma.creatinine,col=dataset$Groups)

#Splitting the data into training and test data
col<-c(names(dataset))
s <-sample(206,124)
a_train <- dataset[s,col]
a_test <- dataset[-s,col]

#Classification Algorithms

#1. Applying support vector Machine Classification Algorithm

svmfit <- svm(Groups ~ ., data=a_train , kernel="radial",cost=2,type="C-classification")
print(svmfit)

#cross validation using tune function, optimal cost = 2 and Gamma = 0.0277, kernel = radial
tune.out = tune(svm, Groups~., data = a_train, kernel = "radial", ranges = list(cost = c(.01,.1,1,2,2.5,3,10,100)))
summary(tune.out)

#Predict the class of test data
pred<-predict(svmfit,a_test[,col])
print(pred)

#Stroring the result in table
table<-table(pred,a_test[,1])

#Finding the accuracy
mean(pred==a_test[,1])

#Confidence interval and p-value
library(caret) 
confusionMatrix(table)

#10 fold cross validation using svm
library(caret)
ctrl <- trainControl(method = "cv", savePred=T, classProb=T)
mod <- train(Groups~., data=dataset, method = "svmLinear", trControl = ctrl)
print(mod$pred)
head(mod$pred)
print(mod)

#2. Applying K-Nearest Neighbours
#Normalizing data to apply knn algorithm
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x))) }
prc_n <- as.data.frame(lapply(dataset[,c(2:32)], normalize))
prc_train <- prc_n[1:124,]
prc_test <- prc_n[125:206,]
prc_train_target <- dataset[1:124,1]
prc_test_target <- dataset[125:206,1]
require(class)
ml<-knn(train=prc_train,test=prc_test,cl=prc_train_target,k=17)
print(ml)
tab<-table(prc_test_target,ml)
library(caret) 
confusionMatrix(tab)

#Feature importance

#1. Using Random-Forest
library(randomForest)
fit=randomForest(factor(Groups)~., data=dataset)
library(caret) 

varImp(fit)
varImpPlot(fit,type=2)

#2. Using Boruta
library(Boruta)
set.seed(134)
boruta.train <- Boruta(Groups~., data = dataset, doTrace = 2)
print(boruta.train)
final.boruta <- TentativeRoughFix(boruta.train)
print(final.boruta)
getSelectedAttributes(final.boruta, withTentative = F)
boruta.df <- attStats(final.boruta)
class(boruta.df)
print(boruta.df)

#3.Using lvq method
data5 <- read.xlsx("TRY.xlsx",1)
data5<-head(data5,-4)
data5 <- subset(data5, select = -c(1,31,32,33,34,35) )
control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
model <- train(Groups~., data=data5, method="lvq", preProcess="scale", trControl=control)
# estimate variable importance
importance <- varImp(model, scale=FALSE)
# summarize importance
print(importance)
# plot importance
plot(importance)

#4.Using recurssive feature elimination
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
results <- rfe(data5[,2:32], data5[,1], sizes=c(2:32), rfeControl=control)
# summarize the results
print(results)
# list the chosen features
predictors(results)
# plot the results
plot(results, type=c("g", "o"))

#Plotting Histograms to see the trend of important attributes with Groups

hist(group_3$Plasma.renin)
hist(group_1$Plasma.renin)
hist(group_3$SBP_3avg )
hist(group_1$SBP_3avg )

#5. Genetic algorithm

library(caret)
library(doParallel) # parallel processing
library(dplyr) # Used by caret
library(pROC) # plot the ROC curve

set.seed(10)
dim(data5)
head(data5,2)

trainIndex <- createDataPartition(dataset$Groups,p=.9,list=FALSE)
trainData <- dataset[trainIndex,]
testData <- dataset[-trainIndex,]

head(trainData,2)
head(testData,2)

summary(trainData)
summary(testData)

str(trainData)
trainX <-trainData[,-1] # Create training feature data frame
testX <- testData[,-1] # Create test feature data frame 
y=trainData$Groups # Target variable for training
str(trainX)

registerDoParallel(4) # Registrer a parallel backend for train
getDoParWorkers() # check that there are 4 workers

ga_ctrl <- gafsControl(functions = rfGA, # Assess fitness with RF
                       method = "cv",    # 10 fold cross validation
                       genParallel=TRUE, # Use parallel programming
                       allowParallel = TRUE)

set.seed(10)
lev <- c(1,2,3)     # Set the levels

system.time(rf_ga3 <- gafs(x = trainX, y = y,
                           iters = 100, # 100 generations of algorithm
                           popSize = 20, # population size for each generation
                           levels = lev,
                           gafsControl = ga_ctrl))

plot(rf_ga3) + theme_bw()




# Applying SVM on Features selected by GA
svmfitt <- svm(Groups ~ Gender+Weight_Kgs+Dose.Mg+Vasopressin.+Urinary.creatinine+DM+Plasma.creatinine+Plasma.renin.+Hyperlipidaemia+Plasma.sodium+Total.Proteins+Family_HisCAD+Fractional_Na+IVC+LVEF, data=a_train , kernel="radial",cost=2,type="C-classification")


svmfitt <- svm(Groups ~ osmolality+Proteinuria_24Hrs+Albumin+QTc+DBP_3avg+QTinterval_ms+Dose.Mg+Vasopressin.+Plasma.renin.+Fractional_Na+IVC+LVEF+SBP_3avg+VAS+VMA, data=a_train , kernel="radial",cost=2,type="C-classification")
svmfitt <- svm(Groups ~ DBP_3avg+Dose.Mg+Plasma.renin.+Fractional_Na+IVC+LVEF+SBP_3avg+VAS, data=a_train , kernel="radial",cost=2,type="C-classification")
print(svmfitt)
pred<-predict(svmfitt,a_test[,col])
print(pred)
table<-table(pred,a_test[,1])
plot(pred)
#Finding the accuracy
mean(pred==a_test[,1])
library(caret) 
confusionMatrix(table)



