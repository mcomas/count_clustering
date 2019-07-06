#install.packages('DirichletReg')
library(DirichletReg)
inputData <- ArcticLake
set.seed(100)
train <- sample (1:nrow (inputData), round (0.7*nrow (inputData))) 
train
inputData_train <- inputData [train, ] # training Data
inputData_test <- inputData [-train, ] # test Data
inputData$Y <- DR_data (inputData[,1:3])  # prepare the Y's
inputData_train$Y <- DR_data (inputData_train[,1:3])
inputData_test$Y <- DR_data (inputData_test[,1:3])


res1 <- DirichReg(Y ~ depth + I(depth^2), inputData_train)
res2 <- DirichReg(Y ~ depth + I(depth^2) | depth, inputData_train, model="alternative") 
summary(res1)
summary(res2)


predict(res1)
resid(res1)

predicted_res1 <- predict(res1, inputData_test)
predicted_res2 <- predict(res2, inputData_test)

plot(DR_data(predicted_res1))
plot(DR_data(predicted_res2))
plot(inputData$Y)


############
####
# source("https://bioconductor.org/biocLite.R")
# biocLite('DirichletMultinomial')
