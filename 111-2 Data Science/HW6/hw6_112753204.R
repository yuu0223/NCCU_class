library(factoextra)
library(cluster)
library(zoo)
library(caret)
library(randomForest)
library(tibble)


#cmd存取
args = commandArgs(trailingOnly=TRUE)
#判斷flag數量 看有沒有missing flag
if (length(args)==0) {
  stop("missing flag", call.=FALSE)
}

#cmd解析
i<-1
while(i < length(args))
{
  print(args[i])
  if(args[i] == "--train"){
    train<-args[i+1]
    i<-i+2
  }else if(args[i] == "--test"){
    test<-args[i+1]
    i<-i+2
  }else if(args[i] == "--predict"){
    output<-args[i+1]
    i<-i+2
  }else{
    stop(paste("Unknown flag", args[i]), call.=FALSE)
  }
}

###read file
df <- read.csv(train, header = TRUE)
test_df <- read.csv(test, header = TRUE)

###rename columns
colnames(df)[1] <- "ID"
colnames(df)[2] <- "target"
colnames(test_df)[1] <- "ID"
colnames(test_df)[2] <- "target"

# Delete NA
df<- df[complete.cases(df), ]
test_df[is.na(test_df)] <- -999


### Training Start
train_data <- df
test_data <- test_df

### Random Forest
model <- randomForest(target ~ ., data = train_data, ntree = 70)
predictions <- predict(model, newdata = test_data)

### Normalization
probabilities <- (predictions - min(predictions)) / (max(predictions) - min(predictions))

### Output
final_df <- data.frame(probabilities)
colnames(final_df) <- "Probability"
final_df <- tibble::rownames_to_column(final_df, var = "Id")
write.csv(final_df, file = output, row.names=F)


