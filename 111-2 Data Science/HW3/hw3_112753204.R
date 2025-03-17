library('rpart')

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
  if(args[i] == "--fold"){
    fold<-args[i+1]
    k<-as.numeric(fold)
    i<-i+1
  }else if(args[i] == "--input"){
    j<-grep("-", c(args[(i+1):length(args)], "-"))[1]
    input<-args[(i+1):(i+j-1)]
    i<-i+j-1
  }else if(args[i] == "--output"){
    output<-args[i+1]
    i<-i+1
  }else{
    stop(paste("Unknown flag", args[i]), call.=FALSE)
  }
  i<-i+1
}


if (k <= 2){
  stop("k-fold method k must large than 2.", call.=FALSE)
}

#開始
data<-read.csv(input, header = FALSE, sep = ",")

set.seed(101)
#shuffled 不然後面rpart會報錯
shuffled_data= data[sample(1:nrow(data)), ]


# PCA
new_data<-log(shuffled_data[,3:5602])
pca_data<-prcomp(new_data, center = TRUE, scale. = TRUE)
pred_pca <- predict(pca_data, newdata = new_data)
pca_data <- data.frame(pred_pca, label = shuffled_data$V2)
# 切割資料集
folds<-cut(seq(1, nrow(pca_data)), breaks = k, labels = FALSE)



fold <- list()
out_train <- list()
out_vali <- list()
out_test <- list()
#training test & validation
for (a in 1:k){
  
  if (a == k){
    
    vali_label<-which(folds == a)
    test_label<-which(folds == 1)
    
    validation<-pca_data[vali_label,]
    test_data<-pca_data[test_label,]
    train_data<-pca_data[-c(vali_label,test_label),]
    
  }
  else{
    
    vali_label<-which(folds == a)
    test_label<-which(folds == a+1)
    
    validation<-pca_data[vali_label,]
    test_data<-pca_data[test_label,]
    train_data<-pca_data[-c(vali_label,test_label),]
    
  }
  model <- rpart(label ~ ., data=train_data, method = "class")
  
  #validation
  pred_vali <- predict(model, newdata=validation, type="class")
  count_vali <- table(validation$label,pred_vali)
  acc_vali <- sum(diag(count_vali))/sum(count_vali)
  
  #training
  pred_train <- predict(model, newdata=train_data, type="class")
  count_train <- table(train_data$label,pred_train)
  acc_train <- sum(diag(count_train))/sum(count_train)
  
  #testing
  pred_test <- predict(model, newdata=test_data, type="class")
  count_test <- table(test_data$label,pred_test)
  acc_test <- sum(diag(count_test))/sum(count_test)
  
  #整理成好幾個list
  fold<-append(fold, paste("fold", a, sep=''))
  out_train<-append(out_train, round(acc_train, digits = 2))
  out_vali<-append(out_vali, round(acc_vali, digits = 2))
  out_test<-append(out_test, round(acc_test, digits = 2))
  
}

df<- data.frame(fold = unlist(fold),
                train = unlist(out_train),
                vali = unlist(out_vali),
                test = unlist(out_test))

colnames(df)<-c("set","training","validation","test")

#計算taining、validation、testing的平均accuracy
train_sum =0
vali_sum =0
test_sum=0
for (num in 1:nrow(df)){
  train_sum = train_sum+df$training[num]
  vali_sum = vali_sum+df$validation[num]
  test_sum = test_sum+df$test[num]
}

average<-c("ave.", 
           round(train_sum/nrow(df), digits = 2), 
           round(vali_sum/nrow(df), digits = 2),
           round(test_sum/nrow(df), digits = 2))

#最後整合成一個df
df<- rbind(df, average)
write.csv(df, file = output, row.names = FALSE)


