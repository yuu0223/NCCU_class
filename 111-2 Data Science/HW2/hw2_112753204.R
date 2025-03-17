#cmd存取
args = commandArgs(trailingOnly=TRUE)
#判斷flag數量 看有沒有missing flag
grep_args<-grep("--", args, value = T)
if (length(grep_args)!=4) {
  stop("missing flag", call.=FALSE)
}

#cmd解析
target_i<-which(args == "--target")
target<-args[target_i+1]

badthre_i<-which(args == "--badthre")
badthre<-args[badthre_i+1]

output_i<-which(args == "--output")
output<-args[output_i+1]

input_i<-which(args == "--input")
after_input<-grep(pattern = "-", x = c(args[(input_i+1):length(args)]))[1]
input<-array(c(args[(input_i+1):(input_i+after_input-1)]))

#input計算各指標
final_df<-data.frame()
for (file in input){
  #取出檔名
  input_file_name<-tail(strsplit(file, "/")[[1]],1)
  input_file_name<-strsplit(input_file_name, ".csv")[[1]]
  
  data<-read.csv(file = file, header = T, sep = ",")
  
  #將預測分數標記成bad/good
  for (i in 1:length(data$pred.score)){
    if (data$pred.score[i] > badthre){
      data$pred_reference[i]<-"bad"
    }
    else{
      data$pred_reference[i]<-"good"
    }
  }
  #Confusion Matrix 看target是good/bad再來決定誰是TP
  cm<-table(data$pred_reference, data$reference)
  if (target == "bad"){
    TP<-cm[1,1]
    FP<-cm[1,2]
    FN<-cm[2,1]
    TN<-cm[2,2]
    
  }else{
    TP<-cm[2,2]
    FP<-cm[2,1]
    TN<-cm[1,1]
    FN<-cm[1,2]
  }
  
  #sensitivity
  sensitivity<-round(TP/(TP+FN), digits = 2)
  
  #specificity
  specificity<-round(TN/(TN+FP), digits = 2)
  
  #F1-score
  F1<-round((2*TP)/(2*TP+FP+FN), digits = 2)
  
  #Log LikeliHood
  LLH<-round(sum(ifelse(data$reference=="bad", log(data$pred.score), log(1-data$pred.score))), digits = 2)
  #NULL's Log LikeliHood
  pNull<-sum(ifelse(data$reference=="bad",1,0))/dim(data)[[1]]
  LLH_NULL<-sum(ifelse(data$reference=="bad",1,0))*log(pNull) +sum(ifelse(data$reference=="bad",0,1))*log(1-pNull)
  
  #R squared
  S<-0
  R2<-round(1-(-2*(LLH-S))/(-2*(LLH_NULL-S)), digits = 2)
  
  #把所有input整理成一張table
  df<-data.frame(input = input_file_name,
                 sensitivity = sensitivity, 
                 specificity = specificity, 
                 F1_Score = F1,
                 LogLikeliHood = LLH,
                 pseudoR2 = R2)
  
  final_df<-rbind(final_df,df)
  
}
#找出最佳解
compare_method<-c("best",
                  final_df$input[which.max(final_df$sensitivity)],
                  final_df$input[which.max(final_df$specificity)],
                  final_df$input[which.max(final_df$F1_Score)],
                  final_df$input[which.max(final_df$LogLikeliHood)],
                  final_df$input[which.max(final_df$pseudoR2)])

final_df<-rbind(final_df,compare_method)
colnames(final_df)<-c("method","sensitivity","specificity","F1","logLikelihood","pseudoRsquared")

#output file
write.csv(final_df, file = output, fileEncoding = "UTF-8", row.names = FALSE, quote=FALSE)

