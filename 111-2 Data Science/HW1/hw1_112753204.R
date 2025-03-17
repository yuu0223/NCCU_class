#cmd解析
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: Rscript hw1_exam.R no input", call.=FALSE)
} else if (length(args)==4) {
  if (args[[1]] == "--input"){
    input_file <- args[[2]]
    output_file <- args[[4]]
  }
  else if (args[[1]] == "--output"){
    input_file <- args[[4]]
    output_file <- args[[2]]
  }
} else {
    stop("USAGE: Rscript hw1_exam.R input/output not correct.", call.=FALSE)
}

#前面讀取的檔案有路徑，取出檔案名稱
input_file_name<-tail(strsplit(input_file, "/")[[1]],1)

#讀取input檔案
test <- read.csv(input_file,
                 sep = ",",
                 header = T)

weight_max <- round(max(test[['weight']]), digits = 2)
height_max <- round(max(test[['height']]), digits = 2)

#把input file的副檔名.csv拿掉
set_name <- strsplit(input_file_name, ".csv")[[1]]
df <- data.frame(set = set_name, weight = weight_max, height = height_max)

#將cmd中所輸入的output file帶入這邊
write.csv(df, file = output_file, fileEncoding = "UTF-8", row.names = FALSE, quote=FALSE)


