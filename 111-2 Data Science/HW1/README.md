# Familiar with R basics & submit homework on github

#### Name: 陳品伃
#### Student ID: 112753204

## Description
### cmd
```R
Rscript hw1_yourID.R --input input1.csv --output output1.csv

Rscript hw1_yourID.R --output output1.csv --input input1.csv
```

Your R code should output and round the set name with maximum value of weight and height.

### Read an input file

Input data will have other numeric & category columns besides weight and height.

examples = input1.csv

### Output a summary file

Please follow the same format of the result.csv, i.e., round number into two digitals

example =  output1.csv

### Code for reference
```R
# Simple example of extracting input args
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: Rscript hw1_exam.R input", call.=FALSE)
} else if (length(args)==1) {
  i_f <- args[1] 
}
print(i_f)
```
## Score

### 10 testing data (90%)

```R
Rscript hw1_5566.R --input hw1/data/test.1.csv --output hw1/eval/test1/hw1_001.csv
Rscript hw1_5566.R --output hw1/eval/test2/hw1_002.csv --input hw1/data/test.2.csv
```
Correct answer gets 9 points of each testing data.

### Bonus (10%)

- Output format without “: 3 points
- Concise file name without path: 3 points
- Concise file name without .csv extension: 4 points

### Penalty: -2 points of each problem

- Can not detect missing --input/--ouptut flag
- Arguments order cannot change
- Wrong file name
- Wrong column name
- Not round number to 2 digitals

## Note
- **Please use R version 3.6.3**
- Please do not set working directory(setwd) in a fixed folder. For example,
```R
d <- read.csv("D://DataScience/hw1/example/output1.csv")
```
- **execution time: 1 hour maximum**
