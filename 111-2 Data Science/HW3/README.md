# Predict protein subcellular localization
![PredictProtein](/images/img1.png)

#### Name: 陳品伃
#### Student ID: 112753204

## Description
Perform k-fold cross-validation for protein subcellular localization problem.

### cmd
```R
Rscript hw3_studentID.R --fold k --input Archaeal_tfpssm.csv --output performance.csv
```

### K-fold
* Divide the data into k parts, the number of parts used by each data set
  * (training, validation, testing) = (*k*-2, 1, 1)
* The following shows the example of the 5-fold cross validation.

<br> 

![cross-validation](/images/img2.png)

### Input: Archaeal_tfpssm.csv

[📁 Archaeal_tfpssm.csv download](https://drive.google.com/file/d/1L-gv1dPaEonnaASeBtakePT1t3FJwPFI/view?usp=sharing)

This CSV doesn't contain a header. The information of columns as below:

* `V2`: labels of proteins
  * CP: Cytoplasmic
  * CW: Cell Wall
  * EC: Extracellular
  * IM: Inner membrane

* `V3 ~ V5602`: the gapped-dipeptide features of each protein

### Output format: performance.csv

* accuracy = *P*/*N*, average of *k*-fold cross-validation

set|training|validation|test
---|---|---|---
fold1|0.93|0.91|0.88
fold2|0.92|0.91|0.89
fold3|0.94|0.92|0.90
fold4|0.91|0.89|0.87
fold5|0.90|0.92|0.87
ave.|0.92|0.91|0.88


### Code for reference

```R
library('rpart')
# read input data
d <- read.csv(<Path to Archaeal_tfpssm.csv>, header = F)
# label to be predicted
levels(d[,2])
head(d[,5600:5603])
# select subset of the data
tmp <- d[c(seq(1,700,25), seq(700,800,5)),]
# model using decision tree
model <- rpart(V2 ~ V3 + V4 + V5600 + V5601 + V5602,
               data=tmp, control=rpart.control(maxdepth=4),
               method="class")
# make confusion matrix tabel
resultframe <- data.frame(truth=tmp$V2,
                          pred=predict(model, type="class"))
(rtab <- table(resultframe)) 
```

## Score
### 10 testing cmds (100%)
```R
Rscript hw3_studentID.R --fold 5 --input Archaeal_tfpssm.csv --output hw4/your_ID/output1.csv
...
Rscript hw3_studentID.R --fold 10 --input Archaeal_tfpssm.csv --output hw4/your_ID/output6.csv
```
* Each testing cmd gets 10 points.
* We will have different test data beyond Archaeal_tfpssm.csv. These data have the same columns, but the label ratio and the number of data will be different. The number of data will be less than Archaeal_tfpssm.csv.
* You must design your own test cases to ensure that the code can pass the test by TA 💪💪

### Bonus
* Round number to two decimal places: 2 points
* Output format without `"`: 2 points

## Note
- Please use R version 4
- Filename format of your code: `hw3_YourStudentID.R`
- Please do not set input/output in your local path or URL. Otherwise, your code will fail due to fixed path problem.
- Up to 100 points for this assignment. (10 testing cmds + Bonus)
- execution time: 5 minute maximum per cmd
- Only the `rpart` package can be used as a predictive model, no additional packages can be used
- If you cheat, you'll get zero for this assignment, and there's no chance to make up for it. 😠 
  - i.e., Generate the output data yourself or without kfold process or other cheating methods

## References
* Chang, J.-M. M. et al. (2013) [Efficient and interpretable prediction of protein functional classes by correspondence analysis and compact set relations](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0075542). *PLoS ONE* 8, e75542.
* Chang J-M, Su EC-Y, Lo A, Chiu H-S, Sung T-Y, & Hsu W-L (2008) [PSLDoc: Protein subcellular localization prediction based on gapped-dipeptides and probabilistic latent semantic analysis](https://onlinelibrary.wiley.com/doi/full/10.1002/prot.21944). *Proteins: Structure, Function, and Bioinformatics* 72(2):693-710.
