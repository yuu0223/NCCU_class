# Evalute a mode for predicting loans that are about to default

#### Name: 陳品伃
#### Student ID: 112753204

## Description

### cmd
```R
Rscript hw2_yourID.R --target bad/good --badthre <threshold> --input meth1 meth2 ... methx --output result.csv

```

* Read in multiple files
* Positive case defined by “--target” option
* Threshold defined by “--badthre” option
* yourID should be your student ID number, i.e., hw2_106769999.R
* hw2_example.R is for reference only.

### Inputs format
* the last column, pred.score, is the predicted probability of "bad loan".

#### examples: `examples/method1.csv`

|persons|reference|pred.score|
|-------|---------|----------|
|person1|bad      |0.807018548483029|
|person2|bad      |0.740809247596189|
|person3|bad      |0.0944965328089893|
|person4|good     |0.148418645840138|

### Output format
* You should write your own function to calculate metrics.
* Find out which method performs the best regarding the metric. 
* pseudo *R*<sup>2</sup> = 1 - deviance(model)/deviance(null model) for *S*=0.

#### examples: `examples/output1.csv`

|method |sensitivity|specificity|F1     |logLikelihood|pseudoR2|
|-------|-----------|-----------|-------|-------------|--------|
|method1|0.91       |0.96       |0.85   |-132         |0.79    |
|method2|0.99       |0.98       |0.86   |-112         |0.70    |
|best   |method2    |method2    |method2|method2      |method1 |


### Null model
* as being "the obvious guess": always output a constant number no matter the input
* always return the proportion of "bad" loans

### Examples of test command

```R
Rscript hw2_yourID.R --target bad --badthre 0.5 --input examples/method1.csv examples/method2.csv --output examples/output1.csv
Rscript hw2_yourID.R --target bad --badthre 0.4 --input examples/method1.csv examples/method3.csv examples/method5.csv --output examples/output2.csv
Rscript hw2_yourID.R --target good --badthre 0.6 --input examples/method2.csv examples/method4.csv examples/method6.csv --output examples/output3.csv
```

## Scores

### 5 testing data (90%)

Correct answer gets 18 points of each testing data.

### Bonus (10%)

- Output format without `“`: 3 points
- Number in 2 digitals : 3 points
- Concise file name without path in the `set` column: 4 points

*Note: The above regulations are for the format specification in the output data*

### Penalty: -2 points of each problem
- Can not detect missing --input/--ouptut flag

## Note
- **Please use R version 4**
- **Filename format of your code: `hw2_YourStudentID.R`**
- **Please do not set input/output in your local path or URL. Otherwise, your code will fail due to fixed path problem.** 
- **execution time: 1 hour maximum**
