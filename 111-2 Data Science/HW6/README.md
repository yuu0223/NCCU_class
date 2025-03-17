# GiveMeSomeCredit
#### Name: 陳品伃
#### Student ID: 112753204
#### Snapshot:
![](image/112753204_kaggle_snapshot.png)

## Description
Start here! Predict the probability that somebody will experience financial distress in the next two years (column=SeriousDlqin2yrs).

[GiveMeSomeCredit](https://www.kaggle.com/c/GiveMeSomeCredit/overview)

### Steps
1. You could try different models, parameter and data processing methos.
2. Use the best model to make final predictions on test dataset,  and output to `studentID.csv` (--predict)
   
```R
Rscript hw6_studentID.R --train Data/train.csv --test Data/test.csv --predict studentID.csv
```

3. Submit `studentID.csv` to Kaggle system.
4. Make a snapshot of your submissions and update in your README.md.  Like the example above.

## Score
### Private score on kaggle
- 0.86 ~ : 100 points
- 0.80 ~ 0.86: 95 points
- 0.70 ~ 0.80: 90 points
- 0.60 ~ 0.70: 85 points

_base on your best score round down second decimal place. (ex. 0.84724 -> 0.84)_

### Penalty
- -20 points: The results uploaded by the TA do not match the snapshot.

## Note
- Please use R version 4
- execution time: 30 minute maximum
- Filename format of your code: `hw6_YourStudentID.R`
- TA will check your code to make sure the final `studentID.csv` is produced by your code, and upload to Kaggle for verification.
- `sampleEntry.csv` for reference only,  cannot be used directly or indirectly.
- If you cheat, you'll get zero for this assignment, and there's no chance to make up for it. 😠
