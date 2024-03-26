## HW1 - Naive Bayes Classifier 電影情緒正負面分類任務
### Description
使用 Naive Bayes Classifier 來進行電影情緒正負面分類任務，Dataset 使用 **nltk** 的 **movie_reviews**。
目標是希望能將電影分類提升到0.85以上。

### Limit
請修改 Code Lab 提供的 Naive Bayes Classifier，使得電影分類的效能可以提升到 0.85 以上，且越高越好。
注意，只能修改 NaiveBayesClassifier 類別的內容，不要更動 Code Lab 其他部份。此外，必須手動修改程式，不能使用現成的分類器，也不能偷渡非 Naive Bayes 精神的分類器（如 decision tree、logistic regression、neural network 等）。

在 Naive Bayes 模型的基礎上，可能修改的方向：
- word tokenization的方式
- 是否要合併 tokens（忽略大小寫、stemming、lemmatization…）
- 是否要考慮詞彙出現次數，而不只是有無出現
- 是否排除太少出現的詞、不重要的詞？
- 不同的 Additive smoothing 的 k 值
- 不同的 Smoothing 的方式 ( https://dl.acm.org/doi/pdf/10.1145/2187980.2188169 )
- 自行發揮創意

### Note
紀錄在作業上有嘗試過的方法，但未必有用在作業內。
> Delete Stop Words
- **Method 1: nltk.corpus**
  ```Python
  import nltk
  from nltk.corpus import stopwords
  
  nltk.download('stopwords')
  stop_words = set(stopwords.words('english'))
  ```
