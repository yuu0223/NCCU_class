## 電影情緒正負面分類任務

### HW1-3
- [HW1 - Naive Bayes Classifier](https://github.com/yuu0223/NCCU_class/blob/main/112-2%20NLP/HW1/HW1_Naive_Bayes_Classifier.ipynb)
- [HW2 - Simple NN Classifier](https://github.com/yuu0223/NCCU_class/blob/main/112-2%20NLP/HW2/HW2_Simple_NN_Classifier.ipynb)
- [HW3 - Transformer](https://github.com/yuu0223/NCCU_class/blob/main/112-2%20NLP/HW3/HW3_Transformer.ipynb)

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

---
> Delete Stop Words
- **Method 1: nltk.corpus**
  
  ```Python
  import nltk
  from nltk.corpus import stopwords
  
  nltk.download('stopwords')
  stop_words = set(stopwords.words('english'))
  ```

- **Method 2: spacy**
  
  ```Python
  import spacy
  nltk.download('wordnet')
  
  nlp = spacy.load('en_core_web_sm')
  stop_words = spacy.lang.en.stop_words.STOP_WORDS
  ```
- **額外加入標點符號**
  ```Python
  sign_list = [";", ",", ".", "~", ":", "-", "(", ")", "%", "#", "$", "!", "/", "?", "=", "+", "&", "--", "'", '"', '`']
  for sign in sign_list:
      self.stop_words.add(sign)
  ```
---
> Stemmed 詞幹提取 & Lemmatized 詞形還原
  ```Python
  from nltk.stem import PorterStemmer, WordNetLemmatizer
  
  stemmer = PorterStemmer()
  lemmatizer = WordNetLemmatizer()
  
  tokens = [self.lemmatizer.lemmatize(token) for token in tokens]
  tokens = [self.stemmer.stem(token) for token in tokens]
  ```
---
> 忽略大小寫 (全部轉換為小寫)
  ```Python
  tokens = [token.lower() for token in tokens]
  ```
---
> Word Frequency 詞頻計算

  計算每個單詞在這份文件中總共出現的次數，並將次數少的字詞刪除。
  
  ```Python
  features = set()
  min_word_freq = 32 #刪除出現未超過32次的字詞
  features = {word for word in features if word_counts[word] >= min_word_freq }
  ```

---
> N-gram

  但這部分仍有個問題要解決，這個方法的 N-gram 只能處理還未斷詞的句子，使用```join```的話會切成一個一個的字母。
  ```Python
  def generate_ngrams(self, tokens):
      ngram_range = self.ngram_range
      ngrams = []
      for n in range(ngram_range[0], ngram_range[1] + 1):
          for i in range(len(tokens) - n + 1):
              ngram = ' '.join(tokens[i:i+n])
              ngrams.append(ngram)
      return ngrams
  
  def train():
  ...
  ngrams = self.generate_ngrams(tokens)
      for ngram in ngrams:
          self.features.add(ngram)
          self.class_feature_counts[label][ngram] += 1
  
  ```

### Colab Link
- [Google Colab](https://colab.research.google.com/drive/1zYN9doEg8is8gQkc3ewi39L1e-QiKaIz?usp=sharing)

