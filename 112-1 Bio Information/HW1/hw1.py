import pandas as pd
import numpy as np


def generate_pam(x, input_path, output_path):
    # 讀取input檔案 mut.txt 並 轉換成matrix
    mut = pd.read_csv(input_path, delim_whitespace=True, index_col=0, skiprows=[0])
    mut_matrix = mut.values

    # 將mutation/10000
    mut_matrix = mut_matrix.astype(float) / 10000.0
    # 將矩陣相乘x次，PAM250就是乘以250次
    result_matrix = np.linalg.matrix_power(mut_matrix, x)

    # 對照normalized frequency table去除以相對應的f
    result_matrix[0] = result_matrix[0].astype(float) / 0.087  # A
    result_matrix[1] = result_matrix[1].astype(float) / 0.041  # R
    result_matrix[2] = result_matrix[2].astype(float) / 0.040  # N
    result_matrix[3] = result_matrix[3].astype(float) / 0.047  # D
    result_matrix[4] = result_matrix[4].astype(float) / 0.033  # C
    result_matrix[5] = result_matrix[5].astype(float) / 0.038  # Q
    result_matrix[6] = result_matrix[6].astype(float) / 0.050  # E
    result_matrix[7] = result_matrix[7].astype(float) / 0.089  # G
    result_matrix[8] = result_matrix[8].astype(float) / 0.034  # H
    result_matrix[9] = result_matrix[9].astype(float) / 0.037  # I
    result_matrix[10] = result_matrix[10].astype(float) / 0.085  # L
    result_matrix[11] = result_matrix[11].astype(float) / 0.081  # K
    result_matrix[12] = result_matrix[12].astype(float) / 0.015  # M
    result_matrix[13] = result_matrix[13].astype(float) / 0.040  # F
    result_matrix[14] = result_matrix[14].astype(float) / 0.051  # P
    result_matrix[15] = result_matrix[15].astype(float) / 0.070  # S
    result_matrix[16] = result_matrix[16].astype(float) / 0.058  # T
    result_matrix[17] = result_matrix[17].astype(float) / 0.010  # W
    result_matrix[18] = result_matrix[18].astype(float) / 0.030  # Y
    result_matrix[19] = result_matrix[19].astype(float) / 0.065  # V

    # 最後乘上10以及以10為底的log（取對數）
    result_matrix = np.round(10 * np.log10(result_matrix)).astype(int)

    # 整理表格＆輸出
    result_matrix_df = pd.DataFrame(
        result_matrix,
        index=[
            "A",
            "R",
            "N",
            "D",
            "C",
            "Q",
            "E",
            "G",
            "H",
            "I",
            "L",
            "K",
            "M",
            "F",
            "P",
            "S",
            "T",
            "W",
            "Y",
            "V",
        ],
        columns=[
            "A",
            "R",
            "N",
            "D",
            "C",
            "Q",
            "E",
            "G",
            "H",
            "I",
            "L",
            "K",
            "M",
            "F",
            "P",
            "S",
            "T",
            "W",
            "Y",
            "V",
        ],
    )
    result_matrix_df.to_csv(output_path, sep=" ", index=True)
