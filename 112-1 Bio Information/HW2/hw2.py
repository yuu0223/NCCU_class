import pandas as pd
import numpy as np
import subprocess

# 下載pysam套件
package_name = "pysam"
subprocess.check_call(["pip", "install", package_name])

import pysam


def read_fasta(input_path):
    fasta = pysam.FastaFile(input_path)
    seq_name = fasta.references

    seq_list = list()
    for name in seq_name:
        seq_list.append(fasta.fetch(name))

    # 將各seqences兩兩要比對的存成一組
    seq_compare_list = list()
    for i in range(len(seq_list)):
        seq1 = seq_list[i]
        for j in range(i + 1, len(seq_list)):
            seq2 = seq_list[j]
            seq_compare_list.append(np.array([seq1, seq2]))

    return seq_name, seq_compare_list


def read_PAM(pam_txt):
    pam_arr = []
    with open(pam_txt, "r") as file:
        for line in file:
            line = line.strip()

            # 將開頭為#的行數省略
            if not line.startswith("#"):
                row = line.split()
                pam_arr.append(row)

    # 第一行都需要一個空格才能與其他row相同長度
    pam_arr[0] = [""] + pam_arr[0]

    # 開始轉為df
    pam_df = pd.DataFrame(pam_arr)
    pam_df.columns = pam_df.iloc[0]
    pam_df = pam_df.set_index(pam_df.iloc[:, 0])
    pam_df = pam_df.iloc[1:, 1:]

    return pam_df


def calculate_SoP(input_path, score_path, gopen, gextend):
    seq_name, seq_compare_list = read_fasta(input_path)
    pam_df = read_PAM(score_path)

    score_list = []
    for i in range(len(seq_compare_list)):
        seq1 = seq_compare_list[i][0]
        seq2 = seq_compare_list[i][1]

        if len(seq1) == len(seq2):
            for num in range(len(seq1)):
                # 除了第0個，其他有"-"的判斷是extension/open
                if "-" in [seq1[num], seq2[num]] and num != 0:
                    extend = (seq1[num] == seq1[num - 1] == "-") | (
                        seq2[num] == seq2[num - 1] == "-"
                    )
                    if extend:
                        get_score = gextend
                    else:
                        get_score = gopen
                    score_list.append(get_score)
                # 第0個，有出現"-"的為open
                elif num == 0 and (seq1[num] == "-" or seq2[num] == "-"):
                    get_score = gopen
                    score_list.append(get_score)
                # 其餘，使用pam.txt的表進行分數計算
                else:
                    no_1 = seq1[num] not in pam_df.columns
                    no_2 = seq2[num] not in pam_df.columns

                    if no_1:
                        score_list.append(int(pam_df["*"][seq2[num]]))

                    elif no_2:
                        score_list.append(int(pam_df[seq1[num]]["*"]))

                    elif no_1 and no_2:
                        score_list.append(int(pam_df["*"]["*"]))

                    else:
                        score_list.append(int(pam_df[seq1[num]][seq2[num]]))

    return sum(score_list)
