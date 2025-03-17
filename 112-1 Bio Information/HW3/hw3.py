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


def global_alignment(seq_compare_list, pam_df, gap):
    ### Load Data
    seq1 = "-" + str(seq_compare_list[0][0])
    seq2 = "-" + str(seq_compare_list[0][1])
    len_1, len_2 = len(seq1), len(seq2)

    ### 創建兩個matrix 一個紀錄最佳分數 一個紀錄腳步
    score_mat = np.zeros((len_1, len_2), dtype=int)
    where_mat = np.zeros((len_1, len_2), dtype=int)
    gap_score = gap

    ### Dynamic Programming
    for i in range(len_1):
        for j in range(len_2):
            # initial the matrix data
            if i == 0 or j == 0:
                if i == 0 and j != 0:
                    score_mat[i, j] = score_mat[i, j - 1] + gap_score
                    where_mat[i, j] = 100

                elif i != 0 and j == 0:
                    score_mat[i, j] = score_mat[i - 1, j] + gap_score
                    where_mat[i, j] = 100

                else:
                    score_mat[i, j] = 0
                    where_mat[i, j] = 100

            # start to calculate
            else:
                alignment = score_mat[i - 1, j - 1] + int(pam_df[seq1[i]][seq2[j]])
                deletion = score_mat[i, j - 1] + gap_score
                insertion = score_mat[i - 1, j] + gap_score

                select = [alignment, deletion, insertion]
                select_index = np.argmax(select)

                # 將matrix分數更新 並儲存他是從哪一個方法來的
                score_mat[i, j] = select[select_index]
                where_mat[i, j] = int(
                    select_index
                )  # 0=alignment/1=deletion/2=insertion

    ### Trace Back
    road = []

    i = len_1 - 1
    j = len_2 - 1
    while i > 0 or j > 0:
        current = where_mat[i, j]
        road.append(current)

        if current == 0:
            i -= 1
            j -= 1
        elif current == 1:
            j -= 1
        elif current == 2:
            i -= 1

    road.reverse()

    ### 整理Sequences結果
    seq1 = str(seq_compare_list[0][0])
    seq2 = str(seq_compare_list[0][1])

    global_seq1, global_seq2 = "", ""
    s1, s2 = 0, 0
    for current_road in road:
        if current_road == 0:
            global_seq1 += seq1[s1]
            global_seq2 += seq2[s2]
            s1 += 1
            s2 += 1
        elif current_road == 1:
            global_seq1 += "-"
            global_seq2 += seq2[s2]
            s2 += 1
        elif current_road == 2:
            global_seq1 += seq1[s1]
            global_seq2 += "-"
            s1 += 1

    return global_seq1, global_seq2


def local_alignment(seq_compare_list, pam_df, gap):
    ### Load Data
    seq1 = "-" + str(seq_compare_list[0][0])
    seq2 = "-" + str(seq_compare_list[0][1])
    len_1, len_2 = len(seq1), len(seq2)

    ### 創建兩個matrix 一個紀錄最佳分數 一個紀錄腳步
    score_mat = np.zeros((len_1, len_2), dtype=int)
    where_mat = np.zeros((len_1, len_2), dtype=int)
    gap_score = gap

    ### Dynamic Programming
    for i in range(len_1):
        for j in range(len_2):
            # initial the matrix data
            if i == 0 or j == 0:
                if i == 0 and j != 0:
                    score_mat[i, j] = (
                        score_mat[i, j - 1] + gap_score
                        if (score_mat[i, j - 1] + gap_score > 0)
                        else 0
                    )
                    where_mat[i, j] = 100

                elif i != 0 and j == 0:
                    score_mat[i, j] = (
                        score_mat[i - 1, j] + gap_score
                        if (score_mat[i - 1, j] + gap_score > 0)
                        else 0
                    )
                    where_mat[i, j] = 100

                else:
                    score_mat[i, j] = 0
                    where_mat[i, j] = 100

            # start to calculate
            else:
                alignment = score_mat[i - 1, j - 1] + int(pam_df[seq1[i]][seq2[j]])
                deletion = score_mat[i, j - 1] + gap_score
                insertion = score_mat[i - 1, j] + gap_score
                restart = 0

                select = [alignment, deletion, insertion, restart]
                select_index = np.argmax(select)

                # 將matrix分數更新 並儲存他是從哪一個方法來的
                score_mat[i, j] = select[select_index]
                where_mat[i, j] = int(
                    select_index
                )  # 0=alignment/1=deletion/2=insertion/3=restart

    ### Find the max score's index (maybe not the only one)
    max_where = np.where(score_mat == np.max(score_mat))
    max_count = len(max_where[0])

    trace_index = []
    for i in range(max_count):
        trace_index.append([max_where[0][i], max_where[1][i]])
        i += 1

    ### Trace Back
    roads_get = []

    for i, j in trace_index:
        road = []
        end_i, end_j = i, j

        while i > 0 or j > 0:
            current = where_mat[i, j]
            road.append(current)

            if current == 0:
                i -= 1
                j -= 1
            elif current == 1:
                j -= 1
            elif current == 2:
                i -= 1
            else:
                break

        road.reverse()
        roads_get.append([i, j, end_i, end_j, road])

    ### 整理Sequences結果
    seq1 = str(seq_compare_list[0][0])
    seq2 = str(seq_compare_list[0][1])

    # 選出最長的序列回傳
    road_select_max = []
    max_i, max_j, max_tail_i, max_tail_j = 0, 0, 0, 0
    max_road = [0]
    for road in roads_get:
        current_i, current_j, current_tail_i, current_tail_j, current_road = (
            road[0],
            road[1],
            road[2],
            road[3],
            road[4],
        )

        if len(current_road) > len(max_road):
            max_road = current_road
            max_i, max_j = current_i, current_j
            max_tail_i, max_tail_j = current_tail_i, current_tail_j

    road_select_max.append([max_i, max_j, max_tail_i, max_tail_j, max_road])

    # 最長的序列開始返回原始sequences
    local_seq1, local_seq2 = "", ""
    s1, s2 = road_select_max[0][0], road_select_max[0][1]
    road = road_select_max[0][4]

    for current_road in road:
        if current_road == 0 or current_road == 3:
            local_seq1 += seq1[s1]
            local_seq2 += seq2[s2]
            s1 += 1
            s2 += 1
        elif current_road == 1:
            local_seq1 += "-"
            local_seq2 += seq2[s2]
            s2 += 1
        elif current_road == 2:
            local_seq1 += seq1[s1]
            local_seq2 += "-"
            s1 += 1

    return local_seq1, local_seq2


def output_fasta(seq1, seq2, seq_name, output_path):
    seqs = [seq1, seq2]
    output_seq = list(zip(seq_name, seqs))

    with open(output_path, "w") as fasta_file:
        for seq_name, seq in output_seq:
            fasta_file.write(f">{seq_name}\n")
            fasta_file.write(seq + "\n")


def alignment(input_path, score_path, output_path, aln, gap):
    seq_name, seq_compare_list = read_fasta(input_path)
    pam_df = read_PAM(score_path)

    if not gap:
        gap = -10

    if aln == "global":
        seq1, seq2 = global_alignment(seq_compare_list, pam_df, gap)
        output_fasta(seq1, seq2, seq_name, output_path)

    elif aln == "local":
        seq1, seq2 = local_alignment(seq_compare_list, pam_df, gap)
        output_fasta(seq1, seq2, seq_name, output_path)

    else:
        print("There are some errors in your parameters.")
