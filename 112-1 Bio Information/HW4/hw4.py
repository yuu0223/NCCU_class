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


def global_alignment(seq_compare_list, pam_df, gop, gep):
    ### Load Data
    seq1 = "-" + str(seq_compare_list[0][0])
    seq2 = "-" + str(seq_compare_list[0][1])
    len_1, len_2 = len(seq1), len(seq2)

    ### 創建3個matrix，紀錄分數用
    M = np.zeros((len_1, len_2), dtype=int)
    Ix = np.zeros((len_1, len_2), dtype=int)
    Iy = np.zeros((len_1, len_2), dtype=int)
    gop_score = int(gop)
    gep_score = int(gep)

    ### 再創建3個matrix用來記錄footprint
    M_step = np.zeros((len_1, len_2), dtype=int)
    Ix_step = np.zeros((len_1, len_2), dtype=int)
    Iy_step = np.zeros((len_1, len_2), dtype=int)

    ###Start to do dynamic programming
    for i in range(len_1):
        for j in range(len_2):
            # Initial the matrix
            if i == 0 or j == 0:
                if i == 0 and j != 0:
                    if j == 1:
                        M[i, j] = M[i, j - 1] + gop_score
                        Ix[i, j] = Ix[i, j - 1] + gop_score
                        Iy[i, j] = Iy[i, j - 1] + gop_score
                    else:
                        M[i, j] = M[i, j - 1] + gep_score
                        Ix[i, j] = Ix[i, j - 1] + gep_score
                        Iy[i, j] = Iy[i, j - 1] + gep_score
                elif i != 0 and j == 0:
                    if i == 1:
                        M[i, j] = M[i - 1, j] + gop_score
                        Ix[i, j] = Ix[i - 1, j] + gop_score
                        Iy[i, j] = Iy[i - 1, j] + gop_score
                    else:
                        M[i, j] = M[i - 1, j] + gep_score
                        Ix[i, j] = Ix[i - 1, j] + gep_score
                        Iy[i, j] = Iy[i - 1, j] + gep_score

                M_step[i, j] = 999
                Ix_step[i, j] = 999
                Iy_step[i, j] = 999

            # Start to caculate three matrix's score
            else:
                ## M matrix
                pam_score = int(pam_df[seq1[i]][seq2[j]])

                M_m = M[i - 1, j - 1] + pam_score
                M_x = Ix[i - 1, j - 1] + pam_score
                M_y = Iy[i - 1, j - 1] + pam_score

                # Choose the highest score to store in M matrix
                select = [M_m, M_x, M_y]
                M[i, j] = max(select)
                # Record the footprint
                where = np.argmax(select)
                M_step[i, j] = 0 if (where == 0) else (1 if (where == 1) else 2)

                ## Ix Matrix
                Ix_m = M[i - 1, j] + gop_score
                Ix_x = Ix[i - 1, j] + gep_score
                Ix_y = Iy[i - 1, j] + gop_score

                # Choose the highest score to store in M matrix
                select = [Ix_m, Ix_x, Ix_y]
                Ix[i, j] = max(select)
                # Record the footprint
                where = np.argmax(select)
                Ix_step[i, j] = 0 if (where == 0) else (1 if (where == 1) else 2)

                ## Iy Matrix
                Iy_m = M[i, j - 1] + gop_score
                Iy_x = Ix[i, j - 1] + gop_score
                Iy_y = Iy[i, j - 1] + gep_score

                # Choose the highest score to store in M matrix
                select = [Iy_m, Iy_x, Iy_y]
                Iy[i, j] = max(select)
                # Record the foodprint
                where = np.argmax(select)
                Iy_step[i, j] = 0 if (where == 0) else (1 if (where == 1) else 2)

    ### Trace Back
    road = []

    i = len_1 - 1
    j = len_2 - 1

    first_step = np.argmax([M[i, j], Ix[i, j], Iy[i, j]])
    now_matrix = "M" if (first_step == 0) else ("Ix" if (first_step == 1) else "Iy")

    while i > 0 and j > 0:
        if now_matrix == "M":
            next_matrix = M_step[i][j]
            road.append(0)
            i -= 1
            j -= 1
        elif now_matrix == "Ix":
            next_matrix = Ix_step[i][j]
            road.append(1)
            i -= 1
        elif now_matrix == "Iy":
            next_matrix = Iy_step[i][j]
            road.append(2)
            j -= 1
        now_matrix = (
            "M" if (next_matrix == 0) else ("Ix" if (next_matrix == 1) else "Iy")
        )

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
        elif current_road == 2:
            global_seq1 += "-"
            global_seq2 += seq2[s2]
            s2 += 1
        elif current_road == 1:
            global_seq1 += seq1[s1]
            global_seq2 += "-"
            s1 += 1

    return global_seq1, global_seq2


def local_alignment(seq_compare_list, pam_df, gop, gep):
    ### Load Data
    seq1 = "-" + str(seq_compare_list[0][0])
    seq2 = "-" + str(seq_compare_list[0][1])
    len_1, len_2 = len(seq1), len(seq2)

    ### 創建3個matrix，紀錄分數用
    M = np.zeros((len_1, len_2), dtype=int)
    Ix = np.zeros((len_1, len_2), dtype=int)
    Iy = np.zeros((len_1, len_2), dtype=int)
    gop_score = int(gop)
    gep_score = int(gep)

    ### 再創建3個matrix用來記錄footprint
    M_step = np.zeros((len_1, len_2), dtype=int)
    Ix_step = np.zeros((len_1, len_2), dtype=int)
    Iy_step = np.zeros((len_1, len_2), dtype=int)

    ###Start to do dynamic programming
    for i in range(len_1):
        for j in range(len_2):
            # Initial the matrix
            if i == 0 or j == 0:
                M_step[i, j] = 999
                Ix_step[i, j] = 999
                Iy_step[i, j] = 999

            # Start to caculate three matrix's score
            else:
                ## M matrix
                pam_score = int(pam_df[seq1[i]][seq2[j]])

                M_m = M[i - 1, j - 1] + pam_score
                M_x = Ix[i - 1, j - 1] + pam_score
                M_y = Iy[i - 1, j - 1] + pam_score

                # Choose the highest score to store in M matrix
                select = [M_m, M_x, M_y]
                M[i, j] = max(select) if (max(select) > 0) else 0
                # Record the footprint
                where = np.argmax(select) if (max(select) > 0) else 3
                M_step[i, j] = (
                    0
                    if (where == 0)
                    else (1 if (where == 1) else (2 if (where == 2) else 3))
                )  # 3=Restart

                ## Ix Matrix
                Ix_m = M[i - 1, j] + gop_score
                Ix_x = Ix[i - 1, j] + gep_score
                Ix_y = Iy[i - 1, j] + gop_score

                # Choose the highest score to store in M matrix
                select = [Ix_m, Ix_x, Ix_y]
                Ix[i, j] = max(select) if (max(select) > 0) else 0
                # Record the footprint
                where = np.argmax(select) if (max(select) > 0) else 3
                Ix_step[i, j] = (
                    0
                    if (where == 0)
                    else (1 if (where == 1) else (2 if (where == 2) else 3))
                )

                ## Iy Matrix
                Iy_m = M[i, j - 1] + gop_score
                Iy_x = Ix[i, j - 1] + gop_score
                Iy_y = Iy[i, j - 1] + gep_score

                # Choose the highest score to store in M matrix
                select = [Iy_m, Iy_x, Iy_y]
                Iy[i, j] = max(select) if (max(select) > 0) else 0
                # Record the footprint
                where = np.argmax(select) if (max(select) > 0) else 3
                Iy_step[i, j] = (
                    0
                    if (where == 0)
                    else (1 if (where == 1) else (2 if (where == 2) else 3))
                )

    ### Trace Back
    roads_get = []

    ### Find the max score's index (maybe not the only one)
    where_from = np.argmax([np.max(M), np.max(Ix), np.max(Iy)])
    max_where = (
        np.where(M == np.max(M))
        if (where_from == 0)
        else (
            np.where(Ix == np.max(Ix))
            if (where_from == 1)
            else np.where(Iy == np.max(Iy))
        )
    )
    max_count = len(max_where[0])
    now_matrix = "M" if (where_from == 0) else ("Ix" if (where_from == 1) else "Iy")

    trace_index = []
    for i in range(max_count):
        trace_index.append([max_where[0][i], max_where[1][i]])
        i += 1

    for i, j in trace_index:
        road = []
        end_i, end_j = i, j

        while i > 0 or j > 0:
            # print(i, j)
            if i < 2 or j < 2:
                if (i == 0 and j > 0) or (i > 0 and j == 0):
                    break

                elif (i == 1 and j > 0) or (j == 1 and i > 0):
                    if now_matrix == "M":
                        # print(M_step[i, j])
                        if M_step[i, j] == 0:
                            road.append(0)
                            i -= 1
                            j -= 1
                        elif M_step[i, j] == 1:
                            road.append(2)
                            j -= 1
                        elif M_step[i, j] == 2:
                            road.append(1)
                            i -= 1
                        elif M_step[i, j] == 3:
                            break
                    elif now_matrix == "Ix":
                        # print(Ix_step[i, j])
                        if Ix_step[i, j] == 0:
                            road.append(0)
                            i -= 1
                            j -= 1
                        elif Ix_step[i, j] == 1:
                            road.append(2)
                            j -= 1
                        elif Ix_step[i, j] == 2:
                            road.append(1)
                            i -= 1
                        elif Ix_step[i, j] == 3:
                            break
                    elif now_matrix == "Iy":
                        # print(Iy_step[i, j])
                        if Iy_step[i, j] == 0:
                            road.append(0)
                            i -= 1
                            j -= 1
                        elif Iy_step[i, j] == 1:
                            road.append(2)
                            j -= 1
                        elif Iy_step[i, j] == 2:
                            road.append(1)
                            i -= 1
                        elif Iy_step[i, j] == 3:
                            break

            else:
                if now_matrix == "M":
                    if M_step[i, j] == 3:
                        break
                    else:
                        next_matrix = M_step[i, j]
                        road.append(0)
                        i -= 1
                        j -= 1
                elif now_matrix == "Ix":
                    if Ix_step[i, j] == 3:
                        break
                    else:
                        next_matrix = Ix_step[i, j]
                        road.append(1)
                        i -= 1
                elif now_matrix == "Iy":
                    if Iy_step[i, j] == 3:
                        break
                    else:
                        next_matrix = Iy_step[i, j]
                        road.append(2)
                        j -= 1

                now_matrix = (
                    "M"
                    if (next_matrix == 0)
                    else ("Ix" if (next_matrix == 1) else "Iy")
                )

        road.reverse()
        roads_get.append([i, j, end_i, end_j, road])
    # print(roads_get)

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
        elif current_road == 2:
            local_seq1 += "-"
            local_seq2 += seq2[s2]
            s2 += 1
        elif current_road == 1:
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


def alignment(input_path, score_path, output_path, aln, gap_open, gap_extend):
    seq_name, seq_compare_list = read_fasta(input_path)
    pam_df = read_PAM(score_path)

    if aln == "global":
        seq1, seq2 = global_alignment(seq_compare_list, pam_df, gap_open, gap_extend)
        output_fasta(seq1, seq2, seq_name, output_path)

    elif aln == "local":
        seq1, seq2 = local_alignment(seq_compare_list, pam_df, gap_open, gap_extend)
        output_fasta(seq1, seq2, seq_name, output_path)

    else:
        print("There are some errors in your parameters.")


# alignment(
#     "./examples/test.fasta",
#     "./examples/pam250.txt",
#     "result_global.fasta",
#     "global",
#     -10,
#     -2,
# )
