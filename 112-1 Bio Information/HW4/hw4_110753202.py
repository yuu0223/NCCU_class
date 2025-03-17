"""
Author: yen-nan ho
Github: https://github.com/aaron1aaron2
Date: 2022.05.11
"""
import argparse
import itertools
import functools

import os
import numpy as np

def get_args():
    """ 讀取使用者輸入參數 """
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", help = "輸入的 fasta 檔案", type = str, default = 'test.fasta')
    parser.add_argument("--output", help = "result檔案位置", type = str, default = 'result.fasta')

    parser.add_argument("--score", help = "分數表", type = str, default = 'pam250.txt')
    parser.add_argument("--gap_open", help = "gap open penalty", type = int, default = -10)
    parser.add_argument("--gap_extend", help = "gap extend penalty", type = int, default = -2)

    parser.add_argument("--aln", help = "global、local", type = str, default = 'global')

    parser.add_argument("--record_path", help = "記錄計算過程檔案", type = str, default = None)
    parser.add_argument("--output_record", help = "是否要輸出 record", action='store_true')

    args = parser.parse_args()

    return args


def read_fasta(path):
    """ 讀取 fasta 資料 """
    with open(path, 'r') as f:
        read = f.readlines()
        # seq_read = {line1.strip()[1:]:line2.strip() for line1,line2 in itertools.zip_longest(*[f]*2)}

    seq1, seq2, seq1_id, seq2_id = '', '', '', ''
    flag = True
    for i in read:
        if (seq1 != '') & (i.find('>') != -1): 
            flag = False

        if i.find('>') != -1:
            if flag:
                seq1_id = i.strip() # remove \n
            else:
                seq2_id = i.strip() # remove \n
        else:
            if flag:
                seq1 += i.strip() # remove \n
            else:
                seq2 += i.strip()

    seq_read = {
        seq1_id: seq1,
        seq2_id: seq2
    }

    return seq_read


def read_PAM(path):
    """ 讀取 PAM 資料(自動排除開頭 # 註解的 row) """
    with open(path, 'r') as f:
        input_matrix = [[cell for cell in line.split()] for line in f.readlines() if line[0]!='#']

    amino_acids = input_matrix[0]
    score_dt = {(row[0], amino_acids[idx]):i for row in input_matrix[1:] for idx,i in enumerate(row[1:])}

    return score_dt


def global_alignment(seq_pairs:tuple, score_dt:dict, gap_open:int, gap_extend:int, record_path:str, output_record:bool):
    """ global alignment 並輸出結果 

    # 雙 gap 問題，註解中代號的相對位置:

        a  b
        c
    """ 
    seq1, seq2 = seq_pairs
    L1, L2 = len(seq1), len(seq2)

    # step 1: 計算矩陣中每格分數與 alignment 的方式 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # 取分數的邏輯，規定第一欄 or 排的分數，第一個 gap_open 後面都是 gap_extend| 註: i、j 最小只會到 -1
    get_score = lambda i, j, mat: mat[i,j] if ((i>=0) & (j>=0)) else (
        0 if ((i<0) & (j<0)) else (
            gap_open + (gap_extend*j) if i<0 else gap_open + (gap_extend*i)
        )
    )

    # 紀錄分數和 alignment 選擇的矩陣
    score_mat = np.zeros((L1, L2), dtype=int)
    select_type_mat = np.zeros((L1, L2), dtype=int)

    # 紀錄 s1 & s2 gap 同分的情況。同分時，比較與 rollback | 邏輯 -> 程式從上而下，先遇到右邊的在遇到下面的
    select_down_ls =  [] # 在這裡面存的是有更正路線的位置，也就是確定向上開 gap 
    same_score_dt = {} # 儲存個位置與向右的分數

    i = 0
    while i < L1:
        j = 0
        while j < L2:
            # 三種類型 s1_gap(0)、s2_gap(1)、no_gap(2) -------------------
            no_gap = get_score(i-1, j-1, score_mat) + int(score_dt[seq1[i], seq2[j]])
            
            # 看左邊
            s1_gap_s = get_score(i, j-1, score_mat)
            s1_gap = s1_gap_s + gap_extend if ((j>0) & (select_type_mat[i, j-1] == 0)) else (
                s1_gap_s + gap_open
            )

            # 看上面
            s2_gap_s = get_score(i-1, j, score_mat)
            s2_gap = (s2_gap_s + gap_extend) if ((i>0) & (select_type_mat[i-1, j] == 1)) else (
                s2_gap_s + gap_open
            )

            # s1 & s2 gap 同分的情況(處理) -------------------------------
            # (4) a 當確定位置已經要往下 (在 select_down_ls 中)，要控制不往左開 gap(預設)，讓程式往上開 gap。
            if ((i,j) in select_down_ls): s1_gap = s2_gap - 1 # 讓 s1_gap 比較小，後續不會被選擇到

            # 同分時以 list 的前後順序優先，預設順序 -----------------------
            ali_select_ls = [s1_gap, s2_gap, no_gap]
            select_type = np.argmax(ali_select_ls)

            # 儲存分數 & alignment 選擇 ----------------------------------
            score_mat[i, j] = ali_select_ls[select_type]
            select_type_mat[i, j] = select_type

            # s1 & s2 gap 同分的情況(偵測與比較) -------------------------
            # (1) a 確認是不是雙 gap 並將位置加入 same_score_dt -> (s1 gap 與 s2 gap 的分數相同) & (i、j 都要 > 0) & (會選擇gap時)
            if (s1_gap == s2_gap) & (i>0) & (j>0) & (select_type in [0, 1]): 
                same_score_dt.update({(i,j):None}) # 第一次面臨雙 gap 選擇的位置

            # (2) b 往左看，且符合條件，將自己當下的分數更到 same_score_dt -> (如果是雙 gap 的位置) & (有向左 extend 發生)
            if ((i, j-1) in same_score_dt) & (select_type == 0): 
                same_score_dt.update({(i,j-1):score_mat[i, j]}) # 紀錄向右分數
            
            # (3) c 看上面，如有雙 gap 現象，確認記錄任分數，並 rollback 到選擇雙 gap 的位置
            if ((i-1, j) in same_score_dt) & ((i-1, j) not in select_down_ls):
                s2_gap_s_ex = s2_gap_s + gap_extend
                rollback = False

                if same_score_dt[i-1, j] == None:
                    rollback = True # 這邊代表 b 位置沒有選擇對 a extend
                else:
                    if s2_gap_s_ex > same_score_dt[i-1, j]:
                        rollback = True # 這邊代表 c 位置 extend 分數 > b 位置 extend 分數

                if rollback:
                    select_down_ls.append((i-1, j))
                    i -= 1 ; j -= 1
            j += 1
        # 到下一行
        i += 1 
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    # step 2: 由座右下角開始 trace back >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    track_back_mat = np.full((L1, L2), '-')
    route = []
    i, j = L1-1, L2-1 # 走路用位置

    while (i>=0) | (j>=0):
        cur_type = select_type_mat[i, j]
        route.append(cur_type)
        track_back_mat[i, j] = cur_type
        # 控制下一步要走的地方
        if cur_type == 2:
            i -= 1; j -= 1
        elif cur_type == 0:
            j -= 1
        else:
            i -= 1

    route.reverse() # 翻轉 -> 從左上開始
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    # step 3: 整理排序好的序列 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ali_seq1, ali_seq2 = '', '' # 排序結果
    pt1, pt2 = 0, 0 # seq、seq2 指標
    for cur_type in route:
        if cur_type == 2:
            # 配對
            ali_seq1 += seq1[pt1] ; ali_seq2 += seq2[pt2]
            pt1 += 1 ; pt2 += 1
        elif cur_type == 0:
            # s1 gap
            ali_seq1 += '-' ; ali_seq2 += seq2[pt2]
            pt2 += 1
        else:
            # s2 gap
            ali_seq1 += seq1[pt1] ; ali_seq2 += '-'
            pt1 += 1
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    # 檢查用非必要
    if output_record:
        score_dt_sort = {}
        for a1, a2 in score_dt.keys():
            if ((a1, a2) not in score_dt_sort) & ((a2, a1) not in score_dt_sort):
                score_dt_sort.update({(a1, a2):score_dt[a1, a2]})
        score_dt_sort = [['-'.join(k), v] for k, v in score_dt_sort.items() if '*' not in k]
        open(record_path, 'a').write(f'[Pam score]\n{score_dt_sort}\n\n')
        output_matrix(seq1, seq2, score_mat, '[Score matrix]', record_path)
        output_matrix(seq1, seq2, select_type_mat, '[Select type matrix]', record_path)
        output_matrix(seq1, seq2, track_back_mat, '[Track back matrix]', record_path)
        open(record_path, 'a').write(f'[Route 1]\n{route}\n\n')
        open(record_path, 'a').write(f'[Route 1]\n {ali_seq1}\n {ali_seq2}\n\n' + '='*50 + '\n')
        open(record_path, 'a').write(f'[same score detect] None: mean no competition\n {same_score_dt}\n\n')
        open(record_path, 'a').write(f'[select down ls]\n {select_down_ls}\n\n')


    return ali_seq1, ali_seq2

def local_alignment(seq_pairs:tuple, score_dt:dict, gap_open:int, gap_extend:int, record_path:str, output_record:bool):
    """ local alignment 並輸出結果 """ 
    seq1, seq2 = seq_pairs
    L1, L2 = len(seq1), len(seq2)

    # step 1: 計算矩陣中每格分數與 alignment 的方式 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # 規定第一欄 or 排的分數，正常是取對應位置的值
    get_score = lambda i, j, mat: mat[i,j] if ((i>=0) & (j>=0)) else 0

    # 紀錄分數和 alignment 選擇的矩陣
    score_mat = np.zeros((L1, L2), dtype=int)
    select_type_mat = np.zeros((L1, L2), dtype=int)

    # 紀錄 s1 & s2 gap 同分的情況。同分時，比較與 rollback | 邏輯 -> 程式從上而下，先遇到右邊的在遇到下面的
    select_down_ls =  [] # 在這裡面存的是有更正路線的位置，也就是確定向上開 gap 
    same_score_dt = {} # 儲存個位置與向右的分數

    i = 0
    while i < L1:
        j = 0
        while j < L2:
            # 四種類型 s1_gap(0)、s2_gap(1)、no_gap(2)、restart(3) --------------
            no_gap = get_score(i-1, j-1, score_mat) + int(score_dt[seq1[i], seq2[j]])

            # 看左邊
            s1_gap_s = get_score(i, j-1, score_mat)
            s1_gap = s1_gap_s + gap_extend if ((j>0) & (select_type_mat[i, j-1] == 0)) else (
                s1_gap_s + gap_open
            )

            # 看上面
            s2_gap_s = get_score(i-1, j, score_mat)
            s2_gap = (s2_gap_s + gap_extend) if ((i>0) & (select_type_mat[i-1, j] == 1)) else (
                s2_gap_s + gap_open
            )

            # 重新開始
            restart = 0

            # s1 & s2 gap 同分的情況(處理) --------------------------------------
            # (4) a 當確定位置已經要往下 (在 select_down_ls 中)，要控制不往左開 gap(預設)，讓程式往上開 gap。
            if ((i,j) in select_down_ls): s1_gap = s2_gap - 1 # 讓 s1_gap 比較小，後續不會被選擇到


            # 同分時以前面的順序優先(助教說選哪個沒關係) --------------------------
            ali_select_ls = [s1_gap, s2_gap, no_gap, restart]
            select_type = np.argmax(ali_select_ls)

            # 儲存分數 & alignment 選擇 ----------------------------------------
            select_type_mat[i, j] = select_type
            score_mat[i, j] = ali_select_ls[select_type]

            # s1 & s2 gap 同分的情況(偵測與比較) -------------------------------
            # (1) a 確認是不是雙 gap 並將位置加入 same_score_dt -> (s1 gap 與 s2 gap 的分數相同) & (i、j 都要 > 0) & (會選擇gap時)
            if (s1_gap == s2_gap) & (i>0) & (j>0) & (select_type in [0, 1]): 
                same_score_dt.update({(i,j):None}) # 第一次面臨雙 gap 選擇的位置

            # (2) b 往左看，且符合條件，將自己當下的分數更到 same_score_dt -> (如果是雙 gap 的位置) & (有向左 extend 發生)
            if ((i, j-1) in same_score_dt) & (select_type == 0): 
                same_score_dt.update({(i,j-1):score_mat[i, j]}) # 紀錄向右分數
            
            # (3) c 看上面，如有雙 gap 現象，確認記錄任分數，並 rollback 到選擇雙 gap 的位置
            if ((i-1, j) in same_score_dt) & ((i-1, j) not in select_down_ls):
                s2_gap_s_ex = s2_gap_s + gap_extend
                rollback = False

                if same_score_dt[i-1, j] == None:
                    rollback = True # 這邊代表 b 位置沒有選擇對 a extend
                else:
                    if s2_gap_s_ex > same_score_dt[i-1, j]:
                        rollback = True # 這邊代表 c 位置 extend 分數 > b 位置 extend 分數

                if rollback:
                    select_down_ls.append((i-1, j))
                    i -= 1 ; j -= 1
            j += 1

        # 到下一行
        i += 1

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    # step 2: trace back >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # 決定 traceback 起始點(分數 > 長度)

    max_i, max_j = np.where(score_mat == score_mat.max()) # score_mat.argmax(axis=0) # 對 col -> j 他只會返回一個
    total_max = list(zip(max_i, max_j))

    track_back_mat = np.full((L1, L2), '-')
    routes = {}
    for i in total_max:
        route = []
        i, j = tail_i , tail_j =  i

    # 開始 traceback
        while (i>=0) | (j>=0):
            cur_type = select_type_mat[i, j]
            route.append(cur_type)
            track_back_mat[i, j] = cur_type
            # 控制下一步要走的地方
            if (i==0) | (j==0):
                break
            elif cur_type == 2:
                i -= 1; j -= 1
            elif cur_type == 0:
                j -= 1
            elif cur_type == 1:
                i -= 1
            else:
                break

        route.reverse() # 翻轉 -> 從左上開始
        routes.update({(i, j, tail_i , tail_j) : route}) # key: (起始i, 起始j, 結束i , 結束j) | values: route
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    # step 3: 整理排序好的序列 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    max_length = 0
    result_pair = []
    for (head_i, head_j, tail_i, tail_j), route in routes.items():
        ali_seq1, ali_seq2 = '', '' # 排序結果
        pt1, pt2 = head_i, head_j # seq、seq2 指標
        for cur_type in route:
            if cur_type in [2, 3]:
                ali_seq1 += seq1[pt1] ; ali_seq2 += seq2[pt2]
                pt1 += 1 ; pt2 += 1
            elif cur_type == 0:
                ali_seq1 += '-' ; ali_seq2 += seq2[pt2]
                pt2 += 1
            else:
                ali_seq1 += seq1[pt1] ; ali_seq2 += '-'
                pt1 += 1

        if len(ali_seq1) > max_length:
            max_length = len(ali_seq1)

        result_pair.append((ali_seq1, ali_seq2))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    if output_record:
        score_dt_sort = {}
        for a1, a2 in score_dt.keys():
            if ((a1, a2) not in score_dt_sort) & ((a2, a1) not in score_dt_sort):
                score_dt_sort.update({(a1, a2):score_dt[a1, a2]})
        score_dt_sort = [['-'.join(k), v] for k, v in score_dt_sort.items() if '*' not in k]
        open(record_path, 'a').write(f'[Pam score]\n{score_dt_sort}\n\n')
        output_matrix(seq1, seq2, score_mat, '[Score matrix]', record_path)
        output_matrix(seq1, seq2, select_type_mat, '[Select type matrix]', record_path)
        output_matrix(seq1, seq2, track_back_mat, '[Track back matrix]', record_path)
        open(record_path, 'a').write(''.join([f'[Route {idx+1}] ({i1},{j1}) -> {i2, j2}) |score {score_mat[i2, j2]} \n{r}\n\n'
                        for idx, ((i1, j1, i2, j2), r) in enumerate(routes.items())]))
        open(record_path, 'a').write(''.join([f'[Route {idx+1}]\n {s1}\n {s2}\n\n'
                        for idx, (s1, s2) in enumerate(result_pair)]))
        open(record_path, 'a').write('[Result]\n' + '、'.join([f'Route {idx+1}'
                        for idx, (s1, _) in enumerate(result_pair) if len(s1) == max_length])  + '='*50 + '\n')
        open(record_path, 'a').write(f'[same score detect] None: mean no competition\n {same_score_dt}\n\n')
        open(record_path, 'a').write(f'[select down ls]\n {select_down_ls}\n\n')

    # 分數相同的以長度最長的為準，只要是最長的就輸出
    result_pair = [i for i in result_pair if len(i[0]) == max_length]

    return result_pair

def output_matrix(seq1:list, seq2:list, matrix:np.array, msg:str, path:str):
    assert matrix.shape == (len(seq1), len(seq2)), \
        f"[print_matrix] The dimensions of input matrix must be {(len(seq1), len(seq2))}"
    
    import pandas as pd
    df = pd.DataFrame(matrix, columns=[i for i in seq2], index=[i for i in seq1])
    with open(path, 'a') as f:
        f.write(msg + '\n' + df.to_string() + '\n' + '='*50 + '\n\n')


def output_fasta(seqID:str, seqdata:str, path:str):
    with open(path, 'a') as f:
        f.write(f'{seqID}\n')
        f.write(f'{seqdata}\n')

def main():
    args = get_args()

    if args.record_path == None:
        args.record_path = args.output.replace('.fasta', '_record.txt')

    folder, _ = os.path.split(args.output)
    if folder != '':
        os.makedirs(folder, exist_ok=True)

    if args.output_record:
        open(args.record_path, 'a').write('='*50 + '\n[args]\n' + str(args)[10 : -1] + '\n' + '='*50 + '\n')

    # Step 1: read data =============================
    seq_read_dt = read_fasta(args.input) # id & sequence data  | format: {'sequence id': 'sequence data' ...}
    score_dt = read_PAM(args.score) # Save the score that matches letter1 to letter2 | format: {(letter1, letter2): score ...} 

    # Step 2: Get pairs =============================
    pairs_id_ls = list(itertools.combinations(seq_read_dt.keys(), 2))
    pairs_ls = [(seq_read_dt[id1], seq_read_dt[id2]) 
                    for id1, id2 in pairs_id_ls]


    # Step 3: alignment ==========
    if args.aln == 'global':

        result_ls = list(map(
            functools.partial(
                global_alignment, 
                score_dt = score_dt,
                gap_open = args.gap_open,
                gap_extend = args.gap_extend,
                record_path = args.record_path,
                output_record = args.output_record
            ), pairs_ls))

        # 輸出資料
        for task in range(len(pairs_id_ls)):
            for seqID, seq in list(zip(pairs_id_ls[task], result_ls[task])):
                output_fasta(seqID, seq, args.output)

    elif args.aln == 'local':

        result_ls = list(map(
            functools.partial(
                local_alignment, 
                score_dt = score_dt,
                gap_open = args.gap_open,
                gap_extend = args.gap_extend,
                record_path = args.record_path,
                output_record = args.output_record
            ), pairs_ls))

        # 輸出資料
        for task in range(len(pairs_id_ls)):
            for muti_result in result_ls[task]:
                for seqID, seq in list(zip(pairs_id_ls[task], muti_result)):
                    output_fasta(seqID, seq, args.output)
    else:
        raise ValueError("'--aln' must be 'global' or 'local'")



if __name__ == '__main__':
    main()


"""
# 測試範例
python hw4_110753202.py --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result_global(my).fasta


# global
python hw4_110753202.py --input test.fasta --output output/gb_pm100_o-gap10_e-gap2.fasta --aln global --score pam100.txt --gap_open -10 --gap_extend -2 --output_record
python hw4_110753202.py --input test.fasta --output output/gb_pm250_o-gap10_e-gap2.fasta --aln global --score pam250.txt --gap_open -10 --gap_extend -2 --output_record

python hw4_110753202.py --input test.fasta --output output/gb_pm100_o-gap5_e-gap2.fasta --aln global --score pam100.txt --gap_open -5 --gap_extend -2 --output_record
python hw4_110753202.py --input test.fasta --output output/gb_pm250_o-gap5_e-gap2.fasta --aln global --score pam250.txt --gap_open -5 --gap_extend -2 --output_record

# local
python hw4_110753202.py --input test.fasta --output output/lc_pm100_o-gap10_e-gap2.fasta --aln local --score pam100.txt --gap_open -10 --gap_extend -2 --output_record
python hw4_110753202.py --input test.fasta --output output/lc_pm250_o-gap10_e-gap2.fasta --aln local --score pam250.txt --gap_open -10 --gap_extend -2 --output_record

python hw4_110753202.py --input test.fasta --output output/lc_pm100_o-gap5_e-gap2.fasta --aln local --score pam100.txt --gap_open -5 --gap_extend -2 --output_record
python hw4_110753202.py --input test.fasta --output output/lc_pm250_o-gap5_e-gap2.fasta --aln local --score pam250.txt --gap_open -5 --gap_extend -2 --output_record

# test4 (無答案)
python hw4_110753202.py --input test_data/test4.fasta --output output/gb_pm250_o-gap10_e-gap2.test4.fasta --aln global --score pam250.txt --gap_open -10 --gap_extend -2 --output_record
python hw4_110753202.py --input test_data/test4.fasta --output output/lc_pm250_o-gap10_e-gap2.test4.fasta --aln local --score pam250.txt --gap_open -10 --gap_extend -2 --output_record

# test1 
python hw4_110753202.py --input test_data/test1.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output output/gb_pm250_o-gap10_e-gap2.test1.fasta --output_record
"""