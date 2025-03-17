import numpy as np


def read_inputs(input_path, score_path):
    # loading input
    file = open(input_path)
    input = [[], []]
    protein = []
    n = -1  # line count
    for line in file.readlines():
        data = line
        data = [i for i in data if (i and i != "\n")]
        if (
            line[0] != ">"
        ):  # record only even lines (since odd lines are names of protein)
            for x in data:
                input[n].append(x)
        else:
            protein.append(line)
            n += 1  # line count
    file = open(score_path)
    pam = []  # creating Pam
    for line in file.readlines():  # proccessing pam txt
        data = line.strip().split(" ")  # splitting by " "
        data = [i for i in data if i]  # adds i (for each i in data) if i exists
        if (
            data[0] != "#"
        ):  # removing lines that starts with # in the beginning of document
            pam.append(data[1:])  # removing 1st letter and appending
    pam = pam[1:]  # removing top row
    for i in range(len(pam)):
        for j in range(len(pam[i])):
            pam[i][j] = int(pam[i][j])
    return (input, protein, pam)


def create_file(output_path, output, protein_name):
    with open(output_path, "w") as file:
        for i in range(len(output)):
            line = output[i]
            file.write(protein_name[i % 2])
            for j in range(len(line)):
                file.write(str(line[j]))
            file.write("\n")


def find_max_psa(psa, dim):
    max_score = 0
    coord = []
    for i in range(dim[0]):
        for j in range(dim[1]):
            for k in range(3):
                if psa[k][i][j] > max_score:
                    coord = []
                    coord.append([i, j, k])
                    max_score = psa[k][i][j]
                elif psa[k][i][j] == max_score:
                    coord.append([i, j, k])
    return coord


def gen_initial(input, dim):
    psa = np.zeros((3, dim[0], dim[1]))
    path = np.zeros((3, dim[0], dim[1]))
    return psa, path


def psa_gen(input, pam, mapping, gap, gep, aln):
    dim = [len(input[0]) + 1, len(input[1]) + 1]
    psa, path = gen_initial(input, dim)
    for i in range(dim[0]):
        for j in range(dim[1]):
            psa[0][0][0] = 0
            psa[1][0][0] = 0
            psa[2][0][0] = 0
            if i != 0 and j == 0:
                psa[0][i][j] = psa[1][i - 1][j] + gap * int(i == 1) + gep * int(i != 1)
                psa[1][i][j] = psa[1][i - 1][j] + gap * int(i == 1) + gep * int(i != 1)
                psa[2][i][j] = psa[1][i - 1][j] + gap * int(i == 1) + gep * int(i != 1)
            elif i == 0 and j != 0:
                psa[0][i][j] = psa[2][i][j - 1] + gap * int(j == 1) + gep * int(j != 1)
                psa[1][i][j] = psa[2][i][j - 1] + gap * int(j == 1) + gep * int(j != 1)
                psa[2][i][j] = psa[2][i][j - 1] + gap * int(j == 1) + gep * int(j != 1)
            else:  # the rest
                pam_score = pam[mapping[input[0][i - 1]]][mapping[input[1][j - 1]]]
                psa[0][i][j] = max(
                    psa[0][i - 1][j - 1] + pam_score,
                    psa[1][i - 1][j - 1] + pam_score,
                    psa[2][i - 1][j - 1] + pam_score,
                )
                if psa[0][i][j] == psa[0][i - 1][j - 1] + pam_score:
                    path[0][i][j] = 0
                elif psa[0][i][j] == psa[1][i - 1][j - 1] + pam_score:
                    path[0][i][j] = 1
                else:
                    path[0][i][j] = 2
                psa[1][i][j] = max(
                    psa[0][i - 1][j] + gap,
                    psa[1][i - 1][j] + gep,
                    psa[2][i - 1][j] + gap,
                )
                if psa[1][i][j] == psa[0][i - 1][j] + gap:
                    path[1][i][j] = 0
                elif psa[1][i][j] == psa[1][i - 1][j] + gep:
                    path[1][i][j] = 1
                else:
                    path[1][i][j] = 2
                psa[2][i][j] = max(
                    psa[0][i][j - 1] + gap,
                    psa[1][i][j - 1] + gap,
                    psa[2][i][j - 1] + gep,
                )  # M_y extend)
                if psa[2][i][j] == psa[0][i][j - 1] + gap:
                    path[2][i][j] = 0
                elif psa[2][i][j] == psa[1][i][j - 1] + gap:
                    path[2][i][j] = 1
                else:
                    path[2][i][j] = 2
            if aln == "local":
                for n in range(3):
                    if psa[n][i][j] <= 0:
                        psa[n][i][j] = 0
                        path[n][i][j] = -1
    print("M: ", psa[0], "\n X: ", psa[1], "\n Y: ", psa[2])
    return psa, path


def output_gen(input, psa, path, aln):
    dim = [len(input[0]), len(input[1])]
    if aln == "global":
        x = dim[0]
        y = dim[1]
        max_psa = max(psa[0][x][y], psa[1][x][y], psa[2][x][y])
        if max_psa == psa[0][x][y]:
            k = 0
        elif max_psa == psa[1][x][y]:
            k = 1
        else:
            k = 2
        coord = [[x, y, k]]
    elif aln == "local":
        coord = find_max_psa(psa, dim)
    else:
        return "invalid alignment"
    output_fin = [[], []]
    for i in range(len(coord)):
        output = [[], []]
        x = coord[i][0]
        y = coord[i][1]
        k = coord[i][2]
        while x > 0 and y > 0:
            print(k, x, y)
            curr = psa[k][x][y]
            if k == 0 and (curr > 0 or aln == "global"):
                output[0].append(input[0][x - 1])
                output[1].append(input[1][y - 1])
                k = int(path[0][x][y])
                x -= 1
                y -= 1
                continue
            elif k == 1 and (curr > 0 or aln == "global"):  # y gap
                output[0].append(input[0][x - 1])
                output[1].append("-")
                k = int(path[1][x][y])
                x -= 1
                continue
            elif k == 2 and (curr > 0 or aln == "global"):  # x gap
                output[0].append("-")
                output[1].append(input[1][y - 1])
                k = int(path[2][x][y])
                y -= 1
                continue
            else:
                break
        if len(output[0]) >= len(output_fin[0]):
            output_fin[0] = output[0][::-1]
            output_fin[1] = output[1][::-1]
    return output_fin


def alignment(input_path, score_path, output_path, aln, gap_open, gap_extend):
    [input, protein_name, pam] = read_inputs(input_path, score_path)
    mapping = {
        "A": 0,
        "R": 1,
        "N": 2,
        "D": 3,
        "C": 4,
        "Q": 5,
        "E": 6,
        "G": 7,
        "H": 8,
        "I": 9,
        "L": 10,
        "K": 11,
        "M": 12,
        "F": 13,
        "P": 14,
        "S": 15,
        "T": 16,
        "W": 17,
        "Y": 18,
        "V": 19,
        "B": 20,
        "Z": 21,
        "X": 22,
        "*": 23,
    }  # X2 = *
    [psa, path] = psa_gen(input, pam, mapping, gap_open, gap_extend, aln)
    output = output_gen(input, psa, path, aln)
    create_file(output_path, output, protein_name)


def checker():
    file = open("./examples/result_global.fasta")
    input = []
    protein = []
    count = 0  # line count
    for line in file.readlines():
        data = line
        data = [i for i in data if i]
        if (
            count % 2 == 1
        ):  # record only even lines (since odd lines are names of protein)
            input.append(data[0:-1])
        else:
            protein.append(line)
        count += 1  # line count
    create_file("./result.txt", input, protein)


# alignment("./examples/test.fasta", "./examples/pam250.txt", "./output.txt", "global", -10, -2)
# checker()

alignment(
    "./examples/test2.fasta",
    "./examples/pam250.txt",
    "result_global.fasta",
    "global",
    -10,
    -2,
)
