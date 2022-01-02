from Bio.Seq import Seq
import math
from Bio import SeqIO

def codon_list():
    return ["ATT", "ATC", "ATA", "CTT", "CTC", "CTA", "CTG", "TTA", "TTG", "GTT", "GTC", "GTA", "GTG", "TTT", "TTC", "ATG",
                  "TGT", "TGC", "GCT", "GCC", "GCA", "GCG", "GGT", "GGC", "GGA", "GGG", "CCT", "CCC", "CCA", "CCG", "ACT", "ACC",
                  "ACA", "ACG", "TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "TAT", "TAC", "TGG", "CAA", "CAG", "AAT", "AAC", "CAT",
                  "CAC", "GAA", "GAG", "GAT", "GAC", "AAA", "AAG", "CGT", "CGC", "CGA", "CGG", "AGA", "AGG ,TAA", "TAG", "TGA"]

def read_file(filename):
    for seq_record in SeqIO.parse("dataFasta/" + filename + ".fasta", "fasta"):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))
        file_seq = seq_record
    return file_seq


def find_pairs(sequence):
    n = 3
    output = [(sequence.seq[i:i + n]) for i in range(0, len(sequence.seq), n)]
    # print(range(0, len(sequence.seq), n))
    length = len(output)
    # print(length)
    i = 0
    # Suskirstom i start ir stop kodonu segmentus
    seq_pairs = []
    while i < length:
        if output[i] == "ATG":
            tmp_seq = ""
            j = i
            while j < length:
                tmp_seq += output[j]
                if output[j] == "TAG" or output[j] == "TAA" or output[j] == "TGA":
                    seq_pairs.append("".join(str(tmp_seq)))
                    i = j
                    break
                j += 1
        i += 1
    return seq_pairs


def find_reverse_pairs(sequence):
    n = 3
    output = [(sequence.seq.reverse_complement()[i:i + n]) for i in range(0, len(sequence.seq), n)]
    length = len(output)
    i = 0
    # Suskirstom i start ir stop kodonu segmentus
    seq_pairs = []
    while i < length:
        if output[i] == "ATG":
            tmp_seq = ""
            j = i
            while j < length:
                tmp_seq += output[j]
                if output[j] == "TAG" or output[j] == "TAA" or output[j] == "TGA":
                    # if len(str(tmp_seq))>=100: #atrenkam sekas kurios yra daugiau nei 100 elementu
                    seq_pairs.append("".join(str(tmp_seq)))
                    i = j
                    break
                j += 1
        i += 1
    return seq_pairs


def longest_sequence(sequence):
    max_length_tag = 0
    max_length_taa = 0
    max_length_tga = 0
    tag_seq = ""
    taa_seq = ""
    tga_seq = ""
    for i in sequence:
        if len(i) > max_length_tag and str(i).endswith("TAG"):
            max_length_tag = len(i)
            tag_seq = i
        if len(i) > max_length_taa and str(i).endswith("TAA"):
            max_length_taa = len(i)
            taa_seq = i
        if len(i) > max_length_tga and str(i).endswith("TGA"):
            max_length_tga = len(i)
            tga_seq = i
    longest_final = [tag_seq,taa_seq,tga_seq]
    return longest_final


def longer_than(sequence_arr, length):
    new_seq_arr = []
    for sequence in sequence_arr:
        if len(sequence) >= length:
            new_seq_arr.append("".join(sequence))
    return new_seq_arr


def codon_frequency(sequence):
    all_codons = codon_list()
    main_sequence = "".join(sequence)
    individual_seq = [main_sequence[k:k+3] for k in range(0, len(main_sequence), 3)]
    counted_freq = []
    just_freq = []
    for j in all_codons:
        counter = 0
        for i in individual_seq:
            if i == j:
                counter+=1
            ratio = counter/len(individual_seq)
        str_ratio = "{:.5f}".format(ratio)
        counted_freq.append(j+" "+str_ratio)
        just_freq.append(str_ratio)
    return counted_freq, just_freq


def dicodon_frequency(sequence):
    all_codons = codon_list()
    all_dicodons = []
    test_counter = 0
    for i in all_codons:
        for j in all_codons:
            all_dicodons.append(i+j)
            test_counter+=1
    # print(test_counter)
    # print(all_dicodons)
    main_sequence = "".join(sequence)
    individual_seq = [main_sequence[k:k+6] for k in range(0, len(main_sequence), 3)]
    # print(individual_seq)
    # print(len(individual_seq))
    counted_freq = []
    just_freq = []
    for j in all_dicodons:
        counter = 0
        for i in individual_seq:
            if i == j:
                counter+=1
            ratio = counter/len(individual_seq)
        if(ratio > 0):
            str_ratio = "{:.5f}".format(ratio)
            counted_freq.append(j+" "+str_ratio)
            just_freq.append(str_ratio)
    # print(len(counted_freq))
    return counted_freq, just_freq


def distance_calculation(freq1, freq2):
    distance = 0.0
    # Atstumo tarp tasku funkcija  √((x2 - x1 )^2 + (y2 - y1 )^2)...
    distance += math.pow(float(freq1) - float(freq2), 2)
    root = math.sqrt(distance)
    str_root = "{:.5f}".format(root)
    return str_root


def matrix_generator(sequence):
    matrix = [[0.0 for i in range(8)] for j in range(8)]
    # print(matrix)
    for x in range(0,7):
        for y in range(x+1, 8):
            matrix[x][y] = distance_calculation(sequence[x],sequence[y])
            matrix[y][x] = matrix[x][y]
    return matrix


def result_function(filename):
    # PIRMA DALIS
    print("1.1 Kodonu star ir stop sekos (tiesiogine seka)")
    start_stop_seq = find_pairs(read_file(filename))
    print(start_stop_seq)
    print("1.2 Kodonu star ir stop sekos (reverse komplementas)")
    reverse_start_stop_seq = find_reverse_pairs(read_file(filename))
    print(reverse_start_stop_seq)
    # ANTRA DALIS
    print("2.1 Kodonu stop kodono tolimiausias start (tiesiogine seka)")
    print(longest_sequence(start_stop_seq))
    print("2.2 Kodonu stop kodono tolimiausias start (reverse)")
    print(longest_sequence(reverse_start_stop_seq))
    # TRECIA DALIS
    print("3.1 Daugiau arba lygu 100 sekos (tiesiogine)")
    greater_start_stop_seq = longer_than(start_stop_seq,100)
    print(greater_start_stop_seq)
    print("3.2 Daugiau arba lygu 100 sekos (reverse)")
    greater_reverse_start_stop_seq = longer_than(reverse_start_stop_seq,100)
    print(greater_reverse_start_stop_seq)
    # KETVIRTA DALIS
    print("4.1 Kodonu daznis (tiesiogine)")
    with_name1, no_noname1 = codon_frequency(greater_start_stop_seq)
    print(with_name1)
    print("4.2 Dikodonu daznis (tiesiogine)")
    with_name2, no_noname2 = dicodon_frequency(greater_start_stop_seq)
    print(with_name2)
    print("4.3 Kodonu daznis (reverse)")
    with_name3, no_noname3 = codon_frequency(greater_reverse_start_stop_seq)
    print(with_name3)
    print("4.3 Dikodonu daznis (reverse)")
    with_name4, no_noname4 = dicodon_frequency(greater_reverse_start_stop_seq)
    print(with_name4)
    final_codon_freq =no_noname1 + no_noname3
    final_dicodon_freq = no_noname2 + no_noname4

    return final_codon_freq, final_dicodon_freq


# main
if __name__ == '__main__':
    ba1_cod, ba1_di = result_function("bacterial1")
    ba2_cod, ba2_di = result_function("bacterial2")
    ba3_cod, ba3_di = result_function("bacterial3")
    ba4_cod, ba4_di = result_function("bacterial4")
    ma1_cod, ma1_di = result_function("mamalian1")
    ma2_cod, ma2_di = result_function("mamalian2")
    ma3_cod, ma3_di = result_function("mamalian3")
    ma4_cod, ma4_di = result_function("mamalian4")
    cod_full_list = ba1_cod + ba2_cod + ba3_cod + ba4_cod + ma1_cod + ma2_cod + ma3_cod + ma4_cod
    di_full_list = ba1_di + ba2_di + ba3_di + ba4_di + ma1_di + ma2_di + ma3_di + ma4_di
    print("FULL LISTAI")
    print(cod_full_list)
    print(len(cod_full_list))
    print(cod_full_list[1])
    print(float(cod_full_list[1]))
    codon_matrix = matrix_generator(cod_full_list)
    dicodon_matrix = matrix_generator(di_full_list)
    print("Kodonu atstumu matrica")
    for i in codon_matrix:
        print(i)
    print("Dikodonu atstumu matrica")
    for i in dicodon_matrix:
        print(i)