import os
import numpy as np
import pickle


file_dir = os.path.dirname(__file__)


def run_iupred(pairs):
    for pair in pairs:
        os.system("pwd")
        pep_id = pair[1]
        print(pep_id)
        os.system(f"python iupred/iupred2a.py peptide_fastas/{pep_id}.fasta long > iupred/iupred_results/{pep_id}_long.txt")
        os.system(f"python iupred/iupred2a.py -a peptide_fastas/{pep_id}.fasta short > iupred/iupred_results/{pep_id}_short.txt")

        prot_id = pair[0]
        print(prot_id)
        os.system(f"python iupred/iupred2a.py protein_fastas/{prot_id}.fasta long > iupred/iupred_results/{prot_id}_long.txt")
        os.system(f"python iupred/iupred2a.py -a protein_fastas/{prot_id}.fasta short > iupred/iupred_results/{prot_id}_short.txt")

    f = open("full_pep_long.result", "w")
    g = open("full_pep_short.result", "w")

    pep_ids = list(zip(*pairs))[1]

    for pep_id in pep_ids:
        with open(f"iupred/iupred_results/{pep_id}_long.txt", "r") as h:
            f.writelines([">" + pep_id + "\n"])
            f.writelines(h.readlines())
        with open(f"iupred/iupred_results/{pep_id}_short.txt", "r") as h:
            g.writelines([">" + pep_id + "\n"])
            g.writelines(h.readlines())

    f.close()
    g.close()

    f = open("full_prot_long.result", "w")
    g = open("full_prot_short.result", "w")

    prot_ids = list(zip(*pairs))[0]

    for prot_id in prot_ids:
        with open(f"iupred/iupred_results/{prot_id}_long.txt", "r") as h:
            f.writelines([">" + prot_id + "\n"])
            f.writelines(h.readlines())
        with open(f"iupred/iupred_results/{prot_id}_short.txt", "r") as h:
            g.writelines([">" + prot_id + "\n"])
            g.writelines(h.readlines())

    f.close()
    g.close()

    step3_iupred("peptides.fasta", "full_pep_short.result", "full_pep_long.result", "final_pickles/Peptide_Intrinsic_dict_v3")
    step3_iupred("proteins.fasta", "full_prot_short.result", "full_prot_long.result", "final_pickles/Protein_Intrinsic_dict")


def extract_intrinsic_disorder(fasta_filename, disorder_filename):
    raw_fasta_list = []
    with open('./' + fasta_filename, 'r') as f:
        for line in f.readlines():
            line_list = line.strip()
            raw_fasta_list.append(line_list)

    fasta_id_list = [x for x in raw_fasta_list if x[0] == '>']
    fasta_sequence_list = [x for x in raw_fasta_list if x[0] != '>']
    fasta_seq_len_list = [len(x) for x in fasta_sequence_list]

    print(len(fasta_id_list), len(fasta_sequence_list), len(fasta_seq_len_list))

    fasta_dict = {}
    for i in range(len(fasta_id_list)):
        fasta_dict[fasta_id_list[i]] = (fasta_sequence_list[i], fasta_seq_len_list[i])

    # load protein intrinsic disorder result
    raw_result_list = []
    with open(disorder_filename, 'r') as f:
        for line in f.readlines():
            line_list = line.strip()
            if (len(line_list) > 0 and line_list[0] != '#'):
                raw_result_list.append(line_list)

    intrinsic_id_list = [x for x in raw_result_list if x[0] == '>']
    intrinsic_score_list = [x.split('\t') for x in raw_result_list if x[0] != '>']
    print(np.array(intrinsic_score_list).shape)
    print(f"length of score list: {len(intrinsic_score_list)}")
    start_idx = 0
    raw_score_dict = {}
    for idx in range(len(intrinsic_id_list)):
        prot_id = intrinsic_id_list[idx]
        seq_len = fasta_dict[prot_id][1]
        end_idx = start_idx + seq_len
        individual_score_list = intrinsic_score_list[start_idx:end_idx]
        individual_score_list = [x[2:] for x in individual_score_list]
        individual_score_array = np.array(individual_score_list, dtype='float')
        raw_score_dict[prot_id] = individual_score_array
        start_idx = end_idx
    # print(len(fasta_dict.keys()),len(raw_score_dict.keys()))
    return fasta_dict, raw_score_dict


def step3_iupred(fasta_filename, iupred_filename_short, iupred_filename_long, output):
    # long & short
    fasta_dict_long, raw_score_dict_long = extract_intrinsic_disorder(fasta_filename, iupred_filename_long)
    fasta_dict_short, raw_score_dict_short = extract_intrinsic_disorder(fasta_filename, iupred_filename_short)

    Intrinsic_score_long = {}
    for key in fasta_dict_long.keys():
        sequence = fasta_dict_long[key][0]
        seq_len = fasta_dict_long[key][1]
        Intrinsic = raw_score_dict_long[key]
        if Intrinsic.shape[0] != seq_len:
            print(Intrinsic.shape[0], seq_len)
            print('Error!')
        Intrinsic_score_long[sequence] = Intrinsic

    Intrinsic_score_short = {}
    for key in fasta_dict_short.keys():
        sequence = fasta_dict_short[key][0]
        seq_len = fasta_dict_short[key][1]
        Intrinsic = raw_score_dict_short[key]
        if Intrinsic.shape[0] != seq_len:
            print('Error!')
        Intrinsic_score_short[sequence] = Intrinsic

    print(len(Intrinsic_score_short.keys()))

    Intrinsic_score = {}
    for seq in Intrinsic_score_short.keys():
        long_Intrinsic = Intrinsic_score_long[seq][:, 0]
        short_Intrinsic = Intrinsic_score_short[seq]
        concat_Intrinsic = np.column_stack((long_Intrinsic, short_Intrinsic))
        Intrinsic_score[seq] = np.column_stack((long_Intrinsic, short_Intrinsic))

    with open(output, 'wb') as f:  # 'output_intrisic_dict' is the name of the output dict you like
        pickle.dump(Intrinsic_score, f, protocol=2)
