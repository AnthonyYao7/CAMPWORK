import pandas as pd
import os


def aa_ss_concat(aa, ss):
    if len(aa) != len(ss):
        return 'string length error!'
    else:
        new_str = ''
        for i in range(len(aa)):
            concat_str = aa[i] + ss[i] + ','
            new_str = new_str + concat_str
    final_str = new_str[:-1]
    return final_str


def step3_sspro(combined, fasta, output):
    lines = []
    with open(combined, "r") as f:
        for line in f.readlines():
            lines.append([line.strip()])

    df_org = pd.DataFrame(lines)
    # df_org = pd.read_csv(combined, delimiter='>', header=None)  # the generated file by SCRATCH1D SSPro
    print(df_org)
    df_org.columns = ['col_1']

    ss_idx = []
    seq_idx = []
    seq_idx = [2 * x for x in list(range(int(df_org.shape[0] / 2)))]
    ss_idx = [x + 1 for x in seq_idx]

    # subset sequence dataframe and sse dataframe
    df_seq = df_org.iloc[seq_idx]
    df_seq.columns = ['seq_id']
    df_ss = df_org.iloc[ss_idx]
    df_ss.columns = ['seq_ss']

    df_seq = df_seq.reset_index(drop=True)
    df_ss = df_ss.reset_index(drop=True)

    # join sequence & sse together
    df_seq_ss = pd.merge(df_seq, df_ss, left_index=True, right_index=True)

    # load id mapping file
    df_id = pd.read_csv(fasta, sep='\t', header=None)  # the input asta file used for SCRATCH1D SSPro
    df_id.columns = ['col_1']

    ss_idx = []
    seq_idx = []
    seq_idx = [2 * x for x in list(range(int(df_id.shape[0] / 2)))]
    ss_idx = [x + 1 for x in seq_idx]

    # subset sequence dataframe and sse dataframe
    df_seq = df_id.iloc[seq_idx]
    df_seq.columns = ['seq_id']
    df_ss = df_id.iloc[ss_idx]
    df_ss.columns = ['seq']

    df_seq = df_seq.reset_index(drop=True)
    df_ss = df_ss.reset_index(drop=True)

    # join sequence &  sse together
    df_idx = pd.merge(df_seq, df_ss, left_index=True, right_index=True)
    df_output_ss = pd.merge(df_idx, df_seq_ss, left_on=['seq_id'], right_on=['seq_id'])
    df_output_ss['concat_seq'] = df_output_ss.apply(lambda x: aa_ss_concat(x['seq'], x['seq_ss']), axis=1)
    df_output_ss.to_csv(output, encoding='utf-8', index=False,
                        sep='\t')  # 'output_ss_filename' is the name of the output tsv you like
