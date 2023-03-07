import pandas as pd
import argparse
import os


base_dir = os.path.basename(os.path.dirname(__file__))
file_dir = os.path.dirname(__file__)
cwd = os.getcwd()


def preprocess(ppfile):
    directories_to_clear = [
        "final_pickles",
        "iupred/iupred_results",
        "peptide_fastas",
        "protein_fastas",
    ]

    for directory in directories_to_clear:
        os.system(f"rm {directory}/*")

    df = pd.read_csv(ppfile)

    g = open(f"peptides.fasta", "w")
    h = open(f"proteins.fasta", "w")

    for row in df.iterrows():
        with open(f"peptide_fastas/{row[1]['pep_id']}.fasta", "w") as f:
            f.write(f">{row[1]['pep_id']}\n{row[1]['pep_seq']}\n")
            g.write(f">{row[1]['pep_id']}\n{row[1]['pep_seq']}\n")
        with open(f"protein_fastas/{row[1]['prot_id']}.fasta", "w") as f:
            f.write(f">{row[1]['prot_id']}\n{row[1]['prot_seq']}\n")
            h.write(f">{row[1]['prot_id']}\n{row[1]['prot_seq']}\n")

    g.close()
    h.close()


def main():
    parser = argparse.ArgumentParser(description="Run inference on peptide protein pairs")
    parser.add_argument('ppfile',
                        help='path to csv file with headers [optional]<pep_id> <pep_seq> [optional]<prot_id> <prot_seq>. To omit an id, leave blank')
    args = parser.parse_args()

    assert os.path.isfile(args.ppfile), "could not find csv file"
    preprocess(args.ppfile)



if __name__ == "__main__":
    main()
