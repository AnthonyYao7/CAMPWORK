import os
import argparse
import pandas as pd
from iupred.intrinsic_disorder import run_iupred
from sspro.sspro import do_sspro


base_dir = os.path.basename(os.path.dirname(__file__))

def main():
    parser = argparse.ArgumentParser(description="Run inference on peptide protein pairs")
    parser.add_argument('ppfile',
                        help='path to csv file with headers [optional]<pep_id> <pep_seq> [optional]<prot_id> <prot_seq>. To omit an id, leave blank')
    args = parser.parse_args()

    assert os.path.isfile(args.ppfile), "could not find csv file"
    # preprocesses the pep-prot file
    os.system("python preprocess.py {}".format(args.ppfile))

    df = pd.read_csv(args.ppfile)
    pairs = list(zip(df['prot_id'].to_list(), df['pep_id'].to_list()))
    print(pairs)
    # run_iupred(pairs)
    do_sspro(pairs, "sspro/sspro_outputs")


if __name__ == "__main__":
    main()