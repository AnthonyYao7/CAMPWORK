from subprocess import Popen
from contextlib import ExitStack
from sspro.step3_sspro import step3_sspro



def do_sspro(pairs, output_path):
    sspro_path = "sspro/SSpro_5.2/bin/"

    def kill(process):
        if process.poll() is None:
            process.kill()

    files = [f"peptide_fastas/{x}.fasta" for x in list(zip(*pairs))[1]] + [f"protein_fastas/{x}.fasta" for x in list(zip(*pairs))[0]]

    with ExitStack() as stack:
        processes = []
        for file in files:
            pdb_id = file.split('/')[1].split(".")[0]
            print(f"doing {pdb_id} ____________________________________________________")
            if "peptide" in file:
                command = f"{sspro_path}sequence_to_ss.sh {file} {output_path}/{pdb_id}_pep.out 4"
            else:
                command = f"{sspro_path}sequence_to_ss.sh {file} {output_path}/{pdb_id}_prot.out 4"
            processes.append(stack.enter_context(Popen(command, shell=True)))
            stack.callback(kill, processes[-1])
        for process in processes:
            process.wait()

    f = open("sspro/combined_prot.out.ss", "w")
    g = open("sspro/combined_pep.out.ss", "w")

    for pair in pairs:
        with open("sspro/sspro_outputs/{}_pep.out".format(pair[1]), "r") as h:
            g.writelines(h.readlines())

        with open("sspro/sspro_outputs/{}_prot.out".format(pair[0]), "r") as h:
            f.writelines(h.readlines())
    f.close()
    g.close()

    step3_sspro("sspro/combined_prot.out.ss", "proteins.fasta", "sspro/sspro_prot_output.csv")
    step3_sspro("sspro/combined_pep.out.ss", "peptides.fasta", "sspro/sspro_pep_output.csv")


