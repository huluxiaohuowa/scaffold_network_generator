from data import *

with open('data/datasets/scaffolds.smi') as f:
    with open('data/datasets/c_scaffolds.smi', 'w') as w:
        for i, line in enumerate(f.readlines()):
            smi = line.strip().split('\t')[0]
            w.write(
                con_smi_to_c_smi(smi) +
                '\n'
            )
        w.flush()
