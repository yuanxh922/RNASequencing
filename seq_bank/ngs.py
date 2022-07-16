""" main entrance of seq bank
"""
import sys
import os
from datetime import datetime
import click
import json
import pandas as pd
from loguru import logger
from seq import SeqManager, Compound
# PARENT_DIR = os.path.dirname(os.getcwd())
# print(os.getcwd())
# sys.path.append(PARENT_DIR)

def run(seqs, fpath):
    df = pd.read_excel(fpath, 2)
    for idx, seq in df.seqs.iteritems():
        if len(seq) < 32:
            continue
        print(seq)
        seq = seq + 'CCA'

        df_5p = generate_mass_ladder_5p(seq)
        df_3p = generate_mass_ladder_3p(seq)
    return 'aaa'
    # return generate_mass_ladder_5p(seqs[0])

def generate_mass_ladder_5p(seq):
    seq = T2U(seq)
    anchor = 18.0106 + 79.9663

    sm = SeqManager()
    _, cpds = sm.generate_cpds_5p(seq, anchor)
    df = Compound.cpds2df(cpds)
    # print(df.info())
    return df

def generate_mass_ladder_3p(seq):
    seq = T2U(seq)
    rev_anchor = -61.9557

    sm = SeqManager()
    _, cpds = sm.generate_cpds_3p(seq, rev_anchor)
    df = Compound.cpds2df(cpds)
    # print(df.info())
    return df

def T2U(seq):
    """Change input sequence from DNA to RNA stype, 'T' to 'U' 
    """
    # print('input', seq)
    seq = seq.replace('T', 'U').replace('t', 'U')
    # print('output', seq)
    return seq

@click.command()
@click.option("--seqs", "-s", type=str, help="sequences.")
@click.option("--fpath", "-f", type=str, help="path to an Excel file which contains a list of sequences.")
def action(seqs, fpath):
    seqs = seqs.split(',')
    # print(seqs)
    return run(seqs, fpath)

if __name__ == "__main__":
    action()
