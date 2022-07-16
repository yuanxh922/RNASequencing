""" main entrance of seq bank
"""
import sys
import os
from datetime import datetime
import click
import json
from loguru import logger
from seq import SeqManager
PARENT_DIR = os.path.dirname(os.getcwd())
print(os.getcwd())
sys.path.append(PARENT_DIR)

def run(anchor, orientation, length, width, num, output_path, add_noises, seqs):
    """ main entrance
    """
    if not isinstance(anchor, float) or not isinstance(length, int):
        logger.info("wrong input type")
        return

    sm = SeqManager()
    symbol_seqs, cpds = sm.make_cpds(orientation, length, width, num, add_noises, seqs)
    logger.info("got cpds size {}", len(cpds))

    now = datetime.now()
    fpath = "{}_{}_{}.xlsx".format(anchor, length, now.strftime('%Y%m%d%H%M%S'))
    if output_path:
        fpath = os.path.join(output_path, fpath)
    print("generated excel {}".format(fpath))
    cpds = sorted(cpds, key=lambda cpd: cpd.mass)
    df = sm.to_excel(cpds, fpath)

    json_path = fpath[:-4] + "json"
    with open(json_path, 'w') as f:
        json.dump(sm.sequences, f, indent=4)
    print("generated json {}".format(json_path))

def usage(cmd):
    """ prompt usage messages
    """
    print("Usage:")
    print("     {} <anchor> <length>".format(cmd))

@click.command()
@click.option("--anchor", "-a", type=float, default=694.2397, help="Anchor point.")
@click.option("--orientation", "-o", help="Reading orientation.")
@click.option("--length", "-l", default=22, help="Length of sequence.")
@click.option("--width", "-w", default=0, help="Length range of sequence.")
@click.option("--num", "-n", default=1, help="Number of sequences.")
@click.option("--output_dir", "-d", help="output directory")
@click.option("--noises", "-x", type=bool, default=False, help="Add noises.")
@click.option("--seqs", "-s", type=str, help="sequences.")
def action(anchor, orientation, length, width, num, output_dir, noises, seqs):
    if orientation:
        orientation = int(orientation)
    if seqs:
        seqs = seqs.split(',')
    print(seqs)
    run(anchor, orientation, length, width, num, output_dir, noises, seqs)

if __name__ == "__main__":
    action()
