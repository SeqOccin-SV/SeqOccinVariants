#!/usr/bin/env python3

# Import modules here
import sys, os # for example
import argparse

import numpy as np
import pandas as pd
from pathlib import Path

def eprint(*args, **kwargs):
        print(*args,  file=sys.stderr, **kwargs)

def main(sourcedir, samples, output):
    samples = pd.read_csv(samples, names=["sample", "ID", "coverage", "size"],
                          dtype={'coverage': np.int32, 'size': np.int32},
                          skiprows=1, sep="\t")

    gsamples = samples.groupby(by=['sample', 'ID']).agg({'coverage': 'sum', 'size': 'sum'},
                               as_index=False).reset_index()
    animals = []
    reads = []
    for index, row in gsamples.iterrows():
        animal = str(int(row['ID'])) + "_" + row['sample']
        animaldir=os.path.join(sourcedir, animal)
        if not os.path.isdir(animaldir):
            eprint("WARNING %s not and existing dir for %s" % (animaldir, animal))
        subreads = [str(path) for path in Path(animaldir).rglob('*.subreads.bam')]
        animals.append(animal)
        reads.append(",".join(subreads))
    df = pd.DataFrame({"sample":animals, "path": reads})
    df.to_csv(output, sep="\t", index=False)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Adding information to paftools vcf')
    parser.add_argument('-d', '--dir',
                        required=True, help='parental dir')
    parser.add_argument('-s', '--samples',
                        required=True, help='sample tabular file dir')
    parser.add_argument('-o', '--output',
                        required=True, help='output samples.tsv file')

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main(args.dir, args.samples, args.output)
