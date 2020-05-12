#!/usr/bin/env python3

from raiden.read_params import Read_params
rp = Read_params('jiji')
args = rp.set_options()

import re
import os
import sys
import glob
import gzip
import collections
import subprocess as sbp
from raiden.utils import time_stamp, clean_cmd, call_log


class Jiji(object):

    def __init__(self, args):
        rp.check_args(args)
        self.args = args
        self.bed_files = sorted(glob.glob("{}/*.bed".format(self.args.bed)))
        self.N_bed_files = len(self.bed_files)
        self.gff_extension = self.args.gff.split('.')[-1]
        self.check_gff_extension()

        os.mkdir(self.args.out)

    def check_gff_extension(self):
        if self.gff_extension != 'gff' and self.gff_extension != 'gtf':
            print(time_stamp(), 
                  "!!WARNING!! {}'s extension is not 'gff' or 'gtf'\n".format(self.args.gff), 
                  flush=True)
            sys.exit(1)

    def get_transcript_region(self):
        input_annotation = open(self.args.gff)
        output_annotation = open('{0}/transcript.{1}'.format(self.args.out, self.gff_extension), 'w')

        not_comment = re.compile('[^#]')
        for line in input_annotation:
            if not_comment.match(line):
                cols = line.split('\t')

                if cols[2] == 'transcript':
                    output_annotation.write(line)

        input_annotation.close()
        output_annotation.close()

    def filter_VCF(self):
        not_comment = re.compile('[^#]')
        with gzip.open(self.args.vcf, 'rt') as vcf:
            with open('{0}/filtered_markers.bed'.format(self.args.out), 'w') as filtered_markers: 
                for line in vcf:
                    if not_comment.match(line):
                        cols = line.split('\t')

                        contig_name = cols[0]
                        position = int(cols[1])
                        ref = cols[3]
                        alt = cols[4]

                        GTs = [col.split(':')[0] for col in cols[9:]]

                        N_miss = GTs.count('./.')
                        if N_miss > self.args.miss:
                            continue

                        N_inconsistent = GTs.count('0/0')
                        if N_inconsistent > self.args.inconsistent:
                            continue

                        if not self.args.count_het:
                            N_inconsistent_het = sum([1 for GT in GTs if '0' in GT])
                            N_inconsistent += N_inconsistent_het
                            if N_inconsistent > self.args.inconsistent:
                                continue

                        filtered_markers.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(contig_name, 
                                                                                     position - 1,
                                                                                     position, 
                                                                                     ref, 
                                                                                     alt, 
                                                                                     N_miss, 
                                                                                     N_inconsistent))

    def check_mut_annotation(self):
        cmd = 'bedtools intersect -wa \
                                  -a {0}/transcript.{1} \
                                  -b {0}/filtered_markers.bed | \
               sort -u 1> {0}/candidate_genes_from_mutations.{1} \
                       2> {0}/jiji_mut_bedtools.log'.format(self.args.out, 
                                                            self.gff_extension)

        cmd = clean_cmd(cmd)

        try:
            sbp.run(cmd,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.args.out, 'jiji_mut_bedtools', cmd)
            sys.exit(1)

    def check_PA_coverage(self, index):
        cmd = 'bedtools coverage -a {0}/transcript.{1} \
                                 -b {2} \
               1> {0}/candidate_genes_from_PA.{3}.bed \
               2> {0}/jiji_PA_bedtools.log'.format(self.args.out,
                                                   self.gff_extension,
                                                   self.bed_files[index],
                                                   index)

        cmd = clean_cmd(cmd)

        try:
            sbp.run(cmd,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.args.out, 'jiji_PA_bedtools', cmd)
            sys.exit(1)

    def filter_and_count_PA(self):
        lines = []
        for i in range(self.N_bed_files):
            with open('{0}/candidate_genes_from_PA.{1}.bed'.format(self.args.out, i)) as bed:
                for line in bed:
                    line = line.rstrip('\n')
                    cols = line.split('\t')
                    coverage = float(cols[-1])

                    if coverage <= self.args.coverage:
                        line = '\t'.join(cols[:-4])
                        lines.append(line)

                os.remove('{0}/candidate_genes_from_PA.{1}.bed'.format(self.args.out, i))

        N_passed_PA = dict(collections.Counter(lines))

        with open('{0}/candidate_genes_from_PA.gff'.format(self.args.out), 'w') as gff:
            for k, v in N_passed_PA.items():
                if v >= self.N_bed_files - self.args.inconsistent:
                    gff.write('{}\n'.format(k))

    def remove_duplicates(self):
        cmd = 'cat {0}/candidate_genes_from_*.{1} | \
               cut -f 1-9 | \
               sort -u > {0}/all_candidate_genes.{1}'.format(self.args.out, 
                                                             self.gff_extension)

        cmd = clean_cmd(cmd)

        try:
            sbp.run(cmd,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            print(time_stamp(), 
                  '!!ERROR!! {}\n'.format(cmd), 
                  flush=True)
            sys.exit(1)

    def run(self):
        self.get_transcript_region()
        self.filter_VCF()
        self.check_mut_annotation()

        for i in range(self.N_bed_files):
            self.check_PA_coverage(i)

        self.filter_and_count_PA()
        self.remove_duplicates()



def main():
    print(time_stamp(), 'start to filter the causal genes.', flush=True)
    Jiji(args).run()
    print(time_stamp(), 'Filtering process successfully finished.\n', flush=True)

if __name__ == '__main__':
    main()