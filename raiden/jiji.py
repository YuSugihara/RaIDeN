#!/usr/bin/env python3

from raiden.read_params import Read_params
rp = Read_params('jiji')
args = rp.set_options()

import re
import os
import sys
import gzip
import subprocess as sbp
from raiden.utils import time_stamp, clean_cmd, call_log


class Jiji(object):

    def __init__(self, args):
        rp.check_args(args)
        self.args = args
        self.file_extension = self.args.gff.split('.')[-1]

        os.mkdir(self.args.out)

    def filter_VCF(self):
        not_comment = re.compile('[^#]')
        with gzip.open(self.args.vcf, 'rt') as vcf:
            with open('{}/filtered_mut.bed'.format(self.args.out), 'w') as filtered_mut: 
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

                        filtered_mut.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(contig_name, 
                                                                                 position - 1,
                                                                                 position, 
                                                                                 ref, 
                                                                                 alt, 
                                                                                 N_miss, 
                                                                                 N_inconsistent))

    def check_mut_annotation(self):
        file_extension = self.args.gff.split('.')[-1]

        if self.file_extension == 'gff' or self.file_extension == 'gtf':
            cmd = (r"awk -F '\t' "
                   "'{{if ($3=="
                   '"transcript")'
                   "{{print $0}}}}' {0} | \
                    bedtools intersect -wa -a stdin -b {1}/filtered_mut.bed | \
                    sort -u 1> {1}/candidate_genes_from_mutations.{2} \
                            2> {1}/jiji_mut_bedtools.log").format(self.args.gff, 
                                                                  self.args.out, 
                                                                  self.file_extension)

        else:
            print(time_stamp(), 
                  "!!WARNING!! {}'s extension is not 'gff' or 'gtf'\n".format(self.args.gff), 
                  flush=True)
            sys.exit(1)

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

    def check_PA_annotation(self):
        if self.file_extension == 'gff' or self.file_extension == 'gtf':
            cmd = (r"awk -F '\t' "
                   "'{{if ($3=="
                   '"transcript")'
                   "{{print $0}}}}' {0} | \
                   bedtools coverage -hist -a stdin -b {1}/*.bed | "\
                   r"awk -F '\t' "
                   "'{{if ($10<={2} && NF==13){{print $0}}}}' | \
                   sort -u 1> {3}/candidate_genes_from_PA.temp \
                           2> {3}/jiji_PA_bedtools.log").format(self.args.gff, 
                                                                self.args.bed, 
                                                                self.args.inconsistent,
                                                                self.args.out)

        else:
            print(time_stamp(), 
                  "!!WARNING!! {}'s extension is not 'gff' or 'gtf'\n".format(self.args.gff), 
                  flush=True)
            sys.exit(1)

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

    def check_PA_coverage(self):
        temp = open('{0}/candidate_genes_from_PA.temp'.format(self.args.out))
        PA_gff = open('{0}/candidate_genes_from_PA.{1}'.format(self.args.out, 
                                                               self.file_extension), 'w')

        previous_attribute = None
        cumulative_coverage = 0
        for line in temp:
            line = line.rstrip('\n')
            cols = line.split('\t')

            attribute = cols[8]
            feature_coverage = float(cols[12])

            if (previous_attribute != None) and (previous_attribute != attribute):
                if cumulative_coverage <= self.args.coverage:
                    PA_gff.write('{}\n'.format('\t'.join(cols[:9])))
                cumulative_coverage = 0
            else:
                cumulative_coverage += feature_coverage

            previous_attribute = attribute

        os.remove('{0}/candidate_genes_from_PA.temp'.format(self.args.out))

    def remove_duplicates(self):
        cmd = 'cat {0}/candidate_genes_from_*.{1} | \
               cut -f 1-9 | \
               sort -u > {0}/all_candidate_genes.{1}'.format(self.args.out, 
                                                             self.file_extension)

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
        self.filter_VCF()
        self.check_mut_annotation()
        self.check_PA_annotation()
        self.check_PA_coverage()
        self.remove_duplicates()



def main():
    print(time_stamp(), 'start to filter the causal genes.', flush=True)
    Jiji(args).run()
    print(time_stamp(), 'Filtering process successfully finished.\n', flush=True)

if __name__ == '__main__':
    main()