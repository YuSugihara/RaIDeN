#!/usr/bin/env python3

from raiden.read_params import Read_params
rp = Read_params('raiden')
args = rp.set_options()

import os
import sys
import glob
import shutil
import subprocess as sbp
from multiprocessing import Pool
from raiden.index_ref import Index_ref
from raiden.alignment import Alignment
from raiden.annotation import Annotation
from raiden.mpileup import Mpileup
from raiden.utils import time_stamp, clean_cmd, call_log, get_proc_numbers


class RaIDeN(object):

    def __init__(self, args):
        rp.check_args(args)
        args = rp.check_max_threads(args)
        self.args = args

        self.mkdir()
        self.symlink_ref()

    def mkdir(self):
        os.mkdir('{}'.format(self.args.out))
        os.mkdir('{}/log'.format(self.args.out))
        os.mkdir('{}/10_ref'.format(self.args.out))

    def symlink_ref(self):
        path_to_ref = os.path.abspath(self.args.ref)
        ref = os.path.basename(self.args.ref)
        sym_ref = '{}/10_ref/{}'.format(self.args.out, ref)
        os.symlink(path_to_ref, sym_ref)
        self.args.ref = sym_ref

    def index_ref(self):
        ir = Index_ref(args)
        ir.run()

    def alignment(self):
        aln = Alignment(args)
        aln_args = []

        N_files = len(self.args.rna_seq) + len(self.args.whole_genome)
        N_process, each_threads = get_proc_numbers(self.args.threads, N_files)

        for i, fastq in enumerate(self.args.rna_seq):
            fastq1 = fastq.split(',')[0]
            fastq2 = fastq.split(',')[1]
            index = 'RNA-seq.0{:0>3}'.format(i)
            aln_arg = '\t'.join([fastq1, fastq2, index, 'RNA'])
            aln_args.append(aln_arg)

            print(time_stamp(),
                  "{}'s prefix -> {}.".format(fastq1, index),
                  flush=True)

            print(time_stamp(),
                  "{}'s prefix -> {}.".format(fastq2, index),
                  flush=True)

        for i, fastq in enumerate(self.args.whole_genome):
            fastq1 = fastq.split(',')[0]
            fastq2 = fastq.split(',')[1]
            index = 'WGS.0{:0>3}'.format(i)
            aln_arg = '\t'.join([fastq1, fastq2, index, 'Genome'])
            aln_args.append(aln_arg)

            print(time_stamp(),
                  "{}'s prefix -> {}.".format(fastq1, index),
                  flush=True)

            print(time_stamp(),
                  "{}'s prefix -> {}.".format(fastq2, index),
                  flush=True)

        for i in range(N_files):
            aln_args[i] = aln_args[i] + '\t' + str(each_threads[i])

        p = Pool(N_process)
        p.map(aln.run, aln_args)
        p.close()

        print(time_stamp(),
              'alignment successfully finished.',
              flush=True)

    def annotation(self):
        os.mkdir('{}/50_annotation'.format(self.args.out))
        ann = Annotation(self.args)
        ann.run()

    def mpileup(self):
        os.mkdir('{}/60_vcf'.format(self.args.out))
        mp = Mpileup(self.args)
        mp.run()

    def filter_candidates(self):
        cmd = 'jiji -a {0}/50_annotation/annotation.gff \
                    -b {0}/40_bed \
                    -v {0}/60_vcf/raiden.vcf.gz \
                    -o {0}/70_result'.format(self.args.out)

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

        shutil.move('{0}/70_result/jiji_PA_bedtools.log'.format(self.args.out), 
                    '{0}/log/'.format(self.args.out))
        shutil.move('{0}/70_result/jiji_mut_bedtools.log'.format(self.args.out), 
                    '{0}/log/'.format(self.args.out))

    def run(self):
        self.index_ref()
        self.alignment()
        self.annotation()
        self.mpileup()
        self.filter_candidates()




def main():
    print(time_stamp(), 'start to run RaIDeN.', flush=True)
    RaIDeN(args).run()
    print(time_stamp(), 'RaIDeN successfully finished.\n', flush=True)

if __name__ == '__main__':
    main()