import os
import sys
import subprocess as sbp
from raiden.prinseq import Prinseq
from raiden.gzip_fastq import Gzip_fastq
from raiden.get_aligned_region import Get_aligned_region
from raiden.utils import time_stamp, clean_cmd, call_log


class Alignment(object):

    def __init__(self, args):
        self.args = args

        if self.args.disable_RNAseq_trim:
             print(time_stamp(),
                  'disable the trimming of RNA-seq.',
                  flush=True)

        if self.args.disable_WGS_trim:
            print(time_stamp(),
                  'disable the trimming of WGS.',
                  flush=True)

        if self.args.disable_RNAseq_trim and \
           self.args.disable_WGS_trim:
            print(time_stamp(),
                  'start to align reads.',
                  flush=True)
        else:
            print(time_stamp(),
                  'start to trim and align reads.',
                  flush=True)

        os.mkdir('{}/20_fastq'.format(self.args.out))
        os.mkdir('{}/30_bam'.format(self.args.out))

        self.write_json()
        os.mkdir('{0}/40_bed/'.format(self.args.out))

    def write_json(self):
        with open('{0}/log/bamtools.json'.format(self.args.out), 'w') as json:
            json.write(('{"tag" : "NM:<1"}\n' 
                        '###NM=all_miss\n'
                        '###XM=Missonly\n'))

    def trim_RNAseq(self, fastq1, fastq2, index, N_threads):
        cmd = 'FaQCs -1 {0} \
                     -2 {1} \
                     --prefix {2} \
                     -d {3}/20_fastq/FaQCs_{2} \
                     -t {4} \
                     -min_L 50 \
                     -avg_q 20 \
                     --polyA \
                     --adapter \
                     -discard 1 \
                     > {3}/log/FaQCs_{2}.log \
                     2>&1'.format(fastq1,
                                  fastq2,
                                  index, 
                                  self.args.out,
                                  N_threads)

        cmd = clean_cmd(cmd)

        try:
            sbp.run(cmd,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

        except sbp.CalledProcessError:
            call_log(self.args.out, 'FaQCs_{}'.format(index), cmd)
            sys.exit(1)

        fastq1 = "{0}/20_fastq/FaQCs_{1}/{1}.1.trimmed.fastq".format(self.args.out, index)
        fastq2 = "{0}/20_fastq/FaQCs_{1}/{1}.2.trimmed.fastq".format(self.args.out, index)

        return fastq1, fastq2

    def trim_WGS(self, fastq1, fastq2, index, N_threads):
        cmd = 'FaQCs -1 {0} \
                     -2 {1} \
                     --prefix {2} \
                     -d {3}/20_fastq/FaQCs_{2} \
                     -t {4} \
                     --adapter \
                     > {3}/log/FaQCs_{2}.log \
                     2>&1'.format(fastq1,
                                  fastq2,
                                  index, 
                                  self.args.out,
                                  N_threads)

        cmd = clean_cmd(cmd)

        try:
            sbp.run(cmd,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

        except sbp.CalledProcessError:
            call_log(self.args.out, 'FaQCs_{}'.format(index), cmd)
            sys.exit(1)

        prins = Prinseq(self.args, N_threads, index)
        prins.run()

        fastq1 = "{0}/20_fastq/prinseq_{1}/{1}_1.fastq".format(self.args.out, index)
        fastq2 = "{0}/20_fastq/prinseq_{1}/{1}_2.fastq".format(self.args.out, index)

        return fastq1, fastq2

    def align_RNAseq(self, fastq1, fastq2, index, N_threads):
        cmd1 = 'hisat2 -1 {0} \
                       -2 {1} \
                       -x {2} \
                       --no-mixed \
                       --no-discordant \
                       -k 1 \
                       -p {3} \
                       2> {4}/log/hisat2_{5}.log | \
                samtools view -b \
                              -F 004 | \
                samtools sort -@ {3} \
                              -o {4}/30_bam/{5}.bam'.format(fastq1,
                                                            fastq2,
                                                            self.args.ref,
                                                            N_threads,
                                                            self.args.out,
                                                            index)

        cmd2 = 'samtools index {0}/30_bam/{1}.bam \
                               2> {0}/log/samtools_index_{1}.log'.format(self.args.out, 
                                                                         index)

        cmd1 = clean_cmd(cmd1)
        cmd2 = clean_cmd(cmd2)

        try:
            sbp.run(cmd1,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

        except sbp.CalledProcessError:
            call_log(self.args.out, 'hisat2_{}'.format(index), cmd1)
            sys.exit(1)

        try:
            sbp.run(cmd2,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

        except sbp.CalledProcessError:
            call_log(self.args.out, 'samtools_index_{}'.format(index), cmd2)
            sys.exit(1)

    def align_WGS(self, fastq1, fastq2, index, N_threads):
        cmd1 = 'hisat2 -1 {0} \
                       -2 {1} \
                       -x {2} \
                       --no-mixed \
                       --no-discordant \
                       --no-spliced-alignment \
                       -k 1 \
                       -p {3} \
                       2> {4}/log/hisat2_{5}.log | \
                samtools view -b \
                             -F 004 | \
                samtools sort -@ {3} \
                              -o {4}/30_bam/{5}.bam'.format(fastq1,
                                                            fastq2,
                                                            self.args.ref,
                                                            N_threads,
                                                            self.args.out,
                                                            index)

        cmd2 = 'samtools index {0}/30_bam/{1}.bam \
                               2> {0}/log/samtools_index_{1}.log'.format(self.args.out, 
                                                                         index)

        cmd1 = clean_cmd(cmd1)
        cmd2 = clean_cmd(cmd2)

        try:
            sbp.run(cmd1,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

        except sbp.CalledProcessError:
            call_log(self.args.out, 'hisat2_{}'.format(index), cmd1)
            sys.exit(1)

        try:
            sbp.run(cmd2,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

        except sbp.CalledProcessError:
            call_log(self.args.out, 'samtools_index_{}'.format(index), cmd2)
            sys.exit(1)

    def run(self, aln_arg):
        aln_arg = aln_arg.split('\t')
        fastq1 = aln_arg[0]
        fastq2 = aln_arg[1]
        index = aln_arg[2]
        mode = aln_arg[3]
        N_threads = aln_arg[4]

        gf = Gzip_fastq(self.args, N_threads, index)
        gar = Get_aligned_region(self.args)

        if mode == 'RNA':
            if not self.args.disable_RNAseq_trim:
                fastq1, fastq2 = self.trim_RNAseq(fastq1, 
                                                  fastq2, 
                                                  index, 
                                                  N_threads)

            self.align_RNAseq(fastq1, fastq2, index, N_threads)
            gf.gzip_FaQCs()

        if mode == 'Genome':
            if not self.args.disable_WGS_trim:
                fastq1, fastq2 = self.trim_WGS(fastq1, 
                                               fastq2, 
                                               index, 
                                               N_threads)

            self.align_WGS(fastq1, fastq2, index, N_threads)
            gf.gzip_FaQCs()
            gf.gzip_prinseq()
            gar.run(index)
