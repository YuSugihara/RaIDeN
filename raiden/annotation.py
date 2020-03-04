import os
import sys
import shutil
import subprocess as sbp
from raiden.utils import time_stamp, clean_cmd, call_log


class Annotation(object):

    def __init__(self, args):
        self.args = args

    def merge_bam(self):
        if len(self.args.rna_seq) > 1:
            cmd = 'samtools merge {0}/50_annotation/RNA-seq.bam \
                                  {0}/30_bam/RNA-seq.*.bam'.format(self.args.out)

            cmd = clean_cmd(cmd)

            try:
                sbp.run(cmd,
                        stdout=sbp.DEVNULL,
                        stderr=sbp.DEVNULL,
                        shell=True,
                        check=True)

            except sbp.CalledProcessError:
                call_log(self.args.out, 'samtools_merge', cmd)
                sys.exit(1)

        else:
            path_to_bam = os.path.abspath('{0}/30_bam/RNA-seq.0000.bam'.format(self.args.out))
            os.symlink(path_to_bam, 
                       '{0}/50_annotation/RNA-seq.bam'.format(self.args.out))

    def transciptome_assembly(self):
        if self.args.strand == 'None':
            cmd = 'stringtie -p {0} \
                             -m {1} \
                             -o {2}/50_annotation/annotation.gtf \
                             -l annotation \
                             -f 0.9 \
                             {2}/50_annotation/RNA-seq.bam \
                             > {2}/log/stringtie.log \
                             2>&1'.format(self.args.threads,
                                          self.args.minimum_len,
                                          self.args.out)
        else:
            cmd = 'stringtie -p {0} \
                             -m {1} \
                             -o {2}/50_annotation/annotation.gtf \
                             -l annotation \
                             --{3} \
                             -f 0.9 \
                             {2}/50_annotation/RNA-seq.bam \
                             > {2}/log/stringtie.log \
                             2>&1'.format(self.args.threads,
                                          self.args.minimum_len,
                                          self.args.out,
                                          self.args.strand)

        cmd = clean_cmd(cmd)

        try:
            sbp.run(cmd,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

        except sbp.CalledProcessError:
            call_log(self.args.out, 'stringtie', cmd)
            sys.exit(1)

    def make_gff(self):
        cmd = 'gffread -F \
                       -o {0}/50_annotation/annotation.gff \
                       {0}/50_annotation/annotation.gtf \
                       > {0}/log/gffread_gff.log \
                       2>&1'.format(self.args.ref, 
                                    self.args.out)

        cmd = clean_cmd(cmd)

        try:
            sbp.run(cmd,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

        except sbp.CalledProcessError:
            call_log(self.args.out, 'gffread_gff', cmd)
            sys.exit(1)

    def run_gffread(self):
        cmd = 'gffread -g {0} \
                       -w {1}/50_annotation/annotation.fasta \
                       {1}/50_annotation/annotation.gtf \
                       > {1}/log/gffread_fasta.log \
                       2>&1'.format(self.args.ref, 
                                    self.args.out)

        cmd = clean_cmd(cmd)

        try:
            sbp.run(cmd,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

        except sbp.CalledProcessError:
            call_log(self.args.out, 'gffread_fasta', cmd)
            sys.exit(1)

    def run(self):
        print(time_stamp(), 
              'start to annotate the reference genome.', 
              flush=True)
        self.merge_bam()
        self.transciptome_assembly()
        self.make_gff()
        self.run_gffread()