import re
import os
import sys
import glob
import subprocess as sbp
from multiprocessing import Pool
from raiden.utils import time_stamp, clean_cmd, call_log


class Mpileup(object):

    def __init__(self, args):
        self.args = args

    def get_header(self):
        ref = open(self.args.ref, 'r')
        pattern = re.compile('>')
        chr_names = []
        for line in ref:
            if pattern.match(line):
                line = line.rstrip('\n')
                chr_name = re.split('[> ]',line)[1]
                chr_names.append(chr_name)
        return chr_names

    def check_chr_name(self, chr_name):
        if ':' in chr_name:
            chr_name =  chr_name.replace(":", "_")
        return chr_name

    def mpileup(self, chr_name):
        cmd1 = 'bcftools mpileup -a AD,ADF,ADR \
                                 -B \
                                 -q {0}\
                                 -Q {1} \
                                 -O u \
                                 -r {2} \
                                 -f {3} \
                                 --ignore-RG \
                                 {4}/30_bam/WGS.*.bam | \
                bcftools call -vm \
                              -f GQ,GP \
                              -O u | \
                bcftools filter -i "INFO/MQ>={0}" \
                                -O z \
                                -o {4}/60_vcf/raiden.{5}.vcf.gz \
                                > {4}/log/bcftools.{5}.log \
                                2>&1'.format(self.args.min_MQ,
                                             self.args.min_BQ,
                                             chr_name,
                                             self.args.ref,
                                             self.args.out, 
                                             self.check_chr_name(chr_name))

        cmd2 = 'tabix -f \
                      -p vcf \
                      {0}/60_vcf/raiden.{1}.vcf.gz \
                      >> {0}/log/tabix.{1}.log \
                      2>&1'.format(self.args.out, 
                                   self.check_chr_name(chr_name))

        cmd1 = clean_cmd(cmd1)
        cmd2 = clean_cmd(cmd2)

        try:
            sbp.run(cmd1,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.args.out, 
                     'bcftools_{}'.fomrat(self.check_chr_name(chr_name)), 
                     cmd1)
            sys.exit(1)

        try:
            sbp.run(cmd2,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.args.out, 'tabix', cmd2)
            sys.exit(1)

    def concat(self):
        cmd1 = 'cat {0}/log/bcftools.*.log > {0}/log/bcftools.log'.format(self.args.out)
        cmd2 = 'cat {0}/log/tabix.*.log > {0}/log/tabix.log'.format(self.args.out)

        cmd3 = 'bcftools concat -a \
                                -O z \
                                -o {0}/60_vcf/raiden.vcf.gz \
                                {0}/60_vcf/raiden.*.vcf.gz \
                                >> {0}/log/bcftools.log \
                                2>&1'.format(self.args.out)

        cmd4 = 'rm -f {}/60_vcf/raiden.*.vcf.gz'.format(self.args.out)
        cmd5 = 'rm -f {}/60_vcf/raiden.*.vcf.gz.tbi'.format(self.args.out)
        cmd6 = 'rm -f {}/log/bcftools.*.log'.format(self.args.out)
        cmd7 = 'rm -f {}/log/tabix.*.log'.format(self.args.out)

        cmd1 = clean_cmd(cmd1)
        cmd2 = clean_cmd(cmd2)
        cmd3 = clean_cmd(cmd3)
        cmd4 = clean_cmd(cmd4)
        cmd5 = clean_cmd(cmd5)
        cmd6 = clean_cmd(cmd6)
        cmd7 = clean_cmd(cmd7)

        sbp.run(cmd1, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd2, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)

        try:
            sbp.run(cmd3,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.args.out, 'bcftools', cmd3)
            sys.exit(1)

        sbp.run(cmd4, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd5, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd6, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd7, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)

    def mkindex(self):
        cmd = 'tabix -f \
                     -p vcf \
                     {0}/60_vcf/raiden.vcf.gz \
                     >> {0}/log/tabix.log \
                     2>&1'.format(self.args.out)

        cmd = clean_cmd(cmd)

        try:
            sbp.run(cmd,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.args.out, 'tabix', cmd)
            sys.exit(1)    

    def run(self):
        print(time_stamp(), 'start to call variants.', flush=True)
        chr_names = self.get_header()

        p = Pool(self.args.threads)
        p.map(self.mpileup, chr_names)
        p.close()

        self.concat()
        self.mkindex()