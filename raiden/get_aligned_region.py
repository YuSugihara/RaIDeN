import os
import sys
import subprocess as sbp
from raiden.utils import time_stamp, clean_cmd, call_log

class Get_aligned_region(object):

    def __init__(self, args):
        self.args = args

    def run_bamtools(self, index):
        cmd = 'bamtools filter -in {0}/30_bam/{1}.bam \
                               -out {0}/40_bed/{1}.no_error.bam \
                               -script {0}/log/bamtools.json \
                               > {0}/log/bamtools_{1}.log \
                               2>&1'.format(self.args.out,
                                            index)

        cmd = clean_cmd(cmd)

        try:
            sbp.run(cmd,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

        except sbp.CalledProcessError:
            call_log(self.args.out, 'bamtools_{}'.format(index), cmd)
            sys.exit(1)

    def run_bedtools(self, index):
        cmd = 'bedtools bamtobed -i {0}/40_bed/{1}.no_error.bam | \
               bedtools merge 1> {0}/40_bed/{1}.bed \
                              2> {0}/log/bedtools_{1}.log'.format(self.args.out,
                                                                  index)

        cmd = clean_cmd(cmd)

        try:
            sbp.run(cmd,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

        except sbp.CalledProcessError:
            call_log(self.args.out, 'bedtools_{}'.format(index), cmd)
            sys.exit(1)

    def run(self, index):
        self.run_bamtools(index)
        self.run_bedtools(index)

        os.remove('{0}/40_bed/{1}.no_error.bam'.format(self.args.out, 
                                                       index))