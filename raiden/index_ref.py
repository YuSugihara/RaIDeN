import sys
import subprocess as sbp
from raiden.utils import time_stamp, clean_cmd, call_log


class Index_ref(object):

    def __init__(self, args):
        self.args = args

    def run(self):
        print(time_stamp(),
              'start to index reference fasta.',
              flush=True)

        cmd1 = 'hisat2-build -p {0} {1} {1} \
                > {2}/log/hisat2-build.log \
                2>&1'.format(self.args.threads, self.args.ref, self.args.out)

        cmd2 = 'samtools faidx {} \
                > {}/log/samtools_faidx.log \
                2>&1'.format(self.args.ref, self.args.out)

        cmd1 = clean_cmd(cmd1)
        cmd2 = clean_cmd(cmd2)

        print(time_stamp(),
              'hisat2-build...',
              flush=True)

        try:
            sbp.run(cmd1,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

        except sbp.CalledProcessError:
            call_log(self.args.out, 'hisat2-build', cmd1)
            sys.exit(1)

        try:
            sbp.run(cmd2,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

        except sbp.CalledProcessError:
            call_log(self.args.out, 'samtools_faidx', cmd2)
            sys.exit(1)

        print(time_stamp(),
              'indexing of the reference genome successfully finished.',
              flush=True)