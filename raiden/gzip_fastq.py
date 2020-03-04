import os
import sys
import subprocess as sbp
from raiden.utils import time_stamp, clean_cmd, call_log


class Gzip_fastq(object):

    def __init__(self, args, N_threads, index):
        self.args = args
        self.N_threads = N_threads
        self.index = index

    def gzip_FaQCs(self):
        cmd1 = 'pigz -p {0} \
                     {1}/20_fastq/FaQCs_{2}/{2}.1.trimmed.fastq'.format(self.N_threads, 
                                                                        self.args.out, 
                                                                        self.index)
        
        cmd2 = 'pigz -p {0} \
                     {1}/20_fastq/FaQCs_{2}/{2}.2.trimmed.fastq'.format(self.N_threads, 
                                                                        self.args.out, 
                                                                        self.index)

        cmd3 = 'pigz -p {0} \
                     {1}/20_fastq/FaQCs_{2}/{2}.discard.trimmed.fastq'.format(self.N_threads, 
                                                                              self.args.out, 
                                                                              self.index)

        cmd4 = 'pigz -p {0} \
                     {1}/20_fastq/FaQCs_{2}/{2}.unpaired.trimmed.fastq'.format(self.N_threads, 
                                                                               self.args.out, 
                                                                               self.index)

        cmd1 = clean_cmd(cmd1)
        cmd2 = clean_cmd(cmd2)
        cmd3 = clean_cmd(cmd3)
        cmd4 = clean_cmd(cmd4)

        if os.path.isfile('{0}/20_fastq/FaQCs_{1}/{1}.1.trimmed.fastq'.format(self.args.out, 
                                                                              self.index)):
            try:
                sbp.run(cmd1,
                        stdout=sbp.DEVNULL,
                        stderr=sbp.DEVNULL,
                        shell=True,
                        check=True)

            except sbp.CalledProcessError:
                print(time_stamp(), 
                    '!!ERROR!! {}\n'.format(cmd1), 
                    flush=True)
                sys.exit(1)

        if os.path.isfile('{0}/20_fastq/FaQCs_{1}/{1}.2.trimmed.fastq'.format(self.args.out, 
                                                                              self.index)):
            try:
                sbp.run(cmd2,
                        stdout=sbp.DEVNULL,
                        stderr=sbp.DEVNULL,
                        shell=True,
                        check=True)

            except sbp.CalledProcessError:
                print(time_stamp(), 
                    '!!ERROR!! {}\n'.format(cmd2), 
                    flush=True)
                sys.exit(1)

        if os.path.isfile('{0}/20_fastq/FaQCs_{1}/{1}.discard.trimmed.fastq'.format(self.args.out, 
                                                                                    self.index)):
            try:
                sbp.run(cmd3,
                        stdout=sbp.DEVNULL,
                        stderr=sbp.DEVNULL,
                        shell=True,
                        check=True)

            except sbp.CalledProcessError:
                print(time_stamp(), 
                    '!!ERROR!! {}\n'.format(cmd3), 
                    flush=True)
                sys.exit(1)

        if os.path.isfile('{0}/20_fastq/FaQCs_{1}/{1}.unpaired.trimmed.fastq'.format(self.args.out, 
                                                                                     self.index)):
            try:
                sbp.run(cmd4,
                        stdout=sbp.DEVNULL,
                        stderr=sbp.DEVNULL,
                        shell=True,
                        check=True)

            except sbp.CalledProcessError:
                print(time_stamp(), 
                    '!!ERROR!! {}\n'.format(cmd4), 
                    flush=True)
                sys.exit(1)

    def gzip_prinseq(self):
        cmd1 = 'pigz -p {0} \
                        {1}/20_fastq/prinseq_{2}/{2}_1.fastq'.format(self.N_threads, 
                                                                     self.args.out, 
                                                                     self.index)
        
        cmd2 = 'pigz -p {0} \
                        {1}/20_fastq/prinseq_{2}/{2}_2.fastq'.format(self.N_threads, 
                                                                     self.args.out, 
                                                                     self.index)

        cmd3 = 'pigz -p {0} \
                        {1}/20_fastq/prinseq_{2}/{2}_1_singletons.fastq'.format(self.N_threads, 
                                                                                self.args.out, 
                                                                                self.index)

        cmd4 = 'pigz -p {0} \
                        {1}/20_fastq/prinseq_{2}/{2}_2_singletons.fastq'.format(self.N_threads, 
                                                                                self.args.out, 
                                                                                self.index)

        cmd1 = clean_cmd(cmd1)
        cmd2 = clean_cmd(cmd2)
        cmd3 = clean_cmd(cmd3)
        cmd4 = clean_cmd(cmd4)

        try:
            sbp.run(cmd1,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

        except sbp.CalledProcessError:
            print(time_stamp(), 
                  '!!ERROR!! {}\n'.format(cmd1), 
                  flush=True)
            sys.exit(1)

        try:
            sbp.run(cmd2,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

        except sbp.CalledProcessError:
            print(time_stamp(), 
                  '!!ERROR!! {}\n'.format(cmd2), 
                  flush=True)
            sys.exit(1)

        try:
            sbp.run(cmd3,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

        except sbp.CalledProcessError:
            print(time_stamp(), 
                  '!!ERROR!! {}\n'.format(cmd3), 
                  flush=True)
            sys.exit(1)

        try:
            sbp.run(cmd4,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

        except sbp.CalledProcessError:
            print(time_stamp(), 
                  '!!ERROR!! {}\n'.format(cmd4), 
                  flush=True)
            sys.exit(1)