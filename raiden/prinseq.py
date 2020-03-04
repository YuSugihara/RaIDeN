import os
import sys
import shutil
import subprocess as sbp
from multiprocessing import Pool
from raiden.utils import time_stamp, clean_cmd, call_log


class Prinseq(object):

    def __init__(self, args, N_threads, index):
        self.args = args
        self.N_threads = int(N_threads)
        self.index = index

        os.mkdir('{0}/20_fastq/prinseq_{1}'.format(self.args.out, 
                                                   self.index))

    def seqkit_split2(self):
        cmd = 'seqkit split2 -p {0} \
                             -j {0} \
                             -1 {1}/20_fastq/FaQCs_{2}/{2}.1.trimmed.fastq \
                             -2 {1}/20_fastq/FaQCs_{2}/{2}.2.trimmed.fastq \
                             -o {1}/20_fastq/FaQCs_{2} \
                             > {1}/log/seqkit_{2}.log \
                             2>&1'.format(self.N_threads, 
                                          self.args.out, 
                                          self.index)

        cmd = clean_cmd(cmd)

        try:
            sbp.run(cmd,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)

        except sbp.CalledProcessError:
            call_log(self.args.out, 'seqkit_{}'.format(self.index), cmd)
            sys.exit(1)

    def run_prinseq(self):
        cmd = 'seq -f %03g {0} | \
               xargs -P {0} \
                     -I % \
               prinseq-lite.pl -trim_left 5 \
                               -trim_right 20 \
                               -trim_qual_window 10 \
                               -trim_qual_right 20 \
                               -min_len 75 \
                               -min_qual_mean 20 \
                               -fastq {1}/20_fastq/FaQCs_{2}/{2}.1.trimmed.fastq.split/{2}.1.trimmed.part_%.fastq \
                               -fastq2 {1}/20_fastq/FaQCs_{2}/{2}.1.trimmed.fastq.split/{2}.2.trimmed.part_%.fastq \
                               -out_good {1}/20_fastq/prinseq_{2}/{2}.part_% \
                               -out_bad null \
                               > {1}/log/prinseq_{2}.log \
                               2>&1'.format(self.N_threads, 
                                            self.args.out,
                                            self.index)

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

    def merge_fastq(self): 
        cmd1 = 'cat {0}/20_fastq/prinseq_{1}/{1}.part_*_1.fastq \
                  > {0}/20_fastq/prinseq_{1}/{1}_1.fastq'.format(self.args.out, 
                                                                 self.index)

        cmd2 = 'cat {0}/20_fastq/prinseq_{1}/{1}.part_*_2.fastq > \
                    {0}/20_fastq/prinseq_{1}/{1}_2.fastq'.format(self.args.out, 
                                                                 self.index)

        cmd3 = 'cat {0}/20_fastq/prinseq_{1}/{1}.part_*_1_singletons.fastq \
                  > {0}/20_fastq/prinseq_{1}/{1}_1_singletons.fastq'.format(self.args.out, 
                                                                            self.index)

        cmd4 = 'cat {0}/20_fastq/prinseq_{1}/{1}.part_*_2_singletons.fastq \
                  > {0}/20_fastq/prinseq_{1}/{1}_2_singletons.fastq'.format(self.args.out, 
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

    def rm_temp_fastq(self):
        for n in range(1, self.N_threads + 1):
            os.remove('{0}/20_fastq/prinseq_{1}/{1}.part_{2:0>3}_1.fastq'.format(self.args.out, self.index, n))
            os.remove('{0}/20_fastq/prinseq_{1}/{1}.part_{2:0>3}_2.fastq'.format(self.args.out, self.index, n))
            os.remove('{0}/20_fastq/prinseq_{1}/{1}.part_{2:0>3}_1_singletons.fastq'.format(self.args.out, self.index, n))
            os.remove('{0}/20_fastq/prinseq_{1}/{1}.part_{2:0>3}_2_singletons.fastq'.format(self.args.out, self.index, n))

        shutil.rmtree('{0}/20_fastq/FaQCs_{1}/{1}.1.trimmed.fastq.split/'.format(self.args.out, self.index))

    def run(self):
        self.seqkit_split2()
        self.run_prinseq()
        self.merge_fastq()
        self.rm_temp_fastq()
