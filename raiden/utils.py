import os
import sys
from datetime import datetime


def time_stamp():
    return '[RaIDeN:{}]'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

def clean_cmd(cmd):
    return ' '.join(cmd.split())

def call_log(out_dir, name, cmd):
    print(time_stamp(), 
          '!!ERROR!! {}\n'.format(cmd), 
          flush=True)

    print('please check below:\n')

    with open('{}/log/{}.log'.format(out_dir, name)) as log:
        for line in log:
            print(line, end='')

def get_proc_numbers(N_threads, N_files):
    if N_threads < N_files:
        return N_threads, [1]*N_files
    else:
        each_threads = [N_threads//N_files]*N_files
        for i in range(N_threads%N_files):
            each_threads[i] += 1
        return N_files, each_threads

