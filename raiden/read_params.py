import os
import sys
import argparse
from multiprocessing import Pool
import multiprocessing as multi
from raiden.utils import time_stamp
from raiden.__init__ import __version__


class Read_params(object):

    def __init__(self, program_name):
        self.program_name = program_name

    def set_options(self):
        if self.program_name == 'raiden':
            parser = self.raiden_options()
        elif self.program_name == 'jiji':
            parser = self.jiji_options()

        if len(sys.argv) == 1:
            args = parser.parse_args(['-h'])
        else:
            args = parser.parse_args()
        return args

    def raiden_options(self):
        parser = argparse.ArgumentParser(description='RaIDeN version {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('raiden -r <FASTA>\n'
                        '              -a <FASTQ1>,<FASTQ2> \n'
                        '              -w <FASTQ1>,<FASTQ2>\n'
                        '              -o <OUT_DIR> \n'
                        '             [-t <INT>]')

        # set options
        parser.add_argument('-r',
                            '--ref',
                            action='store',
                            required=True,
                            type=str,
                            help='Reference fasta. You may use De novo assembled contigs.',
                            metavar='')

        parser.add_argument('-a',
                            '--rna-seq',
                            action='append',
                            required=True,
                            type=str,
                            help=('Fastqs of RNA-seq for genome annotation on input reference genome.\n'
                                  'RNA-seq must contain the seqeunce reads of a causal gene.\n'
                                  'You can use this optiion multiple times for each sample\n'
                                  'Please separate paired fastqs by comma (e.g. -a fastq1,fastq2).\n'),
                            metavar='')

        parser.add_argument('-w',
                            '--whole-genome',
                            action='append',
                            required=True,
                            type=str,
                            help=('Fastqs of whole-genome sequence to filter causal genes.\n'
                                  'This sample should have an opposite trait to the sample\n'
                                  'used for the refernce genome. You can use this optiion multiple\n'
                                  'times for each sample. Please separate paired fastqs by comma\n'
                                  '(e.g. -w fastq1,fastq2).'),
                            metavar='')

        parser.add_argument('-o',
                            '--out',
                            action='store',
                            required=True,
                            type=str,
                            help=('Output directory. Specified name must not exist.'),
                            metavar='')

        parser.add_argument('-t',
                            '--threads',
                            action='store',
                            default=2,
                            type=int,
                            help=('Number of threads. If you specify the number below one,\n'
                                  'then RaIDeN will use the maximum number of threads. [2]'),
                            metavar='')

        parser.add_argument('-d',
                            '--disable-RNAseq-trim',
                            action='store_true',
                            default=False,
                            help='Disable the trimming of RNA-seq reads by FaQCs.')

        parser.add_argument('-D',
                            '--disable-WGS-trim',
                            action='store_true',
                            default=False,
                            help=('Disable the trimming of whole-genome sequence\n'
                                  'by FaQCs and prinseq-lite.'))

        parser.add_argument('-s',
                            '--strand',
                            action='store',
                            default='None',
                            type=str,
                            help=('Strand library option.\n'
                                  " 'fr' : assume a strand library fr-first strand.\n"
                                  " 'rf' : assume a strand library fr-second strand.\n"
                                  " 'None' : Don't assume a strand library.\n"
                                  "Default is 'None'."),
                            metavar='')

        parser.add_argument('-m',
                            '--minimum-len',
                            action='store',
                            default=200,
                            type=int,
                            help=('Minimum length of the transcripts predicted by\n'
                                  'stringtie. [200]'),
                            metavar='')

        parser.add_argument('-q',
                            '--min-MQ',
                            action='store',
                            default=40,
                            type=int,
                            help='Minimum mapping quality during mpileup. [40]',
                            metavar='')

        parser.add_argument('-Q',
                            '--min-BQ',
                            action='store',
                            default=18,
                            type=int,
                            help='Minimum base quality during mpileup. [18]',
                            metavar='')

        # set version
        parser.add_argument('-v',
                            '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))

        return parser

    def jiji_options(self):
        parser = argparse.ArgumentParser(description='jiji version {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('jiji -a <GFF/GTF>\n'
                        '            -v <VCF>\n'
                        '            -b <BED_DIR>\n'
                        '            -o <OUT_DIR> \n'
                        '             ')

        # set options
        parser.add_argument('-a',
                            '--gff',
                            action='store',
                            required=True,
                            type=str,
                            help='GFF/GTF containing gene annotation.',
                            metavar='')

        parser.add_argument('-v',
                            '--vcf',
                            action='store',
                            required=True,
                            type=str,
                            help=('VCF generated by RaIDeN.'),
                            metavar='')

        parser.add_argument('-b',
                            '--bed',
                            action='store',
                            required=True,
                            type=str,
                            help=('Directory of BED files generated by RaIDeN.'),
                            metavar='')

        parser.add_argument('-o',
                            '--out',
                            action='store',
                            required=True,
                            type=str,
                            help=('Output directory. Specified name must not exist.'),
                            metavar='')

        parser.add_argument('--inconsistent',
                            action='store',
                            default=0,
                            type=int,
                            help='Number of the allowed inconsistent markers for SNPs and indels. [0]',
                            metavar='')

        parser.add_argument('--miss',
                            action='store',
                            default=0,
                            type=int,
                            help='Number of the allowed missing markers for SNPs and indels. [0]',
                            metavar='')

        parser.add_argument('--count-het',
                            action='store_true',
                            default=False,
                            help='Count heterozygotes as causal mutation.')

        parser.add_argument('--coverage',
                            action='store',
                            default=0.5,
                            type=float,
                            help='Skip a gene whose coverage more than this parameter for P/A marker. [0.5]', 
                            metavar='')


        # set version
        parser.add_argument('--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))

        return parser

    def check_args(self, args):
        if os.path.isdir(args.out):
            sys.stderr.write(('  Output directory already exist.\n'
                              '  Please rename output directory or '
                                'remove existing directory\n\n'))
            sys.exit(1)

    def check_max_threads(self, args):
        max_cpu = multi.cpu_count()
        print(time_stamp(),
              'maximum number of threads which you can use is up to {}.'.format(max_cpu),
              flush=True)
        if max_cpu <= args.threads:
            sys.stderr.write(('!!WARNING!! You can use up to {0} threads. '
                              'This program will use {0} threads.\n').format(max_cpu))
            sys.stderr.flush()
            args.threads = max_cpu
        elif args.threads < 1:
            args.threads = max_cpu
        return args
