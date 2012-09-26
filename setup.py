from distutils.core import setup, Extension

tactmod = Extension('tactmod',
        sources = ['src/base.c', 'src/column.c', 'src/fasta.c', 'src/bam.c',
                   'src/tactmod_init.c'],

        include_dirs = ['include/samtools'],
        libraries = ['bam'])

setup (name = 'tactmod',
        version = '0.1',
        description = 'TACT',
        ext_modules = [tactmod])
