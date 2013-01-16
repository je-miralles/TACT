from distutils.core import setup, Extension
samtools = ['samtools/kstring.c','samtools/sam_header.c','samtools/razf.c',
            'samtools/bgzf.c','samtools/kaln.c','samtools/kprobaln.c',
            'samtools/bam.c','samtools/faidx.c','samtools/sam.c',
            'samtools/bam_pileup.c','samtools/bam_index.c',
            'samtools/bam_aux.c','samtools/bam_import.c',
            'samtools/bam_md.c','samtools/errmod.c']

tactmod = Extension('tactmod',
        sources = samtools + ['src/base.c','src/column.c','src/fasta.c','src/tact_bam.c',
                   'src/tactmod_init.c', 'src/tact_fastq.c', 'src/tact_vcf.c'],

        library_dirs = [],
        include_dirs = ['samtools/','samtools/bcftools','include/samtools'],
        libraries = ['m', 'z'])

setup (name = 'tactmod',
        version = '0.1',
        description = 'TACT',
        ext_modules = [tactmod])
