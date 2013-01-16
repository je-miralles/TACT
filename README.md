TACT
====
Python extension for tandem genome analysis
-------------------------------------------

TACT is a Python extension that exposes sequencing data through a pythonic
object system.

Presently, sequencing data is handled by Samtools.

Examples
--------

If we want to find positions with disagreeing genotypes in a pair of samples:

    import tactmod

    bam1 = tactmod.Bam("sample1.bam")
    bam2 = tactmod.Bam("sample2.bam")
    
    bam1_iter = tactmod.bam1.iterate(lambda x: x.genotype())
    bam2_iter = tactmod.bam2.iterate(lambda x: x.genotype())

    for genotypes in tactmod.multi(bam1_iter, bam2_iter):
        if genotypes[0] is not genotypes[1]:
            print genotypes.position

Note that the position attribute points to the current position; The multi
iterator returns a tuple of appropriate objects with this additional value.

Note also that an arbitrary lambda expression can be passed to an iterator.
For each iteration this callback function will be evaluated.

A simple single nucleotide variant caller:

    import tactmod

    bam = tactmod.Bam("sample.bam")
    ref = tactmod.Fasta("reference.fasta")
    
    for (column, ref_base) in tactmod.multi(bam.iterate(), ref.iterate()):
        print column.base_count(ref_base) / column.depth

Glossary
--------

Included here are some terms as they are relevant to the API, generally
corrosponding to tactmod types.

### Alignment
A Read that is mapped to a reference sequence.  Inherits from the tactmod 
Sequence class.  Includes information about the quality of the mapping, strand
direction, and cigar operations and describe the alignment.

### Bam
A compressed binary format for storing aligned reads.  A tactmod Bam object
serves as an interface for this file type and implements slicing, jumping,
and iteration.

### Base
A member of the set A, C, G, T.  The underlying struct that represents these
nucleotides encodes the identity with 4 bits.  In this way the IUPAC ambiguity
codes can also be represented and tested.

### Column
An unordered sequence of read bases covering a position.  Inherits from the 
tactmod Sequence class.

### Contig
A contiguous sequence of which a group typically constitute a genome.

### Fasta
The fasta format is used to store a set of reference sequences (contigs).
This class implements an iterator over buffered ranges of a fasta file.

### Fastq
A Fastq is an unordered set of Reads from a sequencer.

### Pileup
For a position in a mapped sequence there will be an unordered set of aligned
reads that pile over it.

### Read
An ordered set of bases that have been read using sequencing technology.

### ReadBase
Holds a pointer to a base object.  An instance of a unique base that has
been read on a sequencing platform.

### Reference
A set of reference sequences against which a set of reads is aligned. The
position that a read aligns against a reference is refered to as its
canonical position.

### Sequence
Inherits from the python sequence class.  Adds an integer value to the
object structure that holds the canonical position of the sequence with
respect to a sequence.  Generalizes the relationship between Columns, Reads, and
Alignments.

### VCF
The Variant Call Format is an annotated list of positions.  This object
implements an iterator over this information that can be run in tandem
with another iterator.
