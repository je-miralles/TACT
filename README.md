TACT
====
Python extension for tandem genome analysis
-------------------------------------------

TACT is a Python extension that exposes sequencing data through an abstract
type system

Glossary
--------

Included here are some terms as they are relevant to the API, generally
corrosponding to tactmod types.

Alignment
    A sequence of bases that is mapped to a reference sequence.  Inherits from
    the tactmod sequence class.  Includes information about the quality of the
    mapping.

Bam
    A compressed binary format for storing aligned reads.  A tactmod Bam object
    serves as an interface for this file type and implements slicing, jumping,
    and iterators.

Base
    A four bit field.  Tactmod initializes built in objects for these bases.

Column
    A sequence of read bases covering a position.  Implements methods of 
    particular interest to coinciding bases.  Inherits from the tactmod 
    Sequence class.

Contig
    Typically, a contiguous sequence beloning to a Reference.

Fasta
    The fasta format is used to store a set of reference sequences.  This 
    class implements an iterator over buffered ranges of a fasta file.

Pileup
    For a position in a reference sequence there will be a number of aligned
    reads that pile over it.

ReadBase
    Holds a pointer to a base object.  An instance of a unique base that has
    been read on a sequencing platform.

Reference
    A set of reference sequences against which a set of reads is aligned.

Sequence
    Inherits from the python sequence class.  Adds an integer value to the
    object structure that holds the canonical position of the sequence with
    respect to a canonical sequence.
