#ifndef _fasta_h
#define _fasta_h
#include "debug.h"
#include <faidx.h>

// Fasta file object
typedef struct {
    PyObject_HEAD
    long int position;
    PyObject *contig;
    PyObject *contigs;
    faidx_t *fd;
} tactmod_FastaObject;

// Fasta file object functions
static PyObject *Fasta_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int Fasta_init(tactmod_FastaObject *self, PyObject *args, PyObject *kwds);
PyObject *Fasta_jump(tactmod_FastaObject *self, PyObject *args);
PyObject *Fasta_slice(tactmod_FastaObject *self, PyObject *args);
void Fasta_dealloc(tactmod_FastaObject *self);

PyObject *Fasta_enter(PyObject *self);
PyObject *Fasta_exit(PyObject *self);

// Fasta file object implementation details
PyMethodDef Fasta_methods[];

PyMemberDef Fasta_members[];

PyTypeObject tactmod_FastaType;

// The iterator
typedef struct {
    PyObject_HEAD
    long int position;
    long int length;
    long int i;
    char *sequence;
    tactmod_FastaObject *fasta;
} tactmod_FastaIter;

PyObject *tactmod_FastaIter_iter(PyObject *self);
PyObject *tactmod_FastaIter_next(PyObject *self);

PyTypeObject tactmod_FastaIterType;

#endif
