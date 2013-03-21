#ifndef _fasta_h
#define _fasta_h
#include "debug.h"
#include <faidx.h>

#define LOG_4 1.386294 
// Fasta file object
typedef struct {
    PyObject_HEAD
    uint32_t position;
    PyObject *contig;
    PyObject *contigs;
    char *sequence;
    uint16_t counts[4];
    uint32_t length;
    uint32_t l, r, old_position;
    double entropy, gc;
    PyObject *return_value;
    faidx_t *fd;
} tactmod_FastaObject;

// Fasta file object functions
static PyObject *Fasta_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int Fasta_init(tactmod_FastaObject *self, PyObject *args, PyObject *kwds);
PyObject *Fasta_jump(tactmod_FastaObject *self, PyObject *args);
PyObject *Fasta_slice(tactmod_FastaObject *self, PyObject *args);
PyObject *Fasta_load(tactmod_FastaObject *self, PyObject *args);
PyObject *Fasta_tuple(tactmod_FastaObject *self, PyObject *args);
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
    uint32_t position;
    uint32_t length;
    uint16_t gc;
    double entropy;
    char *sequence;
    PyObject *base;
    tactmod_FastaObject *fasta;
} tactmod_FastaIter;

PyObject *tactmod_FastaIter_iter(PyObject *self);
PyObject *tactmod_FastaIter_next(PyObject *self);

PyTypeObject tactmod_FastaIterType;

uint16_t h_f(char *sequence, uint32_t position);
uint16_t h_b(char *sequence, uint32_t position);
double gc_window(char *sequence, uint16_t counts[4]);
double entropy_window(char *sequence, uint16_t counts[4]);
#endif
