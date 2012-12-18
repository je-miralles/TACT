#ifndef _multiseq_h
#define _multiseq_h
#include "sam.h"
#include "tact_bam.h"
#include "debug.h"
#include "base.h"
#include "column.h"

typedef struct {
    PyObject_HEAD
    int position;
    uint8_t n_genomes;
    PyObject *contig;
    PyObject *genomes;
} tactmod_MultiSeqObject;

typedef struct {
    PyObject_HEAD
    int position;
    char *contig;
    int start;
    int end;
    int i;
    PyObject *content;
    tactmod_MultiSeqObject *parent;
} tactmod_GenericIterator;

typedef struct {
    PyObject_HEAD
    int position;
    char *contig;
    int start;
    int end;
    int i;
    PyObject *content;
    tactmod_MultiSeqObject *parent;
    PyObject *iterators;
} tactmod_MultiSeqIterObject;

/* MultiSequence traverser object functions */
static PyObject *MultiSeq_new(PyTypeObject *type, PyObject *args,
                                                  PyObject *kwds);

int MultiSeq_init(tactmod_MultiSeqObject *self, PyObject *args,
                                                PyObject *kwds);

PyObject *MultiSeq_jump(tactmod_MultiSeqObject *self, PyObject *args);
PyObject *MultiSeq_genomes(tactmod_MultiSeqObject *self, PyObject *args);
PyObject *MultiSeq_iterate(tactmod_MultiSeqObject *self, PyObject *args);

void MultiSeq_dealloc(tactmod_MultiSeqObject *self);
void MultiSeqIter_dealloc(tactmod_MultiSeqIterObject *self);

PyObject *MultiSeqIter_iter(PyObject *self);
PyObject *MultiSeqIter_next(PyObject *self);

/* MultiSequence traverser object implementation details */
PyMethodDef MultiSeq_methods[];
PyMemberDef MultiSeq_members[];
PyTypeObject tactmod_MultiSeqType;
PyTypeObject tactmod_MultiSeqIterType;

#endif
