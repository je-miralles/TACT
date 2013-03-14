#ifndef _multiseq_h
#define _multiseq_h
#include "sam.h"
#include "tact_bam.h"
#include "debug.h"
#include "base.h"
#include "column.h"

#define PY_ARRAY_UNIQUE_SYMBOL tctm
typedef struct {
    uint32_t _min;
    uint32_t _max;
    uint8_t max; // priority queue only needs to retain index now
    uint8_t min;
} priorityq;

typedef struct {
    PyObject_HEAD
    uint32_t position;
    priorityq *q;
    PyObject *contig;
    PyTupleObject *genomes;
    tactmod_BamIter *iterators[2];
    PyObject *iterations[2]; // this is weird
} tactmod_MultiSeqObject;

typedef struct {
    PyObject_HEAD
    uint32_t start;
    uint32_t end;
    tactmod_MultiSeqObject *parent;
} tactmod_MultiSeqIter;


/* MultiSequence traverser object functions */
static PyObject *MultiSeq_new(PyTypeObject *type, PyObject *args,
                                                  PyObject *kwds);

int MultiSeq_init(tactmod_MultiSeqObject *self, PyObject *args,
                                                PyObject *kwds);

PyObject *MultiSeq_jump(tactmod_MultiSeqObject *self, PyObject *args);
PyObject *MultiSeq_genomes(tactmod_MultiSeqObject *self, PyObject *args);
PyObject *MultiSeq_iterate(tactmod_MultiSeqObject *self, PyObject *args);

void MultiSeq_dealloc(tactmod_MultiSeqObject *self);
void MultiSeqIter_dealloc(tactmod_MultiSeqIter *self);

PyObject *MultiSeqIter_iter(PyObject *self);
PyObject *MultiSeqIter_next(PyObject *self);

PyObject *MultiSeq_enter(PyObject *self);
PyObject *MultiSeq_exit(PyObject *self);
/* MultiSequence traverser object implementation details */
PyMethodDef MultiSeq_methods[];
PyMemberDef MultiSeq_members[];
PyTypeObject tactmod_MultiSeqType;
PyTypeObject tactmod_MultiSeqIterType;

void qupdate(priorityq *q, tactmod_BamIter *iterators[2]);
#endif
