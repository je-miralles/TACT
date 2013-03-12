#ifndef _multiseq_h
#define _multiseq_h
#include "sam.h"
#include "tact_bam.h"
#include "debug.h"
#include "base.h"
#include "column.h"

typedef struct {
    uint32_t min;
    uint32_t max;
    uint8_t size; // if you have more than 256 genomes do something else
    PyTupleObject *iterators;
} priority_heap;

typedef struct {
    PyObject_HEAD
    uint32_t position;
    priority_heap *heap;
    PyObject *contig;
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

void hpush(priority_heap *heap, PyObject *insertion);
PyObject *hpop(priority_heap *heap);
#endif
