#ifndef _fastq_h
#define _fastq_h
#include "debug.h"

// Fastq file object
typedef struct {
    PyObject_HEAD
    void *fd;
} tactmod_FastqObject;

// Fastq file object functions
static PyObject *Fastq_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int Fastq_init(tactmod_FastqObject *self, PyObject *args, PyObject *kwds);
void Fastq_dealloc(tactmod_FastqObject *self);

PyObject *Fastq_enter(PyObject *self);
PyObject *Fastq_exit(PyObject *self);

// Fastq file object implementation details
PyMethodDef Fastq_methods[];

PyMemberDef Fastq_members[];

PyTypeObject tactmod_FastqType;

// The iterator
typedef struct {
    PyObject_HEAD
} tactmod_FastqIter;

PyObject *tactmod_FastqIter_iter(PyObject *self);
PyObject *tactmod_FastqIter_next(PyObject *self);

PyTypeObject tactmod_FastqIterType;

#endif
