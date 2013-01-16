#ifndef _vcf_h
#define _vcf_h
#include "debug.h"

// VCF file object
typedef struct {
    PyObject_HEAD
    void *fd;
} tactmod_VcfObject;

// VCF file object functions
static PyObject *Vcf_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int Vcf_init(tactmod_VcfObject *self, PyObject *args, PyObject *kwds);
void Vcf_dealloc(tactmod_VcfObject *self);

PyObject *Vcf_enter(PyObject *self);
PyObject *Vcf_exit(PyObject *self);

// VCF file object implementation details
PyMethodDef Vcf_methods[];

PyMemberDef Vcf_members[];

PyTypeObject tactmod_VcfType;

// The iterator
typedef struct {
    PyObject_HEAD
} tactmod_VcfIter;

PyObject *tactmod_VcfIter_iter(PyObject *self);
PyObject *tactmod_VcfIter_next(PyObject *self);

PyTypeObject tactmod_VcfIterType;

#endif
