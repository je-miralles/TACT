#include <Python.h>
#include "structmember.h"
#include "tact_vcf.h"

/* VCF file iterator functions */
PyMethodDef Vcf_methods[] = {
    {"__enter__", (PyCFunction)Vcf_enter, METH_VARARGS, "context entry"},
    {"__exit__", (PyCFunction)Vcf_exit, METH_VARARGS, "context exit"},
    {NULL}
};

PyMemberDef Vcf_members[] = {
    {NULL}
};

PyTypeObject tactmod_VcfType = {
    PyObject_HEAD_INIT(NULL)
    0,"tactmod.Vcf", sizeof(tactmod_VcfObject),
    0,
    (destructor)Vcf_dealloc,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    "VCF file object",
    0,0,0,0,0,0,
    Vcf_methods,
    Vcf_members,
    0,0,0,0,0,0,
    (initproc)Vcf_init,
    0,
    Vcf_new
};

PyTypeObject tactmod_VcfIterType = {
    PyObject_HEAD_INIT(NULL)
    0,"tactmod.VcfIter",sizeof(tactmod_VcfIter),
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER,
    "Vcf file iterator object",
    0,0,0,0,
    tactmod_VcfIter_iter,
    tactmod_VcfIter_next,
};

PyObject *
Vcf_enter(PyObject *self)
{
    return self;
}

PyObject *
Vcf_exit(PyObject *self)
{
    return self;
}

PyObject *
tactmod_VcfIter_iter(PyObject *self)
{
    Py_INCREF(self);
    return self;
}

PyObject *
tactmod_VcfIter_next(PyObject *self)
{
    tactmod_VcfIter *iter = (tactmod_VcfIter *) self;
    return NULL;
}

void
Vcf_dealloc(tactmod_VcfObject *self)
{
    // TODO:
    // close(self->fd);
    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
Vcf_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    tactmod_VcfObject *self;
    self = (tactmod_VcfObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->fd = NULL;
    }
    return (PyObject *)self;
}

int
Vcf_init(tactmod_VcfObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *filename=NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    // TODO:
    // self->fd = openfastq(filename);
    Py_DECREF(filename);
    if (!self->fd) {
        PyErr_SetString(PyExc_IOError, "Cannot open fastq file");
        return NULL;
    }
    return 0;
}

