#include <Python.h>
#include "structmember.h"
#include "tact_fastq.h"

/* Fastq file iterator functions */
PyMethodDef Fastq_methods[] = {
    {"__enter__", (PyCFunction)Fastq_enter, METH_VARARGS, "context entry"},
    {"__exit__", (PyCFunction)Fastq_exit, METH_VARARGS, "context exit"},
    {NULL}
};

PyMemberDef Fastq_members[] = {
    {NULL}
};

PyTypeObject tactmod_FastqType = {
    PyObject_HEAD_INIT(NULL)
    0,"tactmod.Fastq", sizeof(tactmod_FastqObject),
    0,
    (destructor)Fastq_dealloc,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    "Fastq file object",
    0,0,0,0,0,0,
    Fastq_methods,
    Fastq_members,
    0,0,0,0,0,0,
    (initproc)Fastq_init,
    0,
    Fastq_new
};

PyTypeObject tactmod_FastqIterType = {
    PyObject_HEAD_INIT(NULL)
    0,"tactmod.FastqIter",sizeof(tactmod_FastqIter),
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER,
    "Fastq file iterator object",
    0,0,0,0,
    tactmod_FastqIter_iter,
    tactmod_FastqIter_next,
};

PyObject *
Fastq_enter(PyObject *self)
{
    return self;
}

PyObject *
Fastq_exit(PyObject *self)
{
    return self;
}

PyObject *
tactmod_FastqIter_iter(PyObject *self)
{
    Py_INCREF(self);
    return self;
}

PyObject *
tactmod_FastqIter_next(PyObject *self)
{
    tactmod_FastqIter *iter = (tactmod_FastqIter *) self;
    return NULL;
}

void
Fastq_dealloc(tactmod_FastqObject *self)
{
    // TODO:
    // close(self->fd);
    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
Fastq_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    tactmod_FastqObject *self;
    self = (tactmod_FastqObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->fd = NULL;
    }
    return (PyObject *)self;
}

int
Fastq_init(tactmod_FastqObject *self, PyObject *args, PyObject *kwds)
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

