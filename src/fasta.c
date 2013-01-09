#include <Python.h>
#include <faidx.h>
#include "structmember.h"
#include "base.h"
#include "fasta.h"

/* Fasta file iterator functions */
PyMethodDef Fasta_methods[] = {
    {"jump", (PyCFunction)Fasta_jump, METH_VARARGS, "jump to a position"},
    {"slice", (PyCFunction)Fasta_slice, METH_VARARGS, "slice a range"},
    {NULL}
};

PyMemberDef Fasta_members[] = {
    {"contig", T_OBJECT_EX, offsetof(tactmod_FastaObject, contig), 0, "name"},
    {NULL}
};

PyTypeObject tactmod_FastaType = {
    PyObject_HEAD_INIT(NULL)
    0,"tactmod.Fasta", sizeof(tactmod_FastaObject),
    0,
    (destructor)Fasta_dealloc,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    "Fasta file object",
    0,0,0,0,0,0,
    Fasta_methods,
    Fasta_members,
    0,0,0,0,0,0,
    (initproc)Fasta_init,
    0,
    Fasta_new
};

PyTypeObject tactmod_FastaIterType = {
    PyObject_HEAD_INIT(NULL)
    0,"tactmod.FastaIter",sizeof(tactmod_FastaIter),
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER,
    "Fasta file iterator object",
    0,0,0,0,
    tactmod_FastaIter_iter,
    tactmod_FastaIter_next,
};

PyObject *
tactmod_FastaIter_iter(PyObject *self)
{
    Py_INCREF(self);
    return self;
}

PyObject *
tactmod_FastaIter_next(PyObject *self)
{
    tactmod_FastaIter *iter = (tactmod_FastaIter *) self;
    if (iter->i < iter->length) {
        PyObject *t = chartobase(iter->sequence[iter->i]);
        (iter->i)++;
        return t;
    }
    else {
        free(iter->sequence);
        PyErr_SetNone(PyExc_StopIteration);
    }
    return NULL;
}

void
Fasta_dealloc(tactmod_FastaObject *self)
{
    Py_XDECREF(self->contig);
    fai_destroy(self->fd);
    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
Fasta_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    tactmod_FastaObject *self;
    self = (tactmod_FastaObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->contig = PyString_FromString("");
        if (self->contig == NULL) {
            Py_DECREF(self);
            return NULL;
        }
        self->position = 0;
        self->fd = NULL;
    }
    return (PyObject *)self;
}

int
Fasta_init(tactmod_FastaObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *filename=NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    self->fd = fai_load(filename);
    Py_DECREF(filename);
    if (!self->fd) {
        PyErr_SetString(PyExc_IOError, "Cannot open fasta file");
        return NULL;
    return 0;
    }
}

PyObject *
Fasta_jump(tactmod_FastaObject *self, PyObject *args)
{
    int start, l;
    char *s;
    PyObject *contig = NULL;

    if (!PyArg_ParseTuple(args, "si", &contig, &start)) {
        return NULL;
    }
    self->position = start;
    //self->contig = contig
    int end = start;

    int str_len = snprintf(NULL, 0, "%s:%d-%d", (char *)contig, start, end);
    const char fetch_str[str_len];
    sprintf(fetch_str, "%s:%d-%d", (char *)contig, start, end);
    s = fai_fetch(self->fd, fetch_str, &l);
    if (s == NULL) {
        PyErr_SetString(PyExc_ValueError, "Contig does not exist");
        return NULL;
    }

    if (s[0] == '\x00') {
        Py_INCREF(Py_None);
        return Py_None;
    }

    free(s);
    return chartobase(s[0]);
}

PyObject *
Fasta_slice(tactmod_FastaObject *self, PyObject *args)
{
    long int start;
    long int end;
    PyObject *contig = NULL;
    int length;
    tactmod_FastaIter *i;

    if (!PyArg_ParseTuple(args, "sll", &contig, &start, &end)) return NULL;

    i = PyObject_New(tactmod_FastaIter, &tactmod_FastaIterType);
    if (!i) return NULL;

    if (!PyObject_Init((PyObject *)i, &tactmod_FastaIterType)) {
        Py_DECREF(i);
        return NULL;
    }

    /* TODO: find a way to address positions in the Fasta file without
       this intermediate string */

    int str_len = snprintf(NULL, 0, "%s:%d-%d", (char *)contig,
                                                (int)start, (int)end);
    const char fetch_str[str_len];
    sprintf(fetch_str, "%s:%d-%d", (char *)contig, (int)start, (int)end);
    i->sequence = fai_fetch(self->fd, fetch_str, &length);
    if (!i->sequence) {
        PyErr_SetString(PyExc_ValueError, "Invalid range");
        return NULL;
    }
    i->length = length;
    i->i = 0;
    Py_INCREF(i);
    return (PyObject *)i;
}

