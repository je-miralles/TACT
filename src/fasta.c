#include <Python.h>
#include <faidx.h>
#include "structmember.h"
#include "base.h"
#include "fasta.h"

/* Fasta file iterator functions */
PyMethodDef Fasta_methods[] = {
    {"slice", (PyCFunction)Fasta_slice, METH_VARARGS, "slice a range"},
    {"load", (PyCFunction)Fasta_load, METH_VARARGS, "load chromosome"},
    {"tuple", (PyCFunction)Fasta_tuple, METH_VARARGS, "fetch tuple"},
    {"__enter__", (PyCFunction)Fasta_enter, METH_VARARGS, "context entry"},
    {"__exit__", (PyCFunction)Fasta_exit, METH_VARARGS, "context exit"},
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
Fasta_enter(PyObject *self)
{
    return self;
}

PyObject *
Fasta_exit(PyObject *self)
{
    return self;
}

PyObject *
tactmod_FastaIter_iter(PyObject *self)
{
    Py_INCREF(self);
    return self;
}

PyObject *
tactmod_FastaIter_next(PyObject *self)
{
    tactmod_FastaIter *s;

    s = (tactmod_FastaIter *)self;
    if (!s->base) {
        Py_DECREF(s->base);
    }
    tactmod_FastaIter *iter = (tactmod_FastaIter *) self;
    if (iter->position < iter->length) {
//        PyObject *t = char_base(iter->sequence[iter->i]);
        iter->position++;
//        return t;
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
//    self->contigs = Py_
    return (PyObject *)self;
}

int
Fasta_init(tactmod_FastaObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *filename=NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    self->return_value = NULL;
    self->sequence = NULL;
    self->fd = fai_load(filename);
    Py_DECREF(filename);
    if (!self->fd) {
        PyErr_SetString(PyExc_IOError, "Cannot open fasta file");
        return NULL;
    }
    return 0;
}

PyObject *
Fasta_load(tactmod_FastaObject *self, PyObject *args)
{
    uint32_t l;
    char *s;

    if (self->sequence != NULL) {
        free(self->sequence);
    }
    if (!PyArg_ParseTuple(args, "s", &s)) return NULL; 

//    int str_len = snprintf(NULL, 0, "%s:%d-%d", (char *)contig, start, end);
//    const char fetch_str[str_len];
//    sprintf(fetch_str, "%s:%d-%d", (char *)contig, start, end);
    self->sequence = fai_fetch(self->fd, s, &l);
    self->length = l;
    self->counts[0] = 0;
    self->counts[1] = 0;
    self->counts[2] = 0;
    self->counts[4] = 0;
    self->old_position = 0;
    if (self->sequence == NULL) {
        PyErr_SetString(PyExc_ValueError, "Contig does not exist");
        return NULL;
    }
    
    return Py_None;
}

PyObject *
Fasta_tuple(tactmod_FastaObject *self, PyObject *args)
{
    base2_t b;
    PyTupleObject *tuple; 
    if (self->return_value != NULL) {
        Py_DECREF(self->return_value);
    }
    tuple = PyTuple_New(5);
    Py_INCREF(tuple);
    uint32_t i, l, r;
    uint32_t pos;
    if (!PyArg_ParseTuple(args, "i", &pos)) return NULL;
    if (pos > self->length)  {
        return NULL;
    }
    pos--; // lower 1-based index to 0-based
    uint16_t x;
    uint16_t y;
    uint16_t range = 1000;
    if ((self->old_position == 0) || (abs(self->old_position - pos) > 300)) {
        self->counts[0] = 0;
        self->counts[1] = 0;
        self->counts[2] = 0;
        self->counts[3] = 0;
        if (((int)pos - 500) < 0) {
            l = 0;
            r = 1000;
        } else if ((pos + 500) > self->length) {
            l = self->length - 1000;
            r = self->length;
        } else {
            l = (pos - 500);
            r = (pos + 500);
        }
         for (i = l; i < r; i++) {
            b = char_base2(self->sequence[i]);
            if (b <= 3) {
                self->counts[b]++;
            }
        }
        self->entropy = entropy_window(self->sequence, self->counts);
        self->gc = gc_window(self->sequence, self->counts);
        self->old_position = pos;
    }

    x = h_f(self->sequence, pos);
    y = h_b(self->sequence, pos);
    PyTuple_SET_ITEM(tuple, 0, PyInt_FromLong(char_base2(self->sequence[pos])));
    PyTuple_SET_ITEM(tuple, 1, PyInt_FromLong(x));
    PyTuple_SET_ITEM(tuple, 2, PyInt_FromLong(y));
    PyTuple_SET_ITEM(tuple, 3, PyFloat_FromDouble(self->gc));
    PyTuple_SET_ITEM(tuple, 4, PyFloat_FromDouble(self->entropy));
    self->return_value = tuple;
    return tuple;
}

PyObject *
Fasta_slice(tactmod_FastaObject *self, PyObject *args)
{
    long int start;
    long int end;
    char *contig = NULL;
    int tid;
    uint32_t length;
    tactmod_FastaIter *i;
    if (!PyArg_ParseTuple(args, "s", &contig)) return NULL;
    i = PyObject_New(tactmod_FastaIter, &tactmod_FastaIterType);
    if (!i) return NULL;
    Py_INCREF(i);

    if (!PyObject_Init((PyObject *)i, &tactmod_FastaIterType)) {
        Py_DECREF(i);
        return NULL;
    }
    /* TODO: find a way to address positions in the Fasta file without
       this intermediate string */
    
//    length = end - start;
    i->sequence = fai_fetch(self->fd, contig, &length);

    if (!i->sequence) {
        PyErr_SetString(PyExc_ValueError, "Invalid range");
        return NULL;
    }
    i->position = 0;
    i->length = length;
    return (PyObject *)i;
}

uint16_t
h_f(char *sequence, uint32_t position) {
    base2_t b;
    base2_t next;
    uint16_t count;
    count = 0;
    b = char_base2(sequence[position + 1]);
    next = char_base2(sequence[position + 2]);
    while ((b == next) && (b <= 3)) {
        count++;
        position++;
        b = char_base2(sequence[position + 1]);
        next = char_base2(sequence[position + 2]);
    }
    return count;

}

uint16_t
h_b(char *sequence, uint32_t position) {
    base2_t b;
    base2_t next;
    uint16_t count;
    count = 0;
    if (position <= 1) {
        return 0;
    }
    b = char_base2(sequence[position - 1]);
    next = char_base2(sequence[position - 2]);
    while ((b == next) && (position >= 2) && (b <= 3)) {
        count++;
        position--;
        b = char_base2(sequence[position - 1]);
        next = char_base2(sequence[position - 2]);
    }
    return count;
}

double gc_window(char *sequence, uint16_t counts[4]) {
    double gc, total;
    gc = counts[1] + counts[2];
    total = counts[0] + counts[1] + counts[2] + counts[3];
    if (total == 0) {
        return 0;
    }
    return gc / total;
}

double entropy_window(char *sequence, uint16_t counts[4]) {
    double entropy, pr, total;
    uint8_t i;
    entropy = 0;
    total = counts[0] + counts[1] + counts[2] + counts[3];
    if (total == 0) {
        return 0;
    }
    for (i = 0; i < 4; i++) {
        pr = (double)counts[i] / total;
        if (pr != 0) {
            entropy += ((log(pr)/log(4)) * pr);
        }
    }
    return -entropy;
}


