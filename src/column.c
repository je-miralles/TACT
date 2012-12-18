#include <Python.h>
#include "structmember.h"
#include <math.h>
#include "debug.h"
#include "column.h"
#include "base.h"

/* Column Object */
PyMethodDef Column_methods[] = {
    {"entropy", (PyCFunction)Column_entropy, METH_VARARGS, "info content"},
    {"genotype", (PyCFunction)Column_genotype, METH_VARARGS, "genotype"},
    {NULL}
};

PyMemberDef Column_members[] = {
    {"depth", T_INT, offsetof(tactmod_ColumnObject, depth), 0, "depth"},
    {"position", T_INT, offsetof(tactmod_ColumnObject, position), 0, "pos"},
    {"bases", T_OBJECT_EX, offsetof(tactmod_ColumnObject, bases), 0, "bases"},
    {NULL}
};

PyTypeObject tactmod_ColumnType = {
    PyObject_HEAD_INIT(NULL)
    0,"tactmod.Column", sizeof(tactmod_ColumnObject),
    0,
    (destructor)Column_dealloc,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    "Column of alignments",
    0,0,0,0,0,0,
    Column_methods,
    Column_members,
    0,0,0,0,0,0,
    (initproc)Column_init,
    0,
    Column_new
};

void
Column_dealloc(tactmod_ColumnObject *self)
{
    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
Column_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    tactmod_ColumnObject *self;
    self = (tactmod_ColumnObject *)type->tp_alloc(type, 0);
    int position = -1;
    if(!PyArg_ParseTuple(args, "i", &position)) {
        return NULL;
    }
 
    if (self != NULL) {
        self->position = position;
        self->bases = PyList_New(0);
        Py_INCREF(self->bases);
        self->base_counts.A = 0;
        self->base_counts.C = 0;
        self->base_counts.T = 0;
        self->base_counts.G = 0;
    }
    return (PyObject *)self;
}

int
Column_init(tactmod_ColumnObject *self, PyObject *args, PyObject *kwds)
{
    return 0;
}

/* 
 * Shannon Entropy
 * - \sum_{}^{bases} pr(base) * log_4 (pr(base))
 */
static PyObject *
Column_entropy(tactmod_ColumnObject *self, PyObject *args)
{
    float entropy = 0; 
    int base_counts[4];

    base_counts[0] = self->base_counts.A;
    base_counts[1] = self->base_counts.C;
    base_counts[2] = self->base_counts.G;
    base_counts[3] = self->base_counts.T;

    float depth = self->depth;
    float pr = 0;
    int i;
    for (i = 0; i < 4; i++) {
        pr = base_counts[i] / depth;
        if (pr != 0) {
            /* Depend on the compiler to precompute logl(4) */
            entropy += pr * (logl(pr) / logl(4));
        }
    }
   
    entropy *= -1;
    return Py_BuildValue("f", entropy); 
}

/* Binomial Genotype Model */
static PyObject *
Column_genotype(tactmod_ColumnObject *self, PyObject *args)
{
    return NULL;
}

static PyObject *
Column_ploidy(tactmod_ColumnObject *self, PyObject *args)
{
    return NULL;
}
