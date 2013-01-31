#include <Python.h>
#include "sam.h"
#include "structmember.h"
#include "multiseq.h"

/* 
 * multiseq.c defines an object that will iterate of a set of sequences
 */

PyMethodDef MultiSeq_methods[] = {
    {"iterate", (PyCFunction)MultiSeq_iterate, METH_VARARGS,
     "iterate over a range"},
    {"jump", (PyCFunction)MultiSeq_jump, METH_VARARGS, "jump to a position"},
    {"__exit__", (PyCFunction)MultiSeq_exit, METH_VARARGS, "exit context"},
    {"__enter__", (PyCFunction)MultiSeq_exit, METH_VARARGS, "entry context"},
    {NULL}
};

PyMemberDef MultiSeq_members[] = {
    {"position", T_INT, offsetof(tactmod_MultiSeqObject, position), 0,
     "position in reference genome"},
    {NULL}
};

PyTypeObject tactmod_MultiSeqType = {
    PyObject_HEAD_INIT(NULL)
    0,"tactmod.MultiSeq", sizeof(tactmod_MultiSeqObject),
    0,
    (destructor)MultiSeq_dealloc,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    "Multi-sequence traverser object",
    0,0,0,0,0,0,
    MultiSeq_methods,
    MultiSeq_members,
    0,0,0,0,0,0,
    (initproc)MultiSeq_init,
    0,
    MultiSeq_new
};

PyTypeObject tactmod_MultiSeqIterType = {
    PyObject_HEAD_INIT(NULL)
    0,"tactmod.MultiSeqIter", sizeof(tactmod_MultiSeqIter),
    0,
    (destructor)MultiSeqIter_dealloc,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER,
    "Multi-sequence iterator",
    0,0,0,0,
    MultiSeqIter_iter,
    MultiSeqIter_next,
};

PyObject *
MultiSeq_enter(PyObject *self)
{
    return self;
}

PyObject *
MultiSeq_exit(PyObject *self)
{
    return self;
}

PyObject *
MultiSeqIter_iter(PyObject *self)
{
    /* Recursively create iterator objects from
       internal genomic sequences */
    Py_INCREF(self);
    return self;
}

void
MultiSeq_dealloc(tactmod_MultiSeqObject *self)
{
    self->ob_type->tp_free((PyObject*)self);
}

void
MultiSeqIter_dealloc(tactmod_MultiSeqIter *self)
{
    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
MultiSeq_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    tactmod_MultiSeqObject *self;
    self = (tactmod_MultiSeqObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->position = 0;
    }
    return (PyObject *)self;
}

int
MultiSeq_init(tactmod_MultiSeqObject *self, PyObject *args, PyObject *kwds)
{
    PyObject *genomes;
    if (!PyArg_ParseTuple(args, "O", &genomes)) {
        return NULL;
    }
    self->n_genomes = (int)PySequence_Length(genomes);
    self->genomes = PySequence_Fast(genomes, "Error parsing list of genomes");
    return 0;
}

PyObject *
MultiSeq_jump(tactmod_MultiSeqObject *self, PyObject *args)
{
    int start = 0; 
    int end = 0;
    int tid;
    tid = 0;
    char *seq_string;
    PyObject *contig;
    PyObject *callback = NULL;
    //Py_INCREF(args);
    if (!PyArg_ParseTuple(args, "s|O", &seq_string, &callback)) {
        return NULL;
    }
    if (callback != NULL) {
        if (!PyCallable_Check(callback)) {
            PyErr_SetString(PyExc_TypeError, "callback is not callable"); 
            return NULL;
        }
    }
    else {
        callback = Py_None;
    }
    int i;
    PyObject *item;
    PyObject *ret_tuple = PyList_New(self->n_genomes);
    /* TODO: Repair reference counting */
    for (i = 0; i < self->n_genomes; i++) {
        item = PySequence_Fast_GET_ITEM(self->genomes, i);
        Py_INCREF(item);
        PyObject *seq = NULL;
        seq = PyObject_CallMethod(item, "jump", "(s)", seq_string);
        Py_INCREF(seq);
        if (!seq || (seq == Py_None)) {
            Py_INCREF(Py_None);
            return Py_None;
        }
        PyList_SetItem(ret_tuple, i, seq);
    }
    
    Py_DECREF(seq_string);
    end = start + 1;

    if (callback == Py_None) {
        return ret_tuple;
    }
    PyObject *arglist = Py_BuildValue("(O)", ret_tuple);
    PyObject *r;
    Py_INCREF(r);
    r = PyObject_CallObject(callback, arglist);
    return r;
}

PyObject *
MultiSeq_iterate(tactmod_MultiSeqObject *self, PyObject *args)
{
    tactmod_MultiSeqIter *iter;
    char *range = NULL;
    if(!PyArg_ParseTuple(args, "s", &range)) {
        return NULL;
    }
   
    iter = (tactmod_MultiSeqIter *)PyObject_New(
                        tactmod_MultiSeqIter, &tactmod_MultiSeqIterType);

    if (!iter) return NULL;
    if (!PyObject_Init((PyObject *)iter, &tactmod_MultiSeqIterType)) {
        Py_DECREF(iter);
        return NULL;
    }

    iter->parent = self;
    iter->iterators = PyList_New(self->n_genomes);
    int i;
    for(i = 0; i < self->n_genomes; i++) {
        PyObject *item;
        item = PySequence_Fast_GET_ITEM(self->genomes, i);
        if (!item) {
            return NULL;
        }
        tactmod_BamIter *iterator = PyObject_CallMethod(item,
                                                "slice", "(iii)", 0, 100, 200);

        /* Set multiSeq iterator's position to the internal iterator
           if the internal iterator has decided on one */

        if (iterator->position > -1) {
            iter->position = iterator->position;
        }
        PyList_SetItem(iter->iterators, i, iterator);

    }
    iter->start = -1;
    iter->end = -1;
    self->position = -1;
    return (PyObject *)iter;
}

PyObject *
MultiSeqIter_next(PyObject *self)
{
    tactmod_MultiSeqIter *iter = (tactmod_MultiSeqIter *)self;
    int i;
    int start = 0;
    PyObject *ret_tuple = PyList_New(iter->parent->n_genomes);
    for (i = 0; i < iter->parent->n_genomes; i++) {
        tactmod_BamIter *iterator;
        
        iterator = (tactmod_BamIter *)PySequence_Fast_GET_ITEM(
                                                        iter->iterators, i);
        if (start < iterator->start) {
            start = iterator->start;
        }
        int position = iterator->position;
        PyObject *item = PyIter_Next(iterator);
        if(!item) {
            PyErr_SetNone(PyExc_StopIteration);
            return NULL;
        }
        iter->i++;
        PyList_SetItem(ret_tuple, i, item); 
    }
    iter->parent->position = start + iter->i;
    return ret_tuple;
    if (iter->parent->position < iter->start) {
        iter->parent->position++;
        return Py_None;
    
    }
    PyErr_SetNone(PyExc_StopIteration);
    return NULL;
}

