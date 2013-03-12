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

    tactmod_MultiSeqObject *parent = ((tactmod_MultiSeqIter *)self)->parent; 
    uint8_t i;
    
    // initialize heap
    parent->heap->min = 0;
    parent->heap->max = 0;

    // push iterators onto heap

    for (i = 0; i < parent->heap->size; i++) {
//        PyTuple_SET_ITEM(self->heap->iterators, i, iterator);        
    }
    Py_INCREF(self);
    return self;
}

void
MultiSeq_dealloc(tactmod_MultiSeqObject *self)
{
    free(self->heap);
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
    tactmod_MultiSeqObject *s = (tactmod_MultiSeqObject *)self;

    if (!PyArg_ParseTuple(args, "O", &genomes)) {
        return NULL;
    }
    // initialize heap
    s->heap = malloc(sizeof(priority_heap));
    self->heap->size = (uint8_t)PySequence_Length(genomes);
    self->heap->iterators = PySequence_Fast(genomes, "Error parsing genomes");
    return 0;
}

PyObject *
MultiSeq_jump(tactmod_MultiSeqObject *self, PyObject *args)
{
    return Py_None;
}

PyObject *
MultiSeq_iterate(tactmod_MultiSeqObject *self, PyObject *args)
{
    tactmod_MultiSeqIter *iter;
    char *range = NULL;
    uint8_t i;
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
    for(i = 0; i < self->heap->size; i++) {
        PyObject *item;
        item = PySequence_Fast_GET_ITEM(self->heap->iterators, i);
        if (!item) {
            return NULL;
        }
        tactmod_BamIter *iterator = PyObject_CallMethod(item,
                                                "counts", "(iii)", 0, 100, 200);

        /* Set multiSeq iterator's position to the internal iterator
           if the internal iterator has decided on one */

        if (iterator->position > -1) {
            iter->parent->position = iterator->position;
        }
        PyList_SetItem(iter->parent->heap->iterators, i, iterator);

    }
    iter->start = -1;
    iter->end = -1;
    self->position = -1;
    return (PyObject *)iter;
}

PyObject *
MultiSeqIter_next(PyObject *self)
{
    PyObject *iterator;
    
    PyErr_SetNone(PyExc_StopIteration);
    return NULL;
}

void hpush(priority_heap *heap, PyObject *insertion) {
    
}

PyObject *hpop(priority_heap *heap) {
    return Py_None;
}
