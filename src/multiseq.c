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
    tactmod_MultiSeqObject *parent = ((tactmod_MultiSeqIter *)self)->parent; 
    Py_INCREF(self);
    return self;
}

void
MultiSeq_dealloc(tactmod_MultiSeqObject *self)
{
    free(self->q);
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
    s->q = malloc(sizeof(priorityq));

    self->genomes = PySequence_Fast(genomes, "Error reading genome list");
//    self->heap->size = (uint8_t)PySequence_Length(self->genomes);
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
    tactmod_BamObject *bam;
    tactmod_BamIter *bam_iter;
    tactmod_BamIter *_bam;
    tactmod_MultiSeqObject *s = (tactmod_MultiSeqObject *)self;
    tactmod_BamIter *temp;
    uint8_t i;
    uint32_t start, end;
    uint8_t contig;
    if(!PyArg_ParseTuple(args, "iii", &contig, &start, &end)) {
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
    s->q->_min = 0;
    s->q->_max = 0;
    s->q->max = 0;
    s->q->min = 1;

    for(i = 0; i < 2; i++) {
        bam = (tactmod_BamObject *)PyTuple_GET_ITEM(self->genomes, i); 
        s->iterators[i] = tactmod_BamIter_iter(PyObject_CallMethod(bam, "counts", "iii", 0, 8900, 9100));
        Py_INCREF(s->iterators[i]);
        s->iterations[i] = tactmod_BamIter_next(s->iterators[i]);
        qupdate(s->q, s->iterations[i]);
    }
    s->position = s->q->_max;
    iter->start = 8900;
    iter->end = 9100;
    return (PyObject *)iter;
}

PyObject *
MultiSeqIter_next(PyObject *self)
{
    tactmod_BamIter *iterator;
    tactmod_MultiSeqObject *parent;
    tactmod_MultiSeqIter *s;
    PyObject *ret;
    PyObject *n;
    uint8_t i;
    uint8_t min, max; 
    s = (tactmod_MultiSeqIter*)self;
    parent = s->parent;
    ret = PyTuple_New(2);
    if ((parent->iterations[0] == NULL) || (parent->iterations[1] == NULL)) {
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
    }
    while (parent->q->_max < s->end) {
        ret = PyTuple_New(2);
      
        min = parent->q->_min;
        iterator = parent->iterators[min];
        qupdate(parent->q, iterator);

       if (parent->q->_min < parent->q->_max) {
            if (iterator->position > parent->q->_min) {
                parent->iterations[parent->q->min] = tactmod_BamIter_next(iterator);                
                qupdate(parent->q, parent->iterators);
            }
        }

        if (parent->q->_max == parent->q->_min) {
            PyTuple_SET_ITEM(ret, 0, parent->iterations[0]);
            PyTuple_SET_ITEM(ret, 1, parent->iterations[1]);
            parent->position = parent->q->_max;
            // call next on min iterator before returning
            parent->iterations[parent->q->min] = tactmod_BamIter_next(iterator);
            // push it back onto the queue
            qupdate(parent->q, iterator);
            return ret;
        }
        else {
            qupdate(parent->q, iterator);
            // push back on the popped value
        }
    }
    PyErr_SetNone(PyExc_StopIteration);
    return NULL;
}

void qupdate(priorityq *q, tactmod_BamIter *iterators[2]) {
    tactmod_BamIter *iterator = iterators[0];
    tactmod_BamIter *temp;
//    iterator = (tactmod_BamIter *)PyTuple_GET_ITEM(q->iterators, i);

    if (iterator->position <= q->_min) {
        temp = q->min;
        q->min = iterator;
        q->max = temp;
//        q->_min = q->min->position;
//        q->_max = q->max->position;
    } else {
        temp = q->max;
        q->max = iterator;
        q->min = temp;
//        q->_min = q->min->position;
//        q->_max = q->max->position;
    }

   
}

