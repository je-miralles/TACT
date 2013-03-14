#include <Python.h>
#include <sam.h>
#include <stdlib.h>
#include "structmember.h"
#include "tact_bam.h"

/* 
 * tact_bam.c serves as a driver for the samtools bam library
 *
 */

PyMethodDef Bam_methods[] = {
//    {"column", (PyCFunction)Bam_slice, METH_VARARGS,
//                "columns covering a range"},
//    {"pileup", (PyCFunction)Bam_pileup, METH_VARARGS,
//                "pileups covering a range"},
    {"counts", (PyCFunction)Bam_counts, METH_VARARGS,
                "fixed size tuples"},
    {NULL}
};

PyMemberDef Bam_members[] = {
    {NULL}
};

PyTypeObject tactmod_BamType = {
    PyObject_HEAD_INIT(NULL)
    0,"tactmod.Bam", sizeof(tactmod_BamObject),
    0,
    (destructor)Bam_dealloc,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    "Bam file object",
    0,0,0,0,0,0,
    Bam_methods,
    Bam_members,
    0,0,0,0,0,0,
    (initproc)Bam_init,
    0,
    Bam_new
};

PyTypeObject tactmod_BamIterType = {
    PyObject_HEAD_INIT(NULL)
    0,"tactmod.BamIter",sizeof(tactmod_BamIter),
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER,
    "Bam file iterator object",
    0,0,0,0,
    tactmod_BamIter_iter,
    tactmod_BamIter_next,
};

PyObject *
tactmod_BamIter_iter(PyObject *self) {
    // Initialize pileup buffer
    tactmod_BamIter *s = (tactmod_BamIter *)self;
    s->buffer = queue_init();
    Py_INCREF(self);
    // Check for the rightmost limit of this contig
    return self;
}

PyObject *
tactmod_BamIter_next(PyObject *self) {
    uint32_t start;
    uint32_t end;
    int status;
    int tid;
    uint8_t i, j;
    int loop = 1;
    uint32_t stop = 0;
//    pileup_buffer buffer;
    tactmod_BamIter *iterator = (tactmod_BamIter *)self; 
    tactmod_BamObject *bam = iterator->bam;
    bam_plbuf_t *pileup;
    PyTupleObject *tuple = iterator->return_value;
    PyTupleObject *features;
    PyTupleObject *nested_tuple;
    PyIntObject *value;
    queue *buffer = iterator->buffer;
    column_t column;
    // fall of the end
    if (iterator->position >= iterator->stop) {
        queue_destroy(buffer);
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
    }

    while (((iterator->position >= buffer->end) || (buffer->size == 0)) &&
           (buffer->end < iterator->stop)) {
      
        if (buffer->position == 0) {
            start = iterator->start;
        } else {
            start = buffer->end;
        }

        stop = start + BUFFER_SIZE;

        if (stop > iterator->stop) {
            stop = iterator->stop;
        }

        queue_destroy(buffer);
        buffer = queue_init();
        buffer->fetch_start = start;
        buffer->fetch_stop = stop;
        trace("buffering %d - %d", start, stop);
        pileup = bam_plbuf_init(pileup_func, iterator);
        iterator->pileup = pileup;
        bam_fetch(bam->fd->x.bam, bam->idx, 0,
                  start, stop, (void *)iterator, fetch_f);
        // top off the buffer (as per samtools doc)
        bam_plbuf_push(NULL, pileup);
        bam_plbuf_destroy(pileup);
        loop++;
    }
//
//  Construct the tuple
//    
    column = dequeue(buffer);
    iterator->buffer = buffer; 
    iterator->position = column.position;
//    if (column.features[0] > 30) {
//        trace("lots of reverse strands");
//    }
//    trace("%d", iterator->position);
    if (tuple) {
        Py_DECREF(tuple);
    }
    tuple = PyTuple_New(11);
    PyTuple_SET_ITEM(tuple, 0, PyInt_FromLong((long)column.position));
    for (i = 0; i < 5; i++) {
        nested_tuple = PyTuple_New(6);
        for (j = 0; j < 6; j++) {
            PyTuple_SET_ITEM(nested_tuple, j, PyInt_FromLong((long)column.features[i][j]));
        }
        PyTuple_SET_ITEM(tuple, i + 1, nested_tuple);
    }
    PyTuple_SET_ITEM(tuple, 6, PyInt_FromLong((long)column.major));
    PyTuple_SET_ITEM(tuple, 7, PyInt_FromLong((long)column.minor));
    PyTuple_SET_ITEM(tuple, 8, PyInt_FromLong((long)column.ambiguous));
    PyTuple_SET_ITEM(tuple, 9, PyInt_FromLong((long)column.indels));
    PyTuple_SET_ITEM(tuple, 10, PyFloat_FromDouble(column.entropy));
//    PyTuple_SET_ITEM(tuple, 11, PyFloat_FromDouble(0.0));
//
//        
//        tuple = Py_None;
    Py_INCREF(tuple);
    iterator->return_value = tuple; 
    return tuple;
}

void
Bam_dealloc(tactmod_BamObject *self) {
    bam_index_destroy(self->idx);
    samclose(self->fd);
    //Py_XDECREF(self->text);
    Py_CLEAR(self->contig);
    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *
Bam_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    tactmod_BamObject *self;
    self = (tactmod_BamObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->contig = PyString_FromString("");
        if (self->contig == NULL) {
            Py_CLEAR(self);
            return NULL;
        }
        self->fd = NULL;
    }
    return (PyObject *)self;
}

int
Bam_init(tactmod_BamObject *self, PyObject *args, PyObject *kwds) {
    PyObject *filename = NULL;
    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    self->fd = samopen(filename, "rb", 0);
    self->idx = bam_index_load(filename);

    if (!self->fd) {
        PyErr_SetString(PyExc_IOError, "Cannot open bam file");
        return NULL;
    return 0;
    }
}


PyObject *
Bam_jump(tactmod_BamObject *self, PyObject *args) {
    return Py_None;
}

PyObject *
Bam_slice(tactmod_BamObject *self, PyObject *args) {
    return Py_None;
}

PyObject *
Bam_counts(tactmod_BamObject *self, PyObject *args) {
    tactmod_BamIter *iter;
    int tid, start, stop;
    
    if (!PyArg_ParseTuple(args, "iii", &tid, &start, &stop)) return NULL;

    iter = (tactmod_BamIter *)PyObject_New(tactmod_BamIter,
                                           &tactmod_BamIterType);
//    Py_INCREF(iter); 
    iter->bam = self;
    iter->offset = 0;
    iter->position = start;
    iter->start = start;
    iter->stop = stop;
    return (PyObject *)iter;
}

fetch_f(const bam1_t *b, void *data) {
    // add this alignment to the pileup buffer
    // here would go a low level filter
    if (b->core.n_cigar > 1) {
        return 0;
    }
    tactmod_BamIter *s = (tactmod_BamIter *)data;
    bam_plbuf_t *pileup = (bam_plbuf_t *)s->pileup;
    bam_plbuf_push(b, pileup);
    return 0;
}

static int
pileup_func(uint32_t tid, uint32_t pos, int n,
            const bam_pileup1_t *pl, void *data) {
    int offset, distance, length;
    uint8_t quality, mapping, baq, mapped, reverse, paired, duplicate;
    uint8_t base2;
    long old_value;
    uint8_t i, j;
    bam1_t *b;
    bam_pileup1_t alignment;
    uint16_t base_counts[4] = {0, 0, 0, 0};
    column_t column;
    for (i = 0; i < 5; i++) {
        for (j = 0; j < 6; j++) {
            column.features[i][j] = 0;
        }
    }
    tactmod_BamIter *iterator = (tactmod_BamIter *)data;
    int start, end;
    queue *buffer = iterator->buffer;
    if ((pos < buffer->fetch_start) || (pos > buffer->fetch_stop)) {
        return 0;
    }

    column.indels = 0;
    column.ambiguous = 0;
        //
    column.depth = n;
    for (i = 0; i < n; i++) {
        // append tuple to list
        length = 100; 
        //ops = b->core.n_cigar;
        alignment = pl[i];
        b = alignment.b;
        column.indels += alignment.indel;

        
        //offset = pos - b->core.pos;
        offset = alignment.qpos; 
        base2 = base4_base2(bam1_seqi(bam1_seq(b), offset));
        if (base2 <= 3)  {
            base_counts[base2]++;
        } else {
            column.ambiguous++;
            continue;
        }
        reverse = 1 && (b->core.flag & BAM_FREVERSE);
        quality = bam1_qual(b)[offset];
        mapping = b->core.qual;
        distance = 0;
        if (reverse) {
            distance = length - offset;
        } else {
            distance = offset;
        }

        column.features[base2][0] += 1;
        column.features[base2][1] += quality;
        column.features[base2][2] += mapping;
        column.features[base2][3] += distance;
        column.features[base2][4] += reverse;
        column.features[base2][5] += 0;

        column.features[4][0] += 1;
        column.features[4][1] += quality;
        column.features[4][2] += mapping;
        column.features[4][3] += distance;
        column.features[4][4] += reverse;
        column.features[4][5] += 0;
    }
    column.major = 0;
    uint16_t max = 0;
    for (i = 0; i < 4; i++) {
        if (base_counts[i] > max) {
            column.major = i;
            max = base_counts[i];
        }
    }
    column.minor = column.major;
    max = 0;
    for (i = 0; i < 4; i++) {
        if ((base_counts[i] > max) && (i != column.major)) {
            column.minor = i;
            max = base_counts[i];
        }
    }
//    trace("major:\t%d\tminor:\t%d", column.major, column.minor);
//    uint16_t t = base_counts[column.major] + base_counts[column.minor];
//    column.binomial_lll = binomial_ll(column.major, column.depth, 0.001); 
//    column.binomial_ll = binomial_ll(column.major, column.depth, 0.5);

    column.position = (pos + 1); // HERE BE DRAGONS
    column.entropy = entropy(base_counts, column.depth);

    enqueue(buffer, column, column.position);
    buffer->end = column.position;
    return 0;
}

// filling the things that are empty
int enqueue(queue *list, column_t content, uint32_t pos) {
    queue_node *node;
    node = (queue_node *)malloc(sizeof(queue_node));
    if (list->size == 0) {
        list->position = pos;
        list->head = node;
    } else {
        list->tail->next = node;
    }
     
    list->tail = node;
    node->content = content;
    node->position = pos;
    node->next = NULL;
    (list->size)++;
    return 0;
}

// emptying the things that are full
column_t dequeue(queue *list) {
    column_t content;

    queue_node *next;

    content = list->head->content;
    list->position = list->head->position;
    list->size--;
    next = list->head->next;

    free(list->head);
    list->head = next;
    return content;
}
queue *queue_init(void) {
    queue *buffer;
    buffer = (queue *)malloc(sizeof(queue));
    buffer->size = 0;
    buffer->end = 0;
    buffer->tail = NULL;
    buffer->head = NULL;
    buffer->position = 0;
    return buffer;
}

int queue_destroy(queue *list) {
    column_t item;
    int count = 0;
    while(list->size > 0) {
        item = dequeue(list);
        count++;
    }
    free(list);
    list = NULL;
    return 0;
}

double binomial_ll(uint16_t k, uint16_t n, double mu) {
    uint16_t i;
    uint32_t d;
    unsigned long r = 1;
    d = 0;
    double x;
    for (i = n; i > k; i--) {
        r *= i;
        while (d > 1 && (r % d)) {
            r /= d--;
        }
    }
    x = log(r) + (k * log(mu)) + (d * log(1 - mu));
    //trace("major: %d,\ttotal: %d, score: %f (%f)", k, n, x, log(-x));
    return log(-x);
}

double entropy(uint16_t bases[4], uint16_t depth) {
    uint8_t i;
    double e = 0;
    float pr;
    for (i = 0; i < 4; i++) {
        pr = (float)bases[i] / (float)depth;
        if (pr != 0) {
            e += ((log(pr)/log(4)) * pr); 
        }
    }
    return -e;
}
