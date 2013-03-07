#include <Python.h>
#include <sam.h>
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
    return self;
}

PyObject *
tactmod_BamIter_next(PyObject *self) {
    int start;
    int end;
    int status;
    int tid;
    int i, j;
    int loop = 1;
//    pileup_buffer buffer;
    tactmod_BamIter *iterator = (tactmod_BamIter *)self; 
    tactmod_BamObject *bam = iterator->bam;
    bam_plbuf_t *pileup;
    PyTupleObject *tuple;
    PyTupleObject *features;
    PyIntObject *value;
    queue *buffer = iterator->buffer;

    //trace("--\t(%d)\tsize:\t%d\t%d", iterator->position, buffer->size, buffer->end);
    // fall of the end
    if (iterator->position >= iterator->stop) {
        queue_destroy(buffer);
        PyErr_SetNone(PyExc_StopIteration);
        return NULL;
    }

    while ((iterator->position >= buffer->end) || (buffer->size == 0)) {
        if (buffer->end > 0) {
            start = buffer->end;
        }
        else {
            start = iterator->start;
        }

        queue_destroy(buffer);
        buffer = queue_init();

        pileup = bam_plbuf_init(pileup_func, iterator);
        iterator->pileup = pileup;
        bam_fetch(bam->fd->x.bam, bam->idx, 0,
                  start, start + (100000 * loop), (void *)iterator, fetch_f);

        // top off the buffer (as per samtools doc)
        bam_plbuf_push(NULL, pileup);
        bam_plbuf_destroy(pileup);
        loop++;
    }

    iterator->buffer = buffer; 
    tuple = dequeue(buffer);
//    Py_CLEAR(tuple);
//    tuple = Py_None;
    Py_INCREF(tuple);
    iterator->position = buffer->position;
    return tuple;

}

void
Bam_dealloc(tactmod_BamObject *self) {
    bam_index_destroy(self->idx);
    samclose(self->fd);
    //Py_XDECREF(self->text);
    Py_XDECREF(self->contig);
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
    Py_INCREF(iter); 
    iter->bam = self;
    iter->offset = 0;
    iter->position = start;
    iter->start = start;
    iter->stop = stop;
    return (PyObject *)iter;
}

static int
fetch_column(const bam1_t *b, void *data) {
    int i;
    tactmod_BaseObject *read_base;
    PyIntObject *count;
    pileup_buffer *d = (pileup_buffer*)data;
    uint8_t offset = d->pileup->position - b->core.pos;
    uint8_t *p;
    uint32_t ops;
    uint32_t op;
    uint32_t matches;

    //d->pileup->position = b->core.pos;
    //d->position = b->core.pos
    if (1) {
    //if (d->populate_bases) {
        ops = b->core.n_cigar;
        op = bam1_cigar(b)[1] & 0xF;
        matches = 0;

        // THIS DOESNT WORK
        for(i = 0; i < ops; i++) {
            op = bam1_cigar(b)[i] & 0xF;
            if (op == BAM_CMATCH) {
                matches = bam1_cigar(b)[0] >> 4;
            }
            if (op == BAM_CINS) {
                if (offset > matches) {
                    offset += 1;
                }
            }
            if (op == BAM_CDEL) {
                if (offset > matches) {
                    offset -= 1;
                }
            }
        }
        read_base = char_base(base4_char(bam1_seqi(bam1_seq(b),offset)));
        //tactmod_ReadBaseObject *read_base = int2base(bam1_seqi(bam1_s
        p = bam1_qual(b);
        //read_base->phred = p[offset];    
        d->pileup->depth += 1;
        PyList_Append(d->pileup->bases, (PyObject*)read_base);
        if (read_base == (PyObject*)tact_A) {
            d->pileup->base_counts[0]++;
        } else if (read_base == (PyObject*)tact_C) {
            d->pileup->base_counts[1]++;
        } else if (read_base == (PyObject*)tact_G) {
            d->pileup->base_counts[2]++;
        } else if (read_base == (PyObject*)tact_T) {
            d->pileup->base_counts[3]++;
        }
    }
    bam_plbuf_push(b, d->buf);
    d->pileup->counts = PyTuple_New(4);
    Py_INCREF(d->pileup->counts);
    for(i = 0; i < 4; i++) {
        count = PyInt_FromLong(d->pileup->base_counts[i]);
        Py_INCREF(count);
        PyTuple_SET_ITEM(d->pileup->counts, i, count);
    }
    return 0;
}

static int
fetch_stats(const bam1_t *b, void *data) {
    int i;
//    uint8_t offset = d->pileup->position - b->core.pos;
    uint8_t offset, distance, length;
    uint8_t quality, mapping, baq, mapped, reverse, paired, duplicate;
    uint8_t base2;
    long position;
    long old_value;
    PyTupleObject *features;
    PyTupleObject *tuple;
    PyIntObject *value; 
    length = 100; 
    tuple = (PyTupleObject *)data;
    //ops = b->core.n_cigar;
    position = PyInt_AS_LONG(PyTuple_GET_ITEM(tuple, POSITION));
    if ((position - b->core.pos) < 0) {
        return 1;
    }
    offset = position - b->core.pos;
    
    base2 = base4_base2(bam1_seqi(bam1_seq(b), offset));

    //trace("base: %d", bam1_seqi(bam1_seq(b), offset));
    if (base2 > 3) {
        return 1;
    }
    features = PyTuple_GET_ITEM(tuple, base2);
    
//    value = PyTuple_GET_ITEM(features, 0);
    trace("0 indexed base of seq:\t%d", bam1_seqi(bam1_seq(b), 0)); 

    reverse = 1 && (b->core.flag & BAM_FREVERSE);
    old_value = PyInt_AS_LONG(PyTuple_GET_ITEM(features, TOTAL_INDEX));
    value = PyInt_FromLong(old_value + 1);
    Py_INCREF(value);
    PyTuple_SET_ITEM(features, TOTAL_INDEX, value);
    
    old_value = PyInt_AS_LONG(PyTuple_GET_ITEM(features, QUALITY_INDEX));
    value = PyInt_FromLong(bam1_qual(b)[offset] + old_value);
    Py_INCREF(value);
    PyTuple_SET_ITEM(features, QUALITY_INDEX, value);
    
    old_value = PyInt_AS_LONG(PyTuple_GET_ITEM(features, DIRECTION_INDEX));
    value = PyInt_FromLong(old_value + reverse);
    Py_INCREF(value);
    PyTuple_SET_ITEM(features, DIRECTION_INDEX, value);  
    
    old_value = PyInt_AS_LONG(PyTuple_GET_ITEM(features, MAPPING_INDEX));
    value = PyInt_FromLong(old_value + b->core.qual);
    Py_INCREF(value);
    PyTuple_SET_ITEM(features, MAPPING_INDEX, value);

    old_value = PyInt_AS_LONG(PyTuple_GET_ITEM(features, TAIL_DISTANCE));
    if (reverse) {
        distance = length - offset;
    } else {
        distance = offset;
    }
    value = PyInt_FromLong(distance + old_value);
    Py_INCREF(value);
    PyTuple_SET_ITEM(features, TAIL_DISTANCE, value);
    return 0;
}

static int
fetch_pileup(const bam1_t *b, void *data) {
    //tactmod_ReadObject *read;
    // create a Read object from the string of bam1_seq() bases
    //read = PyObject_CallObject(tactmod_Readbam1_seq(b));
    //pileup_buffer *d = (pileup_buffer*)data;
    //bam_plbuf_push(b, d->buf);
    return 0;
}

fetch_f(const bam1_t *b, void *data) {
    // add this alignment to the pileup buffer
    // here would go a low lever filter
    tactmod_BamIter *s = (tactmod_BamIter *)data;
    bam_plbuf_t *pileup = (bam_plbuf_t *)s->pileup;
    bam_plbuf_push(b, pileup);
    return 0;
}

static int
pileup_func(uint32_t tid, uint32_t pos, int n,
            const bam_pileup1_t *pl, void *data) {
    
    // Create the tuple to be returned for every position in the range
    tactmod_BamIter *iterator = (tactmod_BamIter *)data;
    PyTupleObject *tuple;
    PyTupleObject *features;
    PyIntObject *value;
    int start, end;
    queue *buffer = iterator->buffer;
    int i, j;
    if ((pos < iterator->start) || (pos > iterator->stop)) {
        return 0;
    }
    tuple = PyTuple_New(6);
//    Py_INCREF(tuple);
        
    value = PyInt_FromLong(pos);
//    Py_INCREF(value);
    if (0) {
        PyTuple_SET_ITEM(tuple, POSITION, value); 
        value = PyInt_FromLong(n);
        PyTuple_SET_ITEM(tuple, 5, value);
        for (i = 0; i < 4; i++) {
            features = PyTuple_New(5);
//            Py_INCREF(features);
            value = PyInt_FromLong(0);
//            Py_INCREF(value);
            for (j = 0; j < 5; j++) {
                PyTuple_SET_ITEM(features, j, value);
            }
            PyTuple_SET_ITEM(tuple, i, features);
        }
    }
    // append tuple to list
    
    enqueue(buffer, Py_None, pos);
    buffer->end = pos;
//    trace("buffer size:\t%d", buffer->size);
//    PyList_Append(buffer, (PyObject *)tuple);
    return 0;
}

// filling the things that are empty
int enqueue(queue *list, PyTupleObject *content, uint32_t pos) {
    queue_node *node = malloc(sizeof(queue_node));
    queue_node *tail;
    
    if (list->size > 0) {
        tail = list->tail;
        tail->next = node;
    } else {
        list->position = pos;
        list->head = node;
    }

    list->tail = node;
    node->content = content;
    node->position = pos;
    node->next = NULL;
    (list->size)++;
    return 0;
}

// emptying the things that are full
PyTupleObject *dequeue(queue *list) {
//    trace("deque");
    PyTupleObject *content;

    queue_node *head = list->head;
    content = head->content;
    list->position = head->position;
    list->size--;
    list->head = head->next;
    
    free(head);
    // don't forget to decref content!!!!!
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
//    trace("destroying linked list");
    PyObject *item;
    while(list->size > 0) {
        item = dequeue(list);
        Py_CLEAR(item);
    }
    free(list);
    return 0;
}
