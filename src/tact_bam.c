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
    {"tuple", (PyCFunction)Bam_tuple, METH_VARARGS,
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
//    pileup_buffer buffer;
    tactmod_BamIter *iterator = (tactmod_BamIter *)self; 
    tactmod_BamObject *bam = iterator->bam;
    PyTupleObject *tuple;
    PyTupleObject *features;
    PyIntObject *value;
    while((iterator->position) < iterator->stop) {
        tuple = PyTuple_New(5);
        Py_INCREF(tuple);
        
        value = PyInt_FromLong(iterator->position);
        Py_INCREF(value);
        PyTuple_SET_ITEM(tuple, 0, value); 

        for (i = 1; i < 5; i++) {
            features = PyTuple_New(4);
            Py_INCREF(features);
            value = PyInt_FromLong(0);
            Py_INCREF(value);
            for (j = 0; j < 4; j++) {
                PyTuple_SET_ITEM(features, j, value);
            }
            PyTuple_SET_ITEM(tuple, i, features);
        }

        start = iterator->position + 1;
//        PyObject *arglist = Py_BuildValue("(i)", start);
//        buffer.pileup = PyObject_CallObject((PyObject*)&tactmod_ColumnType,
//                                        arglist);
//        buffer.depth = 0;
//        buffer.buf = bam_plbuf_init(pileup_func, buffer.pileup);
//        iterator->position = buffer.pileup->position;
        iterator->position++;
        end = start + 1;
        tid = 0;
        status = bam_fetch(bam->fd->x.bam, bam->idx, tid,
                       start, end, tuple, fetch_stats);
//        bam_plbuf_push(0, buffer.buf);
//        bam_plbuf_destroy(buffer.buf);
//        if (buffer.pileup->depth == 0) {
//            continue;
//        }
        return tuple;
    }

    PyErr_SetNone(PyExc_StopIteration);
    return NULL;
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
            Py_DECREF(self);
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
Bam_tuple(tactmod_BamObject *self, PyObject *args) {
    return NULL;
}

PyObject *
Bam_jump(tactmod_BamObject *self, PyObject *args) {
    int start, end, tid, status; 

    tid = 0;
    pileup_buffer buffer;
    PyObject *contig;
    buffer.filter == NULL;
    buffer.populate_bases = 1;
    if (!PyArg_ParseTuple(args, "ii|O", &tid, &start, &buffer.filter)) {
        return NULL;
    }
    end = start + 1;

    buffer.pileup = NULL;

    if (buffer.filter != NULL) {
        if (!PyCallable_Check(buffer.filter)) {
            PyErr_SetString(PyExc_TypeError,
                            "filter callback is not callable");
            buffer.filter = Py_None;
        }
    }
    Py_INCREF(buffer.filter);
    if (self->idx == 0) {
        printf("Bam index is not available\n");
        return NULL;
    }
    if (self->fd == 0) {
        printf("Bam could not be loaded\n");
        return NULL;
    }

    if (tid < 0) {
        printf("invalid region\n");
        return NULL;
    }

    PyObject *arglist = Py_BuildValue("(i)", start);
    buffer.pileup = PyObject_CallObject((PyObject*)&tactmod_ColumnType,
                                        arglist);
    buffer.depth = 0;
    buffer.buf = bam_plbuf_init(pileup_func, buffer.pileup);
    status = bam_fetch(self->fd->x.bam, self->idx, tid,
                       start, end, &buffer, fetch_column);
    bam_plbuf_push(0, buffer.buf);
    bam_plbuf_destroy(buffer.buf);
    return buffer.pileup;
}

PyObject *
Bam_slice(tactmod_BamObject *self, PyObject *args) {
    tactmod_BamIter *iter;
    int tid, start, stop;
    pileup_buffer buffer;
    PyObject *contig;
    buffer.filter = NULL;
    buffer.populate_bases = 1;
    if (!PyArg_ParseTuple(args, "iii|O", &tid, &start, &stop, &buffer.filter)) 
    {
        return NULL;
    }
    iter = (tactmod_BamIter *)PyObject_New(tactmod_BamIter, 
                                           &tactmod_BamIterType);
    iter->bam = self;
    iter->offset = 0;
    iter->position = start;
    iter->stop = stop;
    return (PyObject *)iter;
}

PyObject *
Bam_stats(tactmod_BamObject *self, PyObject *args) {
    tactmod_BamIter *iter;
    int tid, start, stop;
    
    if (!PyArg_ParseTuple(args, "iii", &tid, &start, &stop)) return NULL;

    iter = (tactmod_BamIter *)PyObject_New(tactmod_BamIter,
                                           &tactmod_BamIterType);
    Py_INCREF(iter); 
    iter->bam = self;
    iter->offset = 0;
    iter->position = start;
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
    uint8_t offset;
    uint8_t quality, mapping, baq, mapped, reverse, paired, duplicate;
    uint8_t base2;
    long position;
    long old_value;
    PyTupleObject *features;
    PyTupleObject *tuple;
    PyIntObject *value; 
   
    tuple = (PyTupleObject *)data;
    //ops = b->core.n_cigar;
    position = PyInt_AS_LONG(PyTuple_GET_ITEM(tuple, 0));
    offset = position - b->core.pos;
    
    base2 = base4_base2(bam1_seqi(bam1_seq(b), offset));
    if (base2 > 3) {
        return 1;
    }
    features = PyTuple_GET_ITEM(tuple, (base2 + 1));
    
    value = PyTuple_GET_ITEM(features, 0);
    old_value = PyInt_AS_LONG(value);
    value = PyInt_FromLong(old_value++);
    Py_INCREF(value);
    PyTuple_SET_ITEM(features, TOTAL_INDEX, value);

    old_value = PyInt_AS_LONG(PyTuple_GET_ITEM(features, QUALITY_INDEX));
    value = PyInt_FromLong(bam1_qual(b)[offset] + old_value);
    Py_INCREF(value);
    PyTuple_SET_ITEM(features, QUALITY_INDEX, value);
    
    old_value = PyInt_AS_LONG(PyTuple_GET_ITEM(features, DIRECTION_INDEX));
    value = PyInt_FromLong(old_value + b->core.flag & BAM_FREVERSE);
    Py_INCREF(value);
    PyTuple_SET_ITEM(features, DIRECTION_INDEX, value);  
    
    old_value = PyInt_AS_LONG(PyTuple_GET_ITEM(features, MAPPING_INDEX));
    value = PyInt_FromLong(old_value + b->core.qual);
    Py_INCREF(value);
    PyTuple_SET_ITEM(features, MAPPING_INDEX, value);

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

static int
pileup_func(uint32_t tid, uint32_t pos, int n,
            const bam_pileup1_t *pl, void *data) {
    return 0;
}
