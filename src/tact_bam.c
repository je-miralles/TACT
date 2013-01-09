#include <Python.h>
#include <sam.h>
#include "structmember.h"
#include "tact_bam.h"

/* 
 * tact_bam.c serves as a driver for the samtools bam library
 *
 */

PyMethodDef Bam_methods[] = {
    {"jump", (PyCFunction)Bam_jump, METH_VARARGS, "jump to a position"},
    {"slice", (PyCFunction)Bam_slice, METH_VARARGS, "slice a range"},
    {NULL}
};

PyMemberDef Bam_members[] = {
    {"column_callback", T_OBJECT_EX, offsetof(tactmod_BamObject,
                column_callback), 0, "pileup column callback function"},
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
    return NULL;
}

PyObject *
tactmod_BamIter_next(PyObject *self) {
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
        self->position = 0;
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
                       start, end, &buffer, fetch_func);
    bam_plbuf_push(0, buffer.buf);
    bam_plbuf_destroy(buffer.buf);
    return buffer.pileup;
}

PyObject *
Bam_slice(tactmod_BamObject *self, PyObject *args) {
    return NULL;
}

static int
fetch_func(const bam1_t *b, void *data) {
    tactmod_BaseObject *read_base;
    pileup_buffer *d = (pileup_buffer*)data;
    uint8_t offset = d->pileup->position - b->core.pos;
    uint8_t *p;
    if (d->populate_bases) {
        
        read_base = chartobase(inttochar(bam1_seqi(bam1_seq(b),offset)));
        //tactmod_ReadBaseObject *read_base = int2base(bam1_seqi(bam1_s
        p = bam1_qual(b);
        //read_base->phred = p[offset];    

        if (d->filter == Py_None) {
            d->pileup->depth += 1;
            PyList_Append(d->pileup->bases, (PyObject*)read_base);
        } else {
            PyObject *result;
            PyObject *arglist = Py_BuildValue("(O)", read_base);
            if (!arglist) {
                return 1;
            }
            result = PyObject_CallObject(d->filter, arglist);
            if (!result) {
                return 1;
            }
            
            //result = Py_True;

            Py_INCREF(result);
            if (PyBool_Check(result)) {
            } else {
                return 1;
            }
            if (result == Py_True) {
                d->pileup->depth += 1;
                PyList_Append(d->pileup->bases, (PyObject*)read_base);
                if (read_base->nucleotides == (PyObject*)tact_A) {
                    d->pileup->base_counts.A++;
                } else if (read_base->nucleotides == (PyObject*)tact_C) {
                    d->pileup->base_counts.C++;
                } else if (read_base->nucleotides == (PyObject*)tact_G) {
                    d->pileup->base_counts.G++;
                } else if (read_base->nucleotides == (PyObject*)tact_T) {
                    d->pileup->base_counts.T++;
                }
            }
        }
    }
    bam_plbuf_push(b, d->buf);
    return 0;
}

static int
pileup_func(uint32_t tid, uint32_t pos, int n,
            const bam_pileup1_t *pl, void *data) {
    return 0;
}
