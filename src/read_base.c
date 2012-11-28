#include <Python.h>
#include "structmember.h"
#include "debug.h"

#include "base.h"
#include "read_base.h"

/*
 * This object represents a Base that has been read with associated quality
 * score
 */

PyMethodDef ReadBase_methods[] = {
//    {"__eq__", (PyCFunction)ReadBase__eq__, METH_VARARGS, "Compare bases"},
    {"quality", (PyCFunction)phred_quality, METH_VARARGS, 
                                "phred quality scaled into [0..1]"},
    {NULL}
};

PyMemberDef ReadBase_members[] = {
    {"phred", T_INT, offsetof(tactmod_ReadBaseObject, phred), 0, 
                                "phred scaled sequencing quality"},
    {"map_quality", T_INT, offsetof(tactmod_ReadBaseObject, mapq), 0,
                                "phred scaled mapping quality"},
    {NULL}
};

PyTypeObject tactmod_ReadBaseType = {
    PyObject_HEAD_INIT(NULL)
    0,"tactmod.ReadBase", sizeof(tactmod_ReadBaseObject),
    0,
    (destructor)ReadBase_dealloc,
    (printfunc)ReadBase_print,
    0,
    0, //(cmpfunc)ReadBase_cmp,
    (reprfunc)ReadBase_str,
    0,0,0,0,0,0,
    (reprfunc)ReadBase_str, 
    0,0,0,
    Py_TPFLAGS_DEFAULT,
    "ReadBase",
    0,0,
    0,// (richcmpfunc)ReadBase_cmp, 
    0,0,0,
    ReadBase_methods,
    ReadBase_members,
    0,0,0,0,0,0,
    (initproc)ReadBase_init,
    0,
    ReadBase_new
};

void ReadBase_dealloc(tactmod_ReadBaseObject *self) {
    self->ob_type->tp_free((PyObject*)self);
}
PyObject *phred_quality(tactmod_ReadBaseObject *self) {
    return Py_None;
}

PyObject *ReadBase_str(tactmod_ReadBaseObject *self) {
    return PyString_FromFormat("%c", basetochar(self->base->nucleotides));
}

PyObject *ReadBase_print(tactmod_ReadBaseObject *self, FILE *fp, int flags) {
    fprintf(fp, "%c", basetochar(self->base));
    return 0;
}

PyObject *ReadBase_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    tactmod_ReadBaseObject *self;
    self = (tactmod_ReadBaseObject *)type->tp_alloc(type, 0);

    if (self != NULL) {
        self->base = tact_Del;
        self->phred = 0;
        self->mapq = 0;
    }
    return (PyObject *)self;
}

int ReadBase_init(tactmod_ReadBaseObject *self, PyObject *args, PyObject *kwds) {
    return 0;
}

