#include <Python.h>
#include "structmember.h"
#include "debug.h"
#include "base.h"
/*
 * The Base nucleotide object
 */

PyMethodDef Base_methods[] = {
//    {"__eq__", (PyCFunction)Base__eq__, METH_VARARGS, "Compare bases"},
    {NULL}
};

PyMemberDef Base_members[] = {
    {"quality", T_INT, offsetof(tactmod_BaseObject, quality), 0, "base quality"},
    {NULL}
};

tactmod_BaseObject *chartobase(char base) {
    PyObject *arglist = Py_BuildValue("(c)", base);
    tactmod_BaseObject *ret = NULL;
    ret = PyObject_CallObject((PyObject *)&tactmod_BaseType, arglist);
    //Py_DECREF(arglist);
    return ret;
}

tactmod_BaseObject *inttobase(int base) {
    char _b = 0x00;
    switch(base) {
        case BASE_A:
            _b = 'A';
            break;
        case BASE_C:
            _b = 'C';
            break;
        case BASE_G:
            _b = 'G';
            break;
        case BASE_T:
            _b = 'T';
            break;
        case BASE_N:
            _b = 'N';
            break;
    }
    PyObject *arglist = Py_BuildValue("(c)", _b);
    tactmod_BaseObject *ret = NULL;
    ret = PyObject_CallObject((PyObject *)&tactmod_BaseType, arglist);
    return ret;
}

PyTypeObject tactmod_BaseType = {
    PyObject_HEAD_INIT(NULL)
    0,"tactmod.Base", sizeof(tactmod_BaseObject),
    0,
    (destructor)Base_dealloc,
    (printfunc)Base_print,
    0,
    0,//(cmpfunc)Base_cmp,
    (reprfunc)Base_str,
    0,0,0,0,0,0,
    (reprfunc)Base_str, 
    0,0,0,
    Py_TPFLAGS_DEFAULT,
    "Base nucleotide object",
    0,0,
    (richcmpfunc)Base_cmp, 
    0,0,0,
    Base_methods,
    Base_members,
    0,0,0,0,0,0,
    (initproc)Base_init,
    0,
    Base_new
};


PyObject *Base_str(tactmod_BaseObject *self) {
    return PyString_FromFormat("%c", iupactochar(self->nucleotides));
}

PyObject *Base_print(tactmod_BaseObject *self, FILE *fp, int flags) {
    fprintf(fp, "%c", iupactochar(self->nucleotides));
    return 0;
}
void Base_dealloc(tactmod_BaseObject *self) {
    self->ob_type->tp_free((PyObject*)self);
}

PyObject *Base_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    tactmod_BaseObject *self;

    self = (tactmod_BaseObject *)type->tp_alloc(type, 0);

    if (self != NULL) {
        self->nucleotides = 0;
    }
    return (PyObject *)self;
}

PyObject *Base_cmp(PyObject *self, PyObject *other, int op) {
    PyObject *result = NULL;
    tactmod_BaseObject *_self = (tactmod_BaseObject *) self;
    tactmod_BaseObject *_other = (tactmod_BaseObject *) other;
    switch(op) {
        case Py_LT:
            result = Py_False;
            break;
        case Py_LE:
            result = Py_False;
            break;
        case Py_EQ:
            result = (_self->nucleotides & _other->nucleotides) ? Py_True : Py_False;
            break;
        case Py_NE:
            result = (_self->nucleotides & _other->nucleotides) ? Py_True : Py_False;
            break;
    }
    Py_XINCREF(result);
    return result;  

}

int Base_init(tactmod_BaseObject *self, PyObject *args, PyObject *kwds) {
    char iupac = '-';
    if(!PyArg_ParseTuple(args, "c", &iupac)) {
        return NULL;
    }
    self->nucleotides = chartoiupac(iupac);
    if (0) {

    } else {
        // Negative values are unknown
        self->quality = -1;
        self->tail_distance = -1;
        self->forward_strand = -1;
        self->duplicate = -1;
        self->proper_pair = -1;
        self->indels = 0;
    }
    return 0;
}

iupac_base chartoiupac(char code) {
    iupac_base r = 0;
    switch(code) {
        case 'n': case 'N':
            r |= BASE_N;
            break;

        case 'a': case 'A':
            r |= BASE_A;
            break;

        case 'C': case 'c':
            r |= BASE_C;
            break;

        case 'g': case 'G':
            r |= BASE_G;
            break;

        case 't': case 'T':
            r |= BASE_T;
            break;
    }
    return r;
}

char iupactochar(iupac_base base) {
    if (base == BASE_N) {
        return 'N';
    }
    if (base & BASE_A) {
        return 'A';
    }
    if (base & BASE_C) {
        return 'C';
    }
    if (base & BASE_G) {
        return 'G';
    }
    if (base & BASE_T) {
        return 'T';
    }
    return '-';
}
