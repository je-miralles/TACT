#include <Python.h>
#include "structmember.h"
#include "debug.h"
#include "base.h"

/* The Base nucleotide object */

PyMethodDef Base_methods[] = {
    {NULL}
};

PyMemberDef Base_members[] = {
    {NULL}
};

tactmod_BaseObject*
char_base(char base)
{
    switch(base) {
        case 'A': case 'a':
            Py_INCREF(tact_A);
            return tact_A;
        case 'C': case 'c':
            Py_INCREF(tact_C);
            return tact_C;
        case 'G': case 'g':
            Py_INCREF(tact_G);
            return tact_G;
        case 'T': case 't':
            Py_INCREF(tact_T);
            return tact_T;
        case 'N': case 'n':
            Py_INCREF(tact_N);
            return tact_N;
        case '.': case '-':
            return (void*)Py_None;
        case 'M': case 'm':
            return (void*)Py_None;
        case 'W': case 'w':
            return (void*)Py_None;
        case 'S':
            return (void*)Py_None;
        case 'Y':
            return (void*)Py_None;
        case 'K':
            return (void*)Py_None;
    }
    return Py_None;
}

base2_t
base4_base2(base4_t base) {
    switch(base) {
        case 0x1:
            return 0x0;
        case 0x2:
            return 0x1;
        case 0x4:
            return 0x2;
        case 0x8:
            return 0x3;
    }
}

base2_t
char_base2(char base) {
    switch(base) {
        case 'a':
        case 'A':
            return 0x0;
        case 'c':
        case 'C':
            return 0x1;
        case 'g':
        case 'G':
            return 0x2;
        case 't':
        case 'T':
            return 0x3;
    }
    return 0x4;
}

PyObject*
complement_b(tactmod_BaseObject* b)
{
    switch(b->nucleotides) {
        case 0x02:
            Py_INCREF(tact_C);
            return tact_C;
        case 0x04:
            Py_INCREF(tact_G);
            return tact_G;
        case 0x0F:
            Py_INCREF(tact_N);
            return tact_N;
    }
    return Py_None;
}

char
base4_char(base4_t b)
{
    switch(b) {
        case 0x01:
            return 'A';
        case 0x02:
            return 'C';
        case 0x04:
            return 'G';
        case 0x08:
            return 'T';
        case 0x0F:
            return 'N';
    }
    return '.';
}

char
base_char(tactmod_BaseObject *b)
{
    return base4_char(b->nucleotides);
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

PyObject *
Base_str(tactmod_BaseObject *self)
{
    return PyString_FromFormat("%c", base4_char(self->nucleotides));
}

PyObject *
Base_print(tactmod_BaseObject *self, FILE *fp, int flags)
{
    fprintf(fp, "%c", base_char(self));
    return 0;
}

void
Base_dealloc(tactmod_BaseObject *self)
{
    self->ob_type->tp_free((PyObject*)self);
}

PyObject *
Base_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    tactmod_BaseObject *self;
    self = (tactmod_BaseObject *)type->tp_alloc(type, 0);

    if (self != NULL) {
        self->nucleotides = 0;
    }
    return (PyObject *)self;
}

PyObject *
Base_cmp(PyObject *self, PyObject *other, int op)
{
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
            if (_self->nucleotides & _other->nucleotides) {
                result = Py_True;
            }
            else {
                result = Py_False;
            }
            break;
        case Py_NE:
            if (_self->nucleotides & _other->nucleotides) {
                result = Py_False;
            }
            else {
                result = Py_True;
            }
            break;
    }
    Py_XINCREF(result);
    return result;  
}

int
Base_init(tactmod_BaseObject *self, PyObject *args, PyObject *kwds)
{
    char iupac;
    self->nucleotides = 0x00;
    return 0;
}
