#ifndef _column_h
#define _column_h
#include <sam.h>

/*
    A Column is a Sequence of Bases that all pileup at a single position
*/

typedef struct {
    PyObject_HEAD
    /* Every element belonging to a pilup object is a list of length
       equal to the depth of the pileup column */
    long int position;
    unsigned int depth;
    struct {
        unsigned int A;
        unsigned int C;
        unsigned int T;
        unsigned int G;
    } base_counts;
    PyListObject *bases;
} tactmod_ColumnObject;

static PyObject *Column_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

static PyObject *Column_entropy(tactmod_ColumnObject *self, PyObject *args);
static PyObject *Column_genotype(tactmod_ColumnObject *self, PyObject *args);

int Column_setfilter(tactmod_ColumnObject *self, PyObject *args);
int Column_init(tactmod_ColumnObject *self, PyObject *args, PyObject *kwds);
void Column_dealloc(tactmod_ColumnObject *self);

PyMethodDef Column_methods[];
PyMemberDef Column_members[];
PyTypeObject tactmod_ColumnType;

#endif
