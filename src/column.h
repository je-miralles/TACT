#ifndef _column_h
#define _column_h
#include <sam.h>
// The Column class 

#define WT  0
#define HET 1
#define HOM 2

#define b_A 0
#define b_C 1
#define b_G 2
#define b_T 3

typedef struct {
    PyObject_HEAD
    // Every element belonging to a pilup object is a list of length
    // equal to the depth of the pileup column
    long int position;
    unsigned int depth;
    struct {
        unsigned int A;
        unsigned int C;
        unsigned int T;
        unsigned int G;
    } base_counts;
    // By default, a pileup column is interpreted as a column
    // -- there is no relation with this type and an read alignment 
    PyListObject *bases; // Column of readbases
} tactmod_ColumnObject;

PyObject *Column_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

static PyObject *Column_entropy(tactmod_ColumnObject *self, PyObject *args);
static PyObject *Column_genotype(tactmod_ColumnObject *self, PyObject *args);

int Column_setfilter(tactmod_ColumnObject *self, PyObject *args);
int Column_init(tactmod_ColumnObject *self, PyObject *args, PyObject *kwds);
void Column_dealloc(tactmod_ColumnObject *self);

PyMethodDef Column_methods[];
PyMemberDef Column_members[];
PyTypeObject tactmod_ColumnType;


#endif
