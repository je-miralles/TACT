#ifndef _base_h
#define _base_h

/* 
 * The Base Type represents a single nucleotide in a sequence
 * elaborating this type is the ReadBase, which is a Base specifically
 * read from some sequencing or alignment technology and carrying
 * associated information
 */
#define PY_ARRAY_UNIQUE_SYMBOL tctm
typedef uint8_t base4_t;
typedef uint8_t base2_t;

typedef struct {
    PyObject_HEAD
    base4_t nucleotides;
} tactmod_BaseObject;

tactmod_BaseObject *tact_Del;
tactmod_BaseObject *tact_A;
tactmod_BaseObject *tact_C;
tactmod_BaseObject *tact_G;
tactmod_BaseObject *tact_T;
tactmod_BaseObject *tact_N;
tactmod_BaseObject *tact_BaseArray[4];

PyObject *Base_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int Base_init(tactmod_BaseObject *self, PyObject *args, PyObject *kwds);
void Base_dealloc(tactmod_BaseObject *self);

PyObject *Base__eq__(tactmod_BaseObject *self, PyObject *other);
PyObject *Base_cmp(PyObject *self, PyObject *other, int op);

tactmod_BaseObject *chartobase(char b);
char base_char(tactmod_BaseObject *b);
char int_char(base4_t b);

PyObject *Base_str(tactmod_BaseObject *self);
PyObject *Base_print(tactmod_BaseObject *self, FILE *fp, int flags);

uint8_t base4_base2(uint8_t base);

PyMethodDef Base_methods[];
PyMemberDef Base_members[];
PyTypeObject tactmod_BaseType;

#endif
