#ifndef _base_h
#define _base_h

/* 
 * The Base Type represents a single nucleotide in a sequence
 * elaborating this type is the ReadBase, which is a Base specifically
 * read from some sequencing or alignment technology and carrying
 * associated information
 */

typedef unsigned short int iupac_base;

// Definition of IUPAC nucleotide ambiguity codes
// the bases pair with bits in alphabetical order ascending from the
// lowest bit

#define BASE_A  0x01
#define BASE_C  0x02
#define BASE_G  0x04
#define BASE_T  0x08
#define BASE_N  0x0F

iupac_base chartoiupac(char base);
char iupactochar(iupac_base base);

// The Base type
// TODO: Base should be the parent class of a ReadBase object
typedef struct {
    PyObject_HEAD
    iupac_base nucleotides;
    int quality;
    int tail_distance;
    int forward_strand;
    int proper_pair;
    int duplicate;
    int paired;
    int indels;
} tactmod_BaseObject;

PyObject *Base_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int Base_init(tactmod_BaseObject *self, PyObject *args, PyObject *kwds);
void Base_dealloc(tactmod_BaseObject *self);

PyObject *Base__eq__(tactmod_BaseObject *self, PyObject *other);
PyObject *Base_cmp(PyObject *self, PyObject *other, int op);

tactmod_BaseObject *chartobase(char base);
tactmod_BaseObject *inttobase(int base);

PyObject *Base_str(tactmod_BaseObject *self);
PyObject *Base_print(tactmod_BaseObject *self, FILE *fp, int flags);

PyMethodDef Base_methods[];
PyMemberDef Base_members[];
PyTypeObject tactmod_BaseType;

#endif
