#ifndef _readbase_h
#define _readbase_h

// The ReadBase type
// TODO: Base should be the parent class of a ReadBase object ?

typedef struct {
    PyObject_HEAD
    tactmod_BaseObject *base;
    uint8_t phred;
    uint16_t tail_distance;
    int16_t indel;
    uint8_t mapq;
} tactmod_ReadBaseObject;


PyObject *ReadBase_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int ReadBase_init(tactmod_ReadBaseObject *self, PyObject *args, PyObject *kwds);
void ReadBase_dealloc(tactmod_ReadBaseObject *self);

PyObject *ReadBase__eq__(tactmod_ReadBaseObject *self, PyObject *other);
PyObject *ReadBase_cmp(PyObject *self, PyObject *other, int op);

PyObject *ReadBase_str(tactmod_ReadBaseObject *self);
PyObject *ReadBase_print(tactmod_ReadBaseObject *self, FILE *fp, int flags);

PyObject *phred_quality(tactmod_ReadBaseObject *self);
PyMethodDef ReadBase_methods[];
PyMemberDef ReadBase_members[];
PyTypeObject tactmod_ReadBaseType;

#endif
