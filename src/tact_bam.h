#ifndef _tact_bam_h
#define _tact_bam_h
#include <sam.h>
#include "debug.h"
#include "base.h"
#include "column.h"

typedef struct {
    int depth;    
    int populate_bases;
    bam_plbuf_t *buf;
    tactmod_ColumnObject *pileup; 
    PyObject *filter;
} pileup_buffer;

typedef struct {
    PyObject_HEAD
    long int position;
    PyObject *contig;
    samfile_t *fd;
    bam_index_t *idx;
    PyObject *column_callback;
} tactmod_BamObject;

// Bamfile object functions
static PyObject *Bam_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int Bam_init(tactmod_BamObject *self, PyObject *args, PyObject *kwds);
PyObject *Bam_jump(tactmod_BamObject *self, PyObject *args);
PyObject *Bam_slice(tactmod_BamObject *self, PyObject *args);
PyObject *Bam_counts(tactmod_BamObject *self, PyObject *args);
PyObject *Bam_pileup(tactmod_BamObject *self, PyObject *args);
void Bam_dealloc(tactmod_BamObject *self);

// Bamfile object implementation details
PyMethodDef Bam_methods[];
PyMemberDef Bam_members[];
PyTypeObject tactmod_BamType;

// The iterator
typedef struct {
    PyObject_HEAD
    tactmod_BamObject *bam;
    long int position;
    long int offset;
    long int start;
    long int stop;
} tactmod_BamIter;

PyObject *tactmod_BamIter_iter(PyObject *self);
PyObject *tactmod_BamIter_next(PyObject *self);

PyTypeObject tactmod_BamIterType;
//int pileup_func(uint32_t tid, uint32_t pos, int n,
//                      const bam_pileup1_t *pl, void *data);

static int fetch_column(const bam1_t *b, void *data);
static int fetch_pileup(const bam1_t *b, void *data);
static int fetch_counts(const bam1_t *b, void *data);

static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data);
#endif
