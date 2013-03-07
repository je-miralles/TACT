#ifndef _tact_bam_h
#define _tact_bam_h
#include <sam.h>
#include "debug.h"
#include "base.h"
#include "column.h"

#define TOTAL_INDEX 0
#define QUALITY_INDEX   1
#define DIRECTION_INDEX 2
#define MAPPING_INDEX   3
#define TAIL_DISTANCE   4

#define POSITION  4
#define MAX_INDEX   5
#define MINOR_INDEX 6
#define DEPTH 7
#define INDELS  8

typedef struct {
    int depth;    
    int populate_bases;
    bam_plbuf_t *buf;
    tactmod_ColumnObject *pileup; 
    PyObject *filter;
} pileup_buffer;

typedef struct {
    PyObject_HEAD
    PyObject *contig;
    samfile_t *fd;
    bam_index_t *idx;
//    bam_plbuf_t *pileup;
//    PyObject *buffer;
    PyObject *callback;
} tactmod_BamObject;

// Bamfile object functions
static PyObject *Bam_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int Bam_init(tactmod_BamObject *self, PyObject *args, PyObject *kwds);

PyObject *Bam_slice(tactmod_BamObject *self, PyObject *args);
PyObject *Bam_counts(tactmod_BamObject *self, PyObject *args);
PyObject *Bam_pileup(tactmod_BamObject *self, PyObject *args);
PyObject *Bam_stats(tactmod_BamObject *self, PyObject *args);

void Bam_dealloc(tactmod_BamObject *self);

// Bamfile object implementation details
PyMethodDef Bam_methods[];
PyMemberDef Bam_members[];
PyTypeObject tactmod_BamType;

// The iterator
typedef struct {
    PyObject_HEAD
    uint32_t position;
    uint16_t offset;
    uint32_t start;
    uint32_t stop;
    tactmod_BamObject *bam;
    bam_plbuf_t *pileup;
    PyListObject *buffer; // the buffer linked list
} tactmod_BamIter;

typedef struct {
    void *next;
    PyObject *content;
} buffer_node;

typedef struct {
    uint32_t size;
    uint32_t end;
    buffer_node *head;
    buffer_node *tail;
} buffer_ll;


PyObject *tactmod_BamIter_iter(PyObject *self);
PyObject *tactmod_BamIter_next(PyObject *self);

PyTypeObject tactmod_BamIterType;
//int pileup_func(uint32_t tid, uint32_t pos, int n,
//                      const bam_pileup1_t *pl, void *data);

static int fetch_stats(const bam1_t *b, void *data);
static int fetch_column(const bam1_t *b, void *data);
static int fetch_pileup(const bam1_t *b, void *data);
static int fetch_counts(const bam1_t *b, void *data);
static int fetch_f(const bam1_t *b, void *data);

static int pileup_func(uint32_t tid, uint32_t pos, int n,
                       const bam_pileup1_t *pl, void *data);

int ll_push(buffer_ll *list, PyTupleObject *content);
PyTupleObject *ll_pop(buffer_ll *list);
int ll_init(buffer_ll *list);
int ll_destroy(buffer_ll *list);
#endif
