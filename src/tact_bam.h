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

#define BUFFER_SIZE 1000000 // ~ 200MB
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
    PyTupleObject *targets;
    samfile_t *fd;
    bam_index_t *idx;
    bam_header_t *header;
//    bam_plbuf_t *pileup;
//    PyObject *buffer;
    PyDictObject *tids;
    PyObject *callback;
} tactmod_BamObject;

// Bamfile object functions
static tactmod_BamObject *Bam_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
int Bam_init(tactmod_BamObject *self, PyObject *args, PyObject *kwds);

PyObject *Bam_slice(tactmod_BamObject *self, PyObject *args);
PyObject *Bam_pileup(tactmod_BamObject *self, PyObject *args);
PyObject *Bam_stats(tactmod_BamObject *self, PyObject *args);

void Bam_dealloc(tactmod_BamObject *self);

// Bamfile object implementation details
PyMethodDef Bam_methods[];
PyMemberDef Bam_members[];
PyTypeObject tactmod_BamType;

typedef struct {
    uint32_t position;
    uint16_t depth;
    uint16_t features[5][6];
    uint16_t features_f[2];
    uint8_t major;
    uint8_t minor;
    int indels;
    double entropy;
    uint16_t ambiguous;
} column_t;

typedef struct {
    void *next;
    uint32_t position;
    column_t content;
} queue_node;

typedef struct {
    uint32_t size;
    uint32_t fetch_start;
    uint32_t fetch_stop;
    uint32_t end;
    uint32_t position;
    queue_node *head;
    queue_node *tail;
} queue;

// The iterator
typedef struct {
    PyObject_HEAD
    uint32_t position;
    uint16_t offset;
    uint32_t start;
    uint32_t stop;
    int tid;
    tactmod_BamObject *bam;
    PyObject *return_value;
    bam_plbuf_t *pileup;
    queue *buffer; // the buffer linked list
} tactmod_BamIter;


tactmod_BamIter *tactmod_BamIter_iter(tactmod_BamIter *self);
tactmod_BamIter *tactmod_BamIter_next(tactmod_BamIter *self);

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

int enqueue(queue *list, column_t content, uint32_t pos);
column_t dequeue(queue *list);
queue *queue_init(void);
int queue_destroy(queue *list);
tactmod_BamIter *Bam_counts(tactmod_BamObject *self, PyObject *args);

double binomial_ll(uint16_t k, uint16_t n, double mu);
double entropy(uint16_t bases[4], uint16_t depth);
#endif
