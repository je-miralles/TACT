#include <Python.h>
#include "structmember.h"

#include "multiseq.h"
#include "fasta.h"
#include "tact_bam.h"
#include "base.h"
#include "column.h"

/*
 *  TACT
 *  Python extension for genomic datatypes
 */

static PyMethodDef TgtmodMethods[] =
{
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC

inittactmod(void) {
    PyObject *m;
    
    tactmod_FastaType.tp_new = PyType_GenericNew;
    tactmod_MultiSeqIterType.tp_new = PyType_GenericNew;
    
    // Initialize TACT types
    if (PyType_Ready(&tactmod_FastaType) < 0) return;
    if (PyType_Ready(&tactmod_FastaIterType) < 0) return;
    if (PyType_Ready(&tactmod_BamType) < 0) return;
    if (PyType_Ready(&tactmod_BamIterType) < 0) return;
    if (PyType_Ready(&tactmod_BaseType) < 0) return;
    if (PyType_Ready(&tactmod_ColumnType) < 0) return;
    if (PyType_Ready(&tactmod_MultiSeqType) < 0) return;
    if (PyType_Ready(&tactmod_SequenceType) < 0) return;

    m = Py_InitModule("tactmod", TgtmodMethods);

    Py_INCREF(&tactmod_FastaType);
    Py_INCREF(&tactmod_FastaIterType);
    Py_INCREF(&tactmod_BamType);
    Py_INCREF(&tactmod_BamIterType);
    Py_INCREF(&tactmod_BaseType);
    Py_INCREF(&tactmod_ColumnType);
    Py_INCREF(&tactmod_MultiSeqType);
    Py_INCREF(&tactmod_MultiSeqIterType);

    PyModule_AddObject(m, "Fasta", (PyObject *)&tactmod_FastaType);
    PyModule_AddObject(m, "FastaIter", (PyObject *)&tactmod_FastaIterType);
    PyModule_AddObject(m, "Bam", (PyObject *)&tactmod_BamType);
    PyModule_AddObject(m, "BamIter", (PyObject *)&tactmod_BamIterType);
    PyModule_AddObject(m, "Base", (PyObject *)&tactmod_BaseType);
    PyModule_AddObject(m, "Column", (PyObject *)&tactmod_ColumnType);
    PyModule_AddObject(m, "MultiSeq", (PyObject *)&tactmod_MultiSeqType);
    PyModule_AddObject(m, "MultiSeqIter",
                       (PyObject *)&tactmod_MultiSeqIterType);
}

