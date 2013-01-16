#include <Python.h>
#include "structmember.h"

#include "multiseq.h"
#include "fasta.h"
#include "tact_bam.h"
#include "base.h"
#include "column.h"
#include "tact_vcf.h"
#include "tact_fastq.h"

/* Python extension for genomic datatypes */

static PyMethodDef TgtmodMethods[] =
{
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC

inittactmod(void) {
    PyObject *m;
    
    tactmod_FastaType.tp_new = PyType_GenericNew;
//    tactmod_MultiSeqIterType.tp_new = PyType_GenericNew;
   
    if (PyType_Ready(&tactmod_FastaType) < 0) return;
    if (PyType_Ready(&tactmod_FastaIterType) < 0) return;
    if (PyType_Ready(&tactmod_BamType) < 0) return;
    if (PyType_Ready(&tactmod_BamIterType) < 0) return;
    if (PyType_Ready(&tactmod_BaseType) < 0) return;
    if (PyType_Ready(&tactmod_ColumnType) < 0) return;
    if (PyType_Ready(&tactmod_MultiSeqType) < 0) return;
    if (PyType_Ready(&tactmod_MultiSeqIterType) < 0) return;
    if (PyType_Ready(&tactmod_FastqType) < 0 ) return;
    if (PyType_Ready(&tactmod_FastqIterType) < 0) return;
    if (PyType_Ready(&tactmod_VcfType) < 0) return;
    if (PyType_Ready(&tactmod_VcfIterType) < 0) return;

    m = Py_InitModule("tactmod", TgtmodMethods);

    Py_INCREF(&tactmod_FastaType);
    Py_INCREF(&tactmod_FastaIterType);
    Py_INCREF(&tactmod_BamType);
    Py_INCREF(&tactmod_BamIterType);
    Py_INCREF(&tactmod_BaseType);
    Py_INCREF(&tactmod_ColumnType);
    Py_INCREF(&tactmod_MultiSeqType);
    Py_INCREF(&tactmod_MultiSeqIterType);
    Py_INCREF(&tactmod_FastqType);
    Py_INCREF(&tactmod_FastqIterType);

    PyModule_AddObject(m, "Fasta", (PyObject *)&tactmod_FastaType);
    PyModule_AddObject(m, "FastaIter", (PyObject *)&tactmod_FastaIterType);
    PyModule_AddObject(m, "Bam", (PyObject *)&tactmod_BamType);
    PyModule_AddObject(m, "BamIter", (PyObject *)&tactmod_BamIterType);
    PyModule_AddObject(m, "Base", (PyObject *)&tactmod_BaseType);
    PyModule_AddObject(m, "Column", (PyObject *)&tactmod_ColumnType);
    PyModule_AddObject(m, "MultiSeq", (PyObject *)&tactmod_MultiSeqType);
    PyModule_AddObject(m, "MultiSeqIter", 
                       (PyObject *)&tactmod_MultiSeqIterType);
    PyModule_AddObject(m, "Fastq", (PyObject *)&tactmod_FastqType);
    PyModule_AddObject(m, "FastqIter", (PyObject *)&tactmod_FastqType);
    PyModule_AddObject(m, "Vcf", (PyObject *)&tactmod_VcfType);
    PyModule_AddObject(m, "VcfIter", (PyObject *)&tactmod_VcfIterType);

    Py_INCREF(&tact_A);
    tact_A = PyObject_CallObject((PyObject*)&tactmod_BaseType, NULL);
    tact_A->nucleotides = 0x01;
    Py_INCREF(&tact_C);
    tact_C = PyObject_CallObject((PyObject*)&tactmod_BaseType, NULL);
    tact_C->nucleotides = 0x02;
    Py_INCREF(&tact_G);
    tact_G = PyObject_CallObject((PyObject*)&tactmod_BaseType, NULL);
    tact_G->nucleotides = 0x04;
    Py_INCREF(&tact_T);
    tact_T = PyObject_CallObject((PyObject*)&tactmod_BaseType, NULL);
    tact_T->nucleotides = 0x08;
}

