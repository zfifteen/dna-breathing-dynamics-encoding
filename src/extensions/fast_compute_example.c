/**
 * Example high-performance C extension for DNA breathing dynamics.
 *
 * This is a template for creating optimized C extensions with Apple Silicon
 * support. Uncomment and modify setup.py to build this extension.
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "common.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

/**
 * Example: Calculate GC content of DNA sequence (optimized).
 *
 * This demonstrates basic sequence analysis with SIMD-friendly operations.
 */
static PyObject* fast_gc_content(PyObject* self, PyObject* args) {
    const char* sequence;
    Py_ssize_t seq_len;

    // Parse Python string argument
    if (!PyArg_ParseTuple(args, "s#", &sequence, &seq_len)) {
        return NULL;
    }

    if (seq_len == 0) {
        PyErr_SetString(PyExc_ValueError, "Sequence cannot be empty");
        return NULL;
    }

    // Count GC nucleotides
    uint64_t gc_count = 0;
    for (Py_ssize_t i = 0; i < seq_len; i++) {
        nucleotide_t nuc = encode_nucleotide(sequence[i]);

        if (UNLIKELY(nuc == NUC_INVALID)) {
            PyErr_Format(PyExc_ValueError,
                "Invalid nucleotide '%c' at position %zd",
                sequence[i], i);
            return NULL;
        }

        if (is_gc(nuc)) {
            gc_count++;
        }
    }

    // Calculate GC content as fraction
    double gc_content = (double)gc_count / (double)seq_len;

    return PyFloat_FromDouble(gc_content);
}

/**
 * Example: Validate DNA sequence (fast path).
 */
static PyObject* fast_validate_sequence(PyObject* self, PyObject* args) {
    const char* sequence;
    Py_ssize_t seq_len;

    if (!PyArg_ParseTuple(args, "s#", &sequence, &seq_len)) {
        return NULL;
    }

    // Fast validation loop
    for (Py_ssize_t i = 0; i < seq_len; i++) {
        if (UNLIKELY(!is_valid_dna_nucleotide(sequence[i]))) {
            PyErr_Format(PyExc_ValueError,
                "Invalid nucleotide '%c' at position %zd",
                sequence[i], i);
            return NULL;
        }
    }

    Py_RETURN_TRUE;
}

// Method definitions
static PyMethodDef FastComputeMethods[] = {
    {
        "gc_content",
        fast_gc_content,
        METH_VARARGS,
        "Calculate GC content of DNA sequence (optimized).\n\n"
        "Args:\n"
        "    sequence: DNA sequence string\n\n"
        "Returns:\n"
        "    float: GC content as fraction [0, 1]"
    },
    {
        "validate_sequence",
        fast_validate_sequence,
        METH_VARARGS,
        "Validate DNA sequence (fast path).\n\n"
        "Args:\n"
        "    sequence: DNA sequence string\n\n"
        "Returns:\n"
        "    bool: True if valid\n\n"
        "Raises:\n"
        "    ValueError: If sequence contains invalid characters"
    },
    {NULL, NULL, 0, NULL}  // Sentinel
};

// Module definition
static struct PyModuleDef fastcomputemodule = {
    PyModuleDef_HEAD_INIT,
    "fast_compute",
    "High-performance C extensions for DNA analysis (Apple Silicon optimized)",
    -1,
    FastComputeMethods
};

// Module initialization
PyMODINIT_FUNC PyInit_fast_compute(void) {
    PyObject* module = PyModule_Create(&fastcomputemodule);
    if (module == NULL) {
        return NULL;
    }

    // Add platform information
    #ifdef APPLE_SILICON
        PyModule_AddIntConstant(module, "APPLE_SILICON", 1);
    #else
        PyModule_AddIntConstant(module, "APPLE_SILICON", 0);
    #endif

    return module;
}
