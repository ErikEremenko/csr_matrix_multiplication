#ifndef CSRMATRIX_H
#define CSRMATRIX_H

#include <stdint.h>
#include <stddef.h>

// This file contains the definition for the CSR matrix struct.

/*
The struct Matrix represents a CSR Matrix.

noCols is the number of columns in the matrix.
noRows is the number of rows in the matrix.
values is the array of non-zero values in the matrix.
valuesSize is the size of the 'values' array.
colIndices is an array that contains the index of every value in the matrix.
rowPointers is an array that indicates where a row starts and where it ends.
rowPointersSize is the size of the 'rowPointers' array.

values, colIndices and rowPointers are commonly referred to as a matrix'
subarrays in the documentation.
*/
typedef struct Matrix {
    uint64_t noRows;
    uint64_t noCols;

    float* values;
    uint64_t valuesSize;

    uint64_t* colIndices;
    // no need for colIndicesSize because it is the same as valuesSize

    uint64_t* rowPointers;
    uint64_t rowPointersSize;
} Matrix;

#endif