/*
This file contains the definitions of helper functions used in matrix multiplication.
These functions are called in the main multiplication functions.
*/

// Default C library
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// Intel Intrinsics
#include <immintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
// Threading
#include <pthread.h>

// Our headers
#include "matrixutils.h"
#include "utils.h"
#include "constants.h"


// MULTIPLICATION ALGORITHMS //
// ------------------------- //

void* multiply_main_implementation(void* void_arg) {
    // Multithreaded
    // Extract info from the given argument
    struct MultiplyArg* arg = (struct MultiplyArg*) void_arg;  // to fit thread creation signature
    Matrix* matrix_a = arg->matrix_a;
    Matrix* matrix_b = arg->matrix_b;
    Matrix* matrix_result = arg->matrix_result;
    uint64_t start_row = arg->start_row;
    uint64_t end_row = arg->end_row;

    // Iterating the rows of A / rowA := row index of A
    matrix_result->rowPointers[0] = 0;
    for (uint64_t rowA = start_row; rowA < end_row; rowA++) {
        // Getting the row start and end
        uint64_t rowABeg = matrix_a->rowPointers[rowA];
        uint64_t rowAEnd = matrix_a->rowPointers[rowA + 1];

        // Iterating over the columns of A in line rowA
        for (uint64_t indexA = rowABeg; indexA < rowAEnd; indexA++) {
            float valueA = matrix_a->values[indexA];

            // rowB = columnA th row of B
            uint64_t rowBBeg = matrix_b->rowPointers[matrix_a->colIndices[indexA]];
            uint64_t rowBEnd = matrix_b->rowPointers[matrix_a->colIndices[indexA] + 1];

            float valueB = 0;
            float valueC = 0;

            uint64_t columnB = 0;
            uint64_t indexB = rowBBeg;
            // Iterating over the columns of B                 
            for (indexB = rowBBeg; indexB < rowBEnd; indexB++) {
                valueB = matrix_b->values[indexB];
                // Actual multiplication of the values
                valueC = valueA * valueB;
                columnB = matrix_b->colIndices[indexB];
                // Adding the valueC to the exsting value in Matrix C
                matrix_result->values[rowA * matrix_result->noCols + columnB] += valueC;
                // Updating colIndices
                matrix_result->colIndices[rowA * matrix_b->noCols + columnB] = columnB;
            }
        }
        // Updating row pointer -> end of this row => rowA +1 Eintrag (0. -> 1. , 1. -> 2.)
        matrix_result->rowPointers[rowA + 1] = (rowA + 1) * matrix_result->noCols;
    }

    return NULL;  // this is required for pthread_create()
}

float** csr_to_ordinary(const Matrix* const csr) {
    // Initialize memory for the rows of the 2D array
    float** ordinary = malloc_safe(csr->noRows, sizeof(float*));
    if (ordinary == NULL) {
        return NULL;
    }

    for (uint64_t i = 0; i < csr->noRows; i++) {
        // Initialize a row's columns with zero values
        ordinary[i] = calloc(csr->noCols, sizeof(float));
        if (ordinary[i] == NULL) {
            for (uint64_t j = 0; j < i; j++) {
                free(ordinary[j]);
            }
            free(ordinary);
            return NULL;
        }
    }
    // writing CSR data to the 2D-array
    for (uint64_t i = 0; i < csr->noRows; i++) {
        for (uint64_t j = csr->rowPointers[i]; j < csr->rowPointers[i+1]; j++) {
            uint64_t col = csr->colIndices[j];
            ordinary[i][col] = csr->values[j];
        }
    }
    return ordinary;
}

Matrix* ordinary_to_csr(const uint64_t noRows, const uint64_t noCols, float** ordinary) {
    Matrix* csr = malloc(sizeof(Matrix));
    if (csr == NULL) {
        return NULL;
    }

    uint64_t no_non_zero = 0;
    for (uint64_t i = 0; i < noRows; i++) {
        for (uint64_t j = 0; j < noCols; j++) {
            if (ordinary[i][j] != 0) {
                no_non_zero++;
            }
        }

    }

    csr->noRows = noRows;
    csr->noCols = noCols;

    float* values = malloc_safe(sizeof(float), no_non_zero);
    if (values == NULL) {
        return NULL;
    }
    csr->valuesSize = no_non_zero;

    uint64_t* colIndices = malloc_safe(sizeof(uint64_t), no_non_zero);
    if (colIndices == NULL) {
        free(values);
        return NULL;
    }

    uint64_t* rowPointers = malloc_safe(sizeof(uint64_t), noRows + 1);
    if (rowPointers == NULL) {
        free(values);
        free(colIndices);
        return NULL;
    }

    uint64_t c = 0; // insertion counter
    rowPointers[0] = 0;  // the first row pointer is always 0
    for (uint64_t i = 0; i < noRows; i++) {
        for (uint64_t j = 0; j < noCols; j++) {
            if (ordinary[i][j] != 0) {
                values[c] = ordinary[i][j];
                colIndices[c] = j;
                c++;
            }
        }
        rowPointers[i + 1] = c;
    }

    csr->values = values;
    csr->colIndices = colIndices;
    csr->rowPointers = rowPointers;
    csr->rowPointersSize = noRows + 1;

    return csr;
}

void multiply_V2(
    const Matrix* const restrict matrix_a, 
    const Matrix* const restrict matrix_b, 
    Matrix* const restrict matrix_result
    ) {
    // gustafson ohne prediction
    // Iterating the rows of A / rowA := row index of A
    matrix_result->rowPointers[0] = 0;
    for (uint64_t rowA = 0; rowA < matrix_a->noRows; rowA++) {
        // Getting the row start and end
        uint64_t rowABeg = matrix_a->rowPointers[rowA];
        uint64_t rowAEnd = matrix_a->rowPointers[rowA + 1];

        // Iterating over the columns of A in line rowA
        for (uint64_t indexA = rowABeg; indexA < rowAEnd; indexA++) {
            float valueA = matrix_a->values[indexA];

            // rowB = columnA th row of B
            uint64_t rowBBeg = matrix_b->rowPointers[matrix_a->colIndices[indexA]];
            uint64_t rowBEnd = matrix_b->rowPointers[matrix_a->colIndices[indexA] + 1];

            float valueB = 0;
            float valueC = 0;

            uint64_t columnB = 0;
            uint64_t indexB = rowBBeg;
            // Iterating over the columns of B                 
            for (indexB = rowBBeg; indexB < rowBEnd; indexB++) {
                valueB = matrix_b->values[indexB];
                // Actual multiplication of the values
                valueC = valueA * valueB;
                columnB = matrix_b->colIndices[indexB];
                // Adding the valueC to the exsting value in Matrix C
                matrix_result->values[rowA * matrix_result->noCols + columnB] += valueC;
                // Updating colIndices
                matrix_result->colIndices[rowA * matrix_b->noCols + columnB] = columnB;
            }
        }
        // Updating row pointer -> end of this row => rowA +1 Eintrag (0. -> 1. , 1. -> 2.)
        matrix_result->rowPointers[rowA + 1] = (rowA + 1) * matrix_result->noCols;
    }
}

int multiply_V3(
    const Matrix* const restrict matrix_a, 
    const Matrix* const restrict matrix_b, 
    Matrix* const restrict matrix_result
    ) {
    // Iterating the rows of A / rowA := row index of A
    matrix_result->rowPointers[0] = 0;
    float* values_to_add = malloc(sizeof(float) * 4);  // holds 4 float values for SIMD
    if (values_to_add == NULL) {
        return HEAP_MEMORY_ERROR;
    }

    for (uint64_t rowA = 0; rowA < matrix_a->noRows; rowA++) {
        // Getting the row start and end
        uint64_t rowABeg = matrix_a->rowPointers[rowA];
        uint64_t rowAEnd = matrix_a->rowPointers[rowA + 1];

        uint64_t resultValueIndex;
        uint64_t resultColIndex;
        __m128 valueA;
        uint64_t rowBBeg;
        uint64_t rowBEnd;
        uint64_t columnB;
        uint64_t indexB;
        uint64_t valuesToMultiply;
        __m128 valuesB;

        // Iterating over the columns of A in line rowA
        for (uint64_t indexA = rowABeg; indexA < rowAEnd; indexA++) {
            resultValueIndex = rowA * matrix_result->noCols;
            resultColIndex = rowA * matrix_b->noCols;
            // rowB = columnA th row of B
            rowBBeg = matrix_b->rowPointers[matrix_a->colIndices[indexA]];
            rowBEnd = matrix_b->rowPointers[matrix_a->colIndices[indexA] + 1];

            indexB = rowBBeg;

            valuesToMultiply = rowBEnd - rowBBeg;
            if (valuesToMultiply >= 4) {
                // Load in the last lane
                valueA = _mm_load_ss(matrix_a->values + indexA);
                // Copy it in all 4 lanes
                valueA = _mm_shuffle_ps(valueA, valueA, 0x00);

                do {
                    // Multiply b values with a
                    valuesB = _mm_loadu_ps(matrix_b->values + indexB);
                    valuesB = _mm_mul_ps(valuesB, valueA);
                    _mm_storeu_ps(values_to_add, valuesB);
                    for (int i = 0; i < 4; i++) {
                        columnB = matrix_b->colIndices[indexB++];
                        matrix_result->values[resultValueIndex + columnB] += values_to_add[i];
                        matrix_result->colIndices[resultColIndex + columnB] = columnB;
                    }
                    valuesToMultiply -= 4;
                } while (valuesToMultiply >= 4);
            }

            float valueA = matrix_a->values[indexA];
            float valueB;
            float valueC;

            // Iterating over the columns of B                 
            for (; indexB < rowBEnd; indexB++) {
                valueB = matrix_b->values[indexB];
                // Actual multiplication of the values
                valueC = valueA * valueB;
                columnB = matrix_b->colIndices[indexB];
                // Adding the valueC to the exsting value in Matrix C
                matrix_result->values[resultValueIndex + columnB] += valueC;

                // Updating colIndices
                matrix_result->colIndices[resultColIndex + columnB] = columnB;
            }
        }

        // Updating row pointer -> end of this row => rowA +1 Eintrag (0. -> 1. , 1. -> 2.)
        matrix_result->rowPointers[rowA + 1] = (rowA + 1) * matrix_result->noCols;
    }

    free(values_to_add);  // free SIMD array

    return 0;
}

int multiply_V4(
    const Matrix* const restrict matrix_a, 
    const Matrix* const restrict matrix_b, 
    Matrix* const restrict matrix_result
    ) {
    // Iterating the rows of A / rowA := row index of A
    matrix_result->rowPointers[0] = 0;
    float* values_to_add = malloc(sizeof(float) * 8);  // holds 8 float values for SIMD
    if (values_to_add == NULL) {
        return HEAP_MEMORY_ERROR;
    }

    for (uint64_t rowA = 0; rowA < matrix_a->noRows; rowA++) {
        // Getting the row start and end
        uint64_t rowABeg = matrix_a->rowPointers[rowA];
        uint64_t rowAEnd = matrix_a->rowPointers[rowA + 1];

        uint64_t resultValueIndex;
        uint64_t resultColIndex;
        __m256 valueA;
        uint64_t rowBBeg;
        uint64_t rowBEnd;
        uint64_t columnB;
        uint64_t indexB;
        uint64_t valuesToMultiply;
        __m256 valuesB;


        // Iterating over the columns of A in line rowA
        for (uint64_t indexA = rowABeg; indexA < rowAEnd; indexA++) {
            resultValueIndex = rowA * matrix_result->noCols;
            resultColIndex = rowA * matrix_b->noCols;
            // rowB = columnA th row of B
            rowBBeg = matrix_b->rowPointers[matrix_a->colIndices[indexA]];
            rowBEnd = matrix_b->rowPointers[matrix_a->colIndices[indexA] + 1];

            indexB = rowBBeg;

            valuesToMultiply = rowBEnd - rowBBeg;
            if (valuesToMultiply >= 8) {
                // Load value a in all 8 lanes
                valueA = _mm256_broadcast_ss(matrix_a->values + indexA);

                do {
                    // Multiply b values with a
                    valuesB = _mm256_loadu_ps(matrix_b->values + indexB);
                    valuesB = _mm256_mul_ps(valuesB, valueA);
                    _mm256_storeu_ps(values_to_add, valuesB);
                    for (int i = 0; i < 8; i++) {
                        columnB = matrix_b->colIndices[indexB++];
                        matrix_result->values[resultValueIndex + columnB] += values_to_add[i];
                        matrix_result->colIndices[resultColIndex + columnB] = columnB;
                    }
                    valuesToMultiply -= 8;
                } while (valuesToMultiply >= 8);
            }

            float valueA = matrix_a->values[indexA];
            float valueB;
            float valueC;

            // Iterating over the columns of B                 
            for (; indexB < rowBEnd; indexB++) {
                valueB = matrix_b->values[indexB];
                // Actual multiplication of the values
                valueC = valueA * valueB;
                columnB = matrix_b->colIndices[indexB];
                // Adding the valueC to the exsting value in Matrix C
                matrix_result->values[resultValueIndex + columnB] += valueC;

                // Updating colIndices
                matrix_result->colIndices[resultColIndex + columnB] = columnB;
            }
        }
        // Updating row pointer -> end of this row => rowA +1 Eintrag (0. -> 1. , 1. -> 2.)
        matrix_result->rowPointers[rowA + 1] = (rowA + 1) * matrix_result->noCols;
    }

    free(values_to_add);  // free SIMD array

    return 0;
}

void multiply_V5(
    const Matrix* const restrict matrix_a, 
    const Matrix* const restrict matrix_b, 
    Matrix* const restrict matrix_result
    ) {
    // Gustafson with size prediction
    uint64_t valuesBegPtr = 0;
    uint64_t valuesEndPtr = 0; 
    // Iterating the rows of A / rowA := row index of A
    matrix_result->rowPointers[0] = 0;
    for (uint64_t rowA = 0; rowA < matrix_a->noRows; rowA++) {
        // Getting the row start and end
        uint64_t rowABeg = matrix_a->rowPointers[rowA];
        uint64_t rowAEnd = matrix_a->rowPointers[rowA+1]; 

        // Iterating over the columns of A in line rowA
        for (uint64_t indexA = rowABeg; indexA < rowAEnd; indexA++) {
            float valueA = matrix_a->values[indexA];

            // rowB = columnA th row of B
            uint64_t rowBBeg = matrix_b->rowPointers[matrix_a->colIndices[indexA]];
            uint64_t rowBEnd = matrix_b->rowPointers[matrix_a->colIndices[indexA]+1];

            float valueB = 0;  
            float valueC = 0; 

            uint64_t columnB = 0; 
            uint64_t indexB = rowBBeg;

            // Iterating over the columns of B                 
            for (indexB = rowBBeg; indexB < rowBEnd; indexB++) {
                valueB = matrix_b->values[indexB];
                // Actual multiplication of the values
                valueC = valueA * valueB;
                columnB = matrix_b->colIndices[indexB];

                for (uint64_t i = valuesBegPtr; i < matrix_result->valuesSize; i++) {
                    if (matrix_result->values[i] == 0){
                        // Füge Element ein 
                        matrix_result->values[i] = valueC; 
                        matrix_result->colIndices[i] = columnB; 
                        valuesEndPtr = i+1;  
                        break; 
                    } else if(matrix_result->colIndices[i] == columnB){
                        matrix_result->values[i] += valueC; 
                        break;
                    }       
                }

            }             
        }
        // Updating row pointer -> end of this row => rowA +1 Eintrag (0. -> 1. , 1. -> 2.)
        matrix_result->rowPointers[rowA+1] = valuesEndPtr; 
        // Updating the beg Ptr
        valuesBegPtr = valuesEndPtr; 
    }
}

// OTHER HELPER FUNCTIONS BELOW //
// ---------------------------- //

int start_threads(
    const unsigned int thread_count, pthread_t** threads, struct MultiplyArg*** arguments,
    Matrix* matrix_a, Matrix* matrix_b, Matrix* matrix_result
    ) {
    // Create the multiplication arguments
    if (!thread_count) {  // this is required to avoid a maybe uninitialized warning from gcc
        return THREAD_START_ERROR;
    }

    *arguments = malloc(sizeof(struct MultiplyArg*) * thread_count);  // no need to check for overflow, thread_count is a small number
    if (*arguments == NULL) {
        return HEAP_MEMORY_ERROR;
    }

    uint64_t prev = 0;
    uint64_t step = matrix_a->noRows / thread_count;
    unsigned int rest = matrix_a->noRows % thread_count;  // rest is also in this range so unsigned int

    for (unsigned int i = 0; i < thread_count; i++) {
        struct MultiplyArg* arg = malloc(sizeof(struct MultiplyArg));
        if (arg == NULL) {
            for (unsigned int j = 0; j < i; j++) {
                free((*arguments)[j]);
            }
            free(*arguments);
            return HEAP_MEMORY_ERROR;
        }

        arg->matrix_a = matrix_a;
        arg->matrix_b = matrix_b;
        arg->matrix_result = matrix_result;

        arg->start_row = prev;
        if (i < rest) {
            prev += step + 1;
            arg->end_row = prev;
        } else {
            prev += step;
            arg->end_row = prev;
        }
        

        (*arguments)[i] = arg;
    }
    (*arguments)[thread_count-1]->end_row = matrix_a->noRows;

    // Create the threads
    *threads = malloc(sizeof(pthread_t) * thread_count);
    if (*threads == NULL) {
        return HEAP_MEMORY_ERROR;
    }

    pthread_t thread;
    for (unsigned int i = 0; i < thread_count; i++) {
        if (pthread_create(&thread, NULL, &multiply_main_implementation, (*arguments)[i])) {
            // Error creating thread, join previous threads and clean up memory
            for (unsigned int j = 0; j < i; j++) {
                pthread_join((*threads)[j], NULL);
            }

            for (unsigned int j = 0; j < thread_count; j++) {
                free((*arguments)[j]);
            }
            free(*arguments);

            return THREAD_START_ERROR;
        } else {
            (*threads)[i] = thread;
        }
    }

    return 0;
}

int can_multiply(const Matrix* const a, const Matrix* const b) {
    return (a->noCols == b->noRows) && (a->noCols >= 1 && a->noRows >= 1) && (b->noCols >= 1 && b->noRows >= 1);
}

int clean_up_csr(Matrix* const matrix) {
    uint64_t non_zero_values;  // count the non zero values
    if (_clean_up_matrix_arrays(matrix, &non_zero_values) == HEAP_MEMORY_ERROR) {
        return HEAP_MEMORY_ERROR;
    }

    if (non_zero_values > 0) {
        // Don't free arrays on error, free_csr_matrix() does that in main.c, no need for realloc_safe()
        uint64_t alloc_size;
        if (__builtin_umull_overflow(non_zero_values, sizeof(float), &alloc_size)) {
            return HEAP_MEMORY_ERROR;
        }

        float* tmp_values = matrix->values;
        matrix->values = realloc(matrix->values,  alloc_size);
        if (matrix->values == NULL) {
            matrix->values = tmp_values;
            return HEAP_MEMORY_ERROR;
        }

        
        if (__builtin_umull_overflow(non_zero_values, sizeof(uint64_t), &alloc_size)) {
            return HEAP_MEMORY_ERROR;
        }

        uint64_t* tmp_col_indices = matrix->colIndices;
        matrix->colIndices = realloc(matrix->colIndices, alloc_size);
        if (matrix->colIndices == NULL) {
            // don't free arrays on error, free_csr_matrix() does that in main.c
            matrix->colIndices = tmp_col_indices;
            return HEAP_MEMORY_ERROR;
        }
    }

    matrix->valuesSize = non_zero_values;  // this ensures an empty line is printed when nnz = 0

    return 0;
}

int _clean_up_matrix_arrays(Matrix* matrix, uint64_t* non_zero_values) {
    *non_zero_values = 0;
    uint64_t removed = 0;  // how many zero values were removed?
    uint64_t row_beg = matrix->rowPointers[0];
    uint64_t row_end;
    for (uint64_t row = 1; row < matrix->rowPointersSize; row++) {
        // In the start of the iteration, row_beg is the old row_end
        row_end = matrix->rowPointers[row];
        while (row_beg < row_end) {
            // printf("row = %lu\n", row_beg);
            if (matrix->values[row_beg]) {
                // Move non zero value to next free spot
                matrix->values[*non_zero_values] = matrix->values[row_beg];
                matrix->colIndices[(*non_zero_values)++] = matrix->colIndices[row_beg];
            } else {
                // Zero value, this will be removed/overwritten
                removed++;
            }

            row_beg++;
        }
        // Calculate what the new row pointer value should be
        matrix->rowPointers[row] = row_end - removed;
    }

    return 0;
}

int predict_values_dimension(
    const Matrix* const restrict matrix_a, 
    const Matrix* const restrict matrix_b, 
    uint64_t* const values_size
    ) {
    *values_size = 0;
    uint64_t rowBBeg;
    uint64_t rowBEnd;
    uint64_t maxValueSize = matrix_a->noRows * matrix_b->noCols;

    // noCols of a = noRows of b
    uint64_t* numValuesInColsA = calloc(matrix_a->noCols, sizeof(uint64_t));
    if (numValuesInColsA == NULL) {
        return HEAP_MEMORY_ERROR;
    }
    uint64_t* numValuesInRowB = calloc(matrix_b->noRows, sizeof(uint64_t));
    if (numValuesInRowB == NULL) {
        free(numValuesInColsA);
        return HEAP_MEMORY_ERROR;
    }
    uint64_t numsSize = matrix_a->noCols;

    for (uint64_t rowB = 0; rowB < matrix_b->noRows; rowB++) {
        rowBBeg = matrix_b->rowPointers[rowB];
        rowBEnd = matrix_b->rowPointers[rowB + 1];
        // Calculating the nnz in rowB -> = rowEnd - rowBeg because CSR
        numValuesInRowB[rowB] = (rowBEnd - rowBBeg);
    }

    // numValuesInColsA: Iteration of the nnz values in a  
    for (uint64_t valuesIndex = 0; valuesIndex < matrix_a->valuesSize; valuesIndex++) {
        uint64_t column = matrix_a->colIndices[valuesIndex];
        // Updating the column count
        numValuesInColsA[column]++;;
    }

    for (uint64_t index = 0; index < numsSize; index++) {
        *values_size += numValuesInRowB[index] * numValuesInColsA[index];
    }

    free(numValuesInColsA);
    free(numValuesInRowB);

    // Check if the prediction is greater than the maximum
    *values_size = *values_size > maxValueSize ? maxValueSize : *values_size;
    return 0;
}

int init_empty_csr_matrix(
    const Matrix* const restrict matrix_a, const Matrix* const restrict matrix_b, 
    Matrix* const restrict result, const int predict_flag
    ) {
    uint64_t valuesSize;
    if (predict_flag) {
        // Predict size mathematically (upper limit)
        if (predict_values_dimension(matrix_a, matrix_b, &valuesSize) == HEAP_MEMORY_ERROR) {
            return HEAP_MEMORY_ERROR;
        }
    } else {
        // Set values array size to maximum number of elements in the matrix
        valuesSize = matrix_a->noRows * matrix_b->noCols;
    }

    float* values = calloc(valuesSize, sizeof(float));
    if (values == NULL) {
        init_matrix_error:  // set result's subarrays to NULL for safety
        result->values = NULL;
        result->colIndices = NULL;
        result->rowPointers = NULL;
        return HEAP_MEMORY_ERROR;
    }

    uint64_t* colIndices = calloc(valuesSize, sizeof(uint64_t));
    if (colIndices == NULL) {
        free(values);
        goto init_matrix_error;
    }

    uint64_t* rowPointers = malloc_safe(sizeof(uint64_t), matrix_a->noRows + 1);
    if (rowPointers == NULL) {
        free(colIndices);
        free(values);
        goto init_matrix_error;
    }

    result->noRows = matrix_a->noRows;
    result->noCols = matrix_b->noCols;
    result->values = values;
    result->valuesSize = valuesSize;
    result->colIndices = colIndices;
    result->rowPointers = rowPointers;
    result->rowPointersSize = result->noRows + 1;

    return 0;
}

// HELPFUL TESTING FUNCTIONS FOR TUTORS BELOW //
// ------------------------------------------ //

int equals(Matrix* const restrict matrix_a, Matrix* const restrict matrix_b) {
    if (matrix_a == NULL || matrix_b == NULL) {
        return 0;
    }

    if (matrix_a->noCols != matrix_b->noCols ||
        matrix_a->noRows != matrix_b->noRows ||
        matrix_a->valuesSize != matrix_b->valuesSize ||
        matrix_a->rowPointersSize != matrix_b->rowPointersSize
        ) {
        return 0;
    }

    // Ascendingly sort the values in the matrices using the column indices
    for (uint64_t i = 0; i < matrix_a->rowPointersSize - 1; i++) {
        uint64_t rowBegA = matrix_a->rowPointers[i];
        uint64_t rowEndA = matrix_a->rowPointers[i + 1] - 1;
        uint64_t rowBegB = matrix_b->rowPointers[i];
        uint64_t rowEndB = matrix_b->rowPointers[i + 1] - 1;
        sort_matrix(matrix_a, rowBegA, rowEndA);
        sort_matrix(matrix_b, rowBegB, rowEndB);
    }

    // Check if values and column indices are equal
    for (uint64_t i = 0; i < matrix_a->valuesSize; i++) {
        if (
            matrix_a->colIndices[i] != matrix_b->colIndices[i] ||
            matrix_a->values[i] != matrix_b->values[i]
            ) {
            return 0;
        }
    }

    // Check if row pointers are equal
    for (uint64_t i = 0; i < matrix_a->rowPointersSize; i++) {
        if (matrix_a->rowPointers[i] != matrix_b->rowPointers[i]) {
            return 0;
        }
    }

    return 1;
}

void sort_matrix(Matrix* const matrix, const uint64_t beg, const uint64_t end) {
    for (uint64_t i = beg; i <= end; i++) {
        uint64_t col = matrix->colIndices[i];
        uint64_t currentPositionOfValue = i;
        for (uint64_t j = i - 1; j >= beg; j--) {
            if (col < matrix->colIndices[j]) {
                _swap_cols(matrix->colIndices, matrix->valuesSize, currentPositionOfValue, j);
                _swap_values(matrix->values, matrix->valuesSize, currentPositionOfValue, j);
                currentPositionOfValue = j;
            }
            else {
                break;
            }
        }
    }
}

void print_ordinary_matrix(float** matrix, const uint64_t noRows, const uint64_t noCols) {
    for (uint64_t i = 0; i < noRows; i++) {
        for (uint64_t j = 0; j < noCols; j++) {
            printf("%.2f\n", matrix[i][j]);
        }
    }
}

void print_csr_matrix(const Matrix* const csr) {
    float** ordinaryMatrix = csr_to_ordinary(csr);
    if (ordinaryMatrix == NULL)
    {
        fprintf(stderr, "Error converting matrix\n");
    }

    print_ordinary_matrix(ordinaryMatrix, csr->noRows, csr->noCols);

    for (uint64_t i = 0; i < csr->noRows; i++)
    {
        free(ordinaryMatrix[i]);
    }
    free(ordinaryMatrix);
}

void _swap_cols(uint64_t* const cols, const uint64_t size, const uint64_t i, const uint64_t j) {
    if (i < size && j < size) {
        uint64_t tmp = cols[j];
        cols[j] = cols[i];
        cols[i] = tmp;
    }
}

void _swap_values(float* const values, const uint64_t size, const uint64_t i, const uint64_t j) {
    if (i < size && j < size) {
        float tmp = values[j];
        values[j] = values[i];
        values[i] = tmp;
    }
}
