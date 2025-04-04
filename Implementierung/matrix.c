// Default C library headers
#include <stdlib.h>
#include <string.h>
#include <errno.h>

// Threading
#include <pthread.h>
// System info for threading
#include <sys/sysinfo.h>
// Intrinsics
#include <immintrin.h>

// Our headers
#include "utils.h"
#include "matrixutils.h"
#include "csrmatrix.h"
#include "constants.h"
#include "matrix.h"


void matr_mult_csr(const void* a, const void* b, void* result) {
    // Multithreading implementation (Hauptimplementierung)
    Matrix* matrix_a = (Matrix*) a;
    if (matrix_a->valuesSize < 10000 || matrix_a->noRows < MIN_THREADS) {
        // Switch to Gustafson with size prediction if threading doesn't pay off
        matr_mult_csr_V5(a, b, result);
        return;
    }

    Matrix* matrix_b = (Matrix*) b;
    Matrix* matrix_result = (Matrix*) result;

    // Check if matrices are compatible
    if (!can_multiply(matrix_a, matrix_b)) {
        errno = MATRIX_DIMENSION_ERROR;
        return;
    }

    // Initialize matrix subarrays with no prediction
    if (init_empty_csr_matrix(matrix_a, matrix_b, matrix_result, NO_PREDICTION) == HEAP_MEMORY_ERROR) {
        errno = HEAP_MEMORY_ERROR;
        return;
    }

    // Deciding on thread count based on number of system threads as well as the number values and rows in matrix
    uint64_t available_threads = get_nprocs();
    uint64_t row_weight = matrix_a->noRows / THREAD_COND_MIN_VALUES;
    uint64_t elem_weight = matrix_a->noRows / THREAD_COND_MIN_ROWS;
    // Taking the minimum of both conditional parameters as thread count
    uint64_t tmp_min = row_weight < elem_weight ? row_weight : elem_weight;
    uint64_t thread_count = tmp_min < available_threads ? tmp_min : available_threads;

    thread_count = thread_count < MIN_THREADS ? MIN_THREADS : thread_count;

    // Start threads
    pthread_t* threads;
    struct MultiplyArg** arguments;
    int start_result = start_threads(
        thread_count, &threads, &arguments,
        matrix_a, matrix_b, matrix_result
    );
    if (start_result == THREAD_START_ERROR) {
        errno = THREAD_START_ERROR;
        return;
    }

    // Join threads and clean up memory
    for (unsigned int i = 0; i < thread_count; i++) {
        pthread_join(threads[i], NULL);
        free(arguments[i]);
    }
    free(arguments);
    free(threads);

    // Clean up nnz in matrix
    if (clean_up_csr(matrix_result) == HEAP_MEMORY_ERROR) {
        errno = HEAP_MEMORY_ERROR;
    }
}

void matr_mult_csr_V1(const void* a, const void* b, void* result) {
    // CSR to 2D array implementation
    Matrix* matrix_a = (Matrix*) a;
    Matrix* matrix_b = (Matrix*) b;
    Matrix* matrix_result = (Matrix*) result;

    // Check if multiplication is mathematically defined
    if (!can_multiply(matrix_a, matrix_b)) {
        errno = MATRIX_DIMENSION_ERROR;
        return;
    }

    // Convert the matrices to 2D array form
    float** ord_a = csr_to_ordinary(matrix_a);
    if (ord_a == NULL) {
        errno = MATRIX_CONVERSION_ERROR;
        return;
    }

    float** ord_b = csr_to_ordinary(matrix_b);
    if (ord_b == NULL) {
        // error, free converted a
        for (uint64_t i = 0; i < matrix_a->noRows; i++) {
            free(ord_a[i]);
        }
        free(ord_a);

        errno = MATRIX_CONVERSION_ERROR;
        return;
    }

    // Allocate memory to the resulting 2D array
    float** ord_result = malloc_safe(matrix_a->noRows, sizeof(float*));
    if (ord_result == NULL) {
        // Error, free converted a and b
        for (uint64_t i = 0; i < matrix_a->noRows; i++) {
            free(ord_a[i]);
        }
        free(ord_a);

        for (uint64_t i = 0; i < matrix_b->noRows; i++) {
            free(ord_b[i]);
        }
        free(ord_b);

        errno = HEAP_MEMORY_ERROR;
        return;
    }

    // Allocate memory for the columns in each row in the result matrix
    for (uint64_t i = 0; i < matrix_a->noRows; i++) {
        ord_result[i] = calloc(matrix_b->noCols, sizeof(float));
        if (ord_result[i] == NULL) {
            for (uint64_t j = 0; j < i; j++) {
                free(ord_result[j]);
            }
            free(ord_result);
            for (uint64_t i = 0; i < matrix_a->noRows; i++) {
                free(ord_a[i]);
            }
            free(ord_a);

            for (uint64_t i = 0; i < matrix_b->noRows; i++) {
                free(ord_b[i]);
            }
            free(ord_b);

            errno = HEAP_MEMORY_ERROR;
            return;
        }
    }

    // Sequential multiplication
    for (uint64_t i = 0; i < matrix_a->noRows; i++) {
        for (uint64_t j = 0; j < matrix_b->noCols; j++) {
            for (uint64_t k = 0; k < matrix_a->noCols; k++) {
                ord_result[i][j] += ord_a[i][k] * ord_b[k][j];
            }
        }
    }

    // Free up memory for the matrices in ordinary format
    for (uint64_t i = 0; i < matrix_a->noRows; i++) {
        free(ord_a[i]);
    }
    free(ord_a);

    for (uint64_t i = 0; i < matrix_b->noRows; i++) {
        free(ord_b[i]);
    }
    free(ord_b);

    // Convert the result to CSR format
    Matrix* matrix_tmp = ordinary_to_csr(matrix_a->noRows, matrix_b->noCols, ord_result);
    if (matrix_tmp == NULL) {
        for (uint64_t i = 0; i < matrix_a->noRows; i++)
        {
            free(ord_result[i]);
        }
        free(ord_result);

        errno = MATRIX_CONVERSION_ERROR;
        return;
    }

    // Copy the attributes of the temporary converted matrix
    matrix_result->noRows = matrix_tmp->noRows;
    matrix_result->noCols = matrix_tmp->noCols;
    matrix_result->colIndices = matrix_tmp->colIndices;
    matrix_result->rowPointers = matrix_tmp->rowPointers;
    matrix_result->rowPointersSize = matrix_tmp->rowPointersSize;
    matrix_result->values = matrix_tmp->values;
    matrix_result->valuesSize = matrix_tmp->valuesSize;

    // Free up memory of the result 2D array and temporary CSR matrix
    for (uint64_t i = 0; i < matrix_a->noRows; i++) {
        free(ord_result[i]);
    }
    free(matrix_tmp);
    free(ord_result);
}

void matr_mult_csr_V2(const void* a, const void* b, void* result) {
    // Sequential Gustavson's with no prediction
    Matrix* matrix_a = (Matrix*) a;
    Matrix* matrix_b = (Matrix*) b;
    Matrix* matrix_result = (Matrix*) result;

    // Check if multiplication is mathematically defined
    if (!can_multiply(matrix_a, matrix_b)) {
        errno = MATRIX_DIMENSION_ERROR;
        return;
    }

    // Initialize the subarrays with no size prediction
    if (init_empty_csr_matrix(matrix_a, matrix_b, matrix_result, NO_PREDICTION) == HEAP_MEMORY_ERROR) {
        errno = HEAP_MEMORY_ERROR;
        return;
    }

    multiply_V2(matrix_a, matrix_b, matrix_result);

    // Clean up the non zero values in the values array
    if (clean_up_csr(matrix_result) == HEAP_MEMORY_ERROR) {
        errno = HEAP_MEMORY_ERROR;
    }
}

void matr_mult_csr_V3(const void* a, const void* b, void* result) {
    // SIMD with SSE, no prediction
    Matrix* matrix_a = (Matrix*) a;
    Matrix* matrix_b = (Matrix*) b;
    Matrix* matrix_result = (Matrix*) result;

    // Check if multiplication is mathematically defined
    if (!can_multiply(matrix_a, matrix_b)) {
        errno = MATRIX_DIMENSION_ERROR;
        return;
    }

    // Initialize the subarrays with no size prediction
    if (init_empty_csr_matrix(matrix_a, matrix_b, matrix_result, NO_PREDICTION) == HEAP_MEMORY_ERROR) {
        errno = HEAP_MEMORY_ERROR;
        return;
    }

    if (multiply_V3(matrix_a, matrix_b, matrix_result) == HEAP_MEMORY_ERROR) {
        errno = HEAP_MEMORY_ERROR;
        return;
    }

    // Clean up the non zero values in the values array
    if (clean_up_csr(matrix_result) == HEAP_MEMORY_ERROR) {
        errno = HEAP_MEMORY_ERROR;
    }
}

void matr_mult_csr_V4(const void* a, const void* b, void* result) {
    // SIMD with AVX, no prediction
    if (!__builtin_cpu_supports("avx")) {
        // CPU doesn't support AVX, default to SSE implementation
        matr_mult_csr_V1(a, b, result);
        return;
    }

    Matrix* matrix_a = (Matrix*) a;
    Matrix* matrix_b = (Matrix*) b;
    Matrix* matrix_result = (Matrix*) result;

    // Check if multiplication is mathematically defined
    if (!can_multiply(matrix_a, matrix_b)) {
        errno = MATRIX_DIMENSION_ERROR;
        return;
    }

    // Initialize the subarrays with no size prediction
    if (init_empty_csr_matrix(matrix_a, matrix_b, matrix_result, NO_PREDICTION) == HEAP_MEMORY_ERROR) {
        errno = HEAP_MEMORY_ERROR;
        return;
    }

    if (multiply_V4(matrix_a, matrix_b, matrix_result) == HEAP_MEMORY_ERROR) {
        errno = HEAP_MEMORY_ERROR;
        return;
    }

    // Clean up the non zero values in the values array
    if (clean_up_csr(matrix_result) == HEAP_MEMORY_ERROR) {
        errno = HEAP_MEMORY_ERROR;
    }
}

void matr_mult_csr_V5(const void* a, const void* b, void* result) {
    // Gustavson with size prediction
    Matrix* matrix_a = (Matrix*) a;
    Matrix* matrix_b = (Matrix*) b;
    Matrix* matrix_result = (Matrix*) result;

    // Check if multiplication is mathematically defined
    if (!can_multiply(matrix_a, matrix_b)) {
        errno = MATRIX_DIMENSION_ERROR;
        return;
    }

    // Initialize the subarrays with values size prediction
    if (init_empty_csr_matrix(matrix_a, matrix_b, matrix_result, PREDICT_SIZE) == HEAP_MEMORY_ERROR) {
        errno = HEAP_MEMORY_ERROR;
        return;
    }

    multiply_V5(matrix_a, matrix_b, matrix_result);

    // Clean up the non zero values in the values array
    if (clean_up_csr(matrix_result) == HEAP_MEMORY_ERROR) {
        errno = HEAP_MEMORY_ERROR;
    }
}
