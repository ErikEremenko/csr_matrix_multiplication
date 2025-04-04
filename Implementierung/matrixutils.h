#ifndef MATRIXUTILS_H
#define MATRIXUTILS_H

// Default c library headers
#include <stdint.h>
#include <pthread.h>

// Our headers
#include "csrmatrix.h"
#include "utils.h"

// Define constants
#define THREAD_START_ERROR -5  // error starting threads (-5 to fit in with matrix.h error codes)

// Flags for size prediction in init_empty_csr_matrix()
#define PREDICT_SIZE 1
#define NO_PREDICTION 0


/*
The MultiplyArg struct is passed to multiply_main_implementation() as an argument. As a 
POSIX thread can only take in a single argument for the function, this struct is used to 
give the function the necessary arguments, such as the matrices and the row information.
*/
struct MultiplyArg {
    Matrix* matrix_a;
    Matrix* matrix_b;
    Matrix* matrix_result;
    uint64_t start_row;
    uint64_t end_row;
};

/*
This is the main implementation of the multiplication algorithm. It is called by each thread 
created in matr_mult_csr. It uses the same algorithm as version 2 (Gustavson's with no size 
prediction), but the function signature is changed to take in a single argument (a struct 
MultiplyArg), so that threads can call this function.
*/
void* multiply_main_implementation(void* void_arg);

/*
Converts a given Matrix (in CSR format) into a 2D-array of its values.

This function is used in version 1 of the matrix multiplication algorithm.

Return value: The matrix if memory was allocated, otherwise NULL.
*/
float** csr_to_ordinary(const Matrix* const csr);

/*
Converts a matrix in 2D-array format to the CSR format.

This function is used in version 1 of the matrix multiplication algorithm.

Return value: The matrix if memory was allocated, otherwise NULL.
*/
Matrix* ordinary_to_csr(const uint64_t noRows, const uint64_t noCols, float** ordinary);

/*
This is version 2 of the CSR Matrix multiplication algorithm, called in matr_mult_csr_V2().

The function uses Gustafson's algorithm without values size prediction.
*/
void multiply_V2(
    const Matrix* const restrict matrix_a, 
    const Matrix* const restrict matrix_b, 
    Matrix* const restrict matrix_result
    );

/*
This is version 3 of the CSR Matrix multiplication algorithm, called in matr_mult_csr_V3().

The function uses Gustafson's algorithm with SIMD, specifically 128-bit SSE registers.
However, unlike the main implementation, the size of the values array is not predicted,
as it does not fit in with the SIMD usage.

Return values:
    0 on success.
    HEAP_MEMORY_ERROR if the array used for SIMD cannot be malloc'ed.
*/
int multiply_V3(
    const Matrix* const restrict matrix_a, 
    const Matrix* const restrict matrix_b, 
    Matrix* const restrict matrix_result
    );

/*
This is version 4 of the CSR Matrix multiplication algorithm, called in matr_mult_csr_V4().

The function uses, similar to version 1, Gustafson's algorithm with SIMD. But the registers
used here are 256 bit AVX registers. The function defaults to version 1 if the computer
architecture does not support AVX. The size of the values array is not predicted.

Return values:
    0 on success.
    HEAP_MEMORY_ERROR if the array used for SIMD cannot be malloc'ed.
*/
int multiply_V4(
    const Matrix* const restrict matrix_a, 
    const Matrix* const restrict matrix_b, 
    Matrix* const restrict matrix_result
    );

/*
This is the version 5 of the CSR Matrix multiplication algorithm, called in matr_mult_csr_V5().

The function uses Gustafson's algorithm to multiply the matrices. The size of the values array 
is also predicted using an algorithm, drastically decreasing memory usage.
*/
void multiply_V5(
    const Matrix* const restrict matrix_a, 
    const Matrix* const restrict matrix_b, 
    Matrix* const restrict matrix_result
    );


/*
This function, called in matr_mult_csr(), starts threads that then multiply the rows
assigned to them. Each thread calls multiply_main_implementation().

Return values:
    0 if the threads all started with no error.
    HEAP_MEMORY_ERROR if memory for the thread array or the argument array could not be allocated.
*/
int start_threads(
    const unsigned int thread_count, pthread_t** threads, struct MultiplyArg*** arguments,
    Matrix* matrix_a, Matrix* matrix_b, Matrix* matrix_result
    );

/*
Checks if multiplication of the given matrices (in CSR format) is mathematically defined.

Return values:
    1 if they can be multiplied,
    0 if not.
*/
int can_multiply(const Matrix* const a, const Matrix* const b);

/*
Cleans up zero values in the given matrix, as the CSR format does not store them.

Return values:
    0 on success.
    HEAP_MEMORY_ERROR when a malloc/realloc call returns a null pointer.
*/
int clean_up_csr(Matrix* const matrix);

/*
This is a helper function is called in _clean_up_csr(). Used to clean up the zero values in the matrix.

Return values:
    0 on success.
    HEAP_MEMORY_ERROR when a malloc call returns a null pointer.
*/
int _clean_up_matrix_arrays(Matrix* matrix, uint64_t* non_zero_values);

/*
Predicts the amount of non-zero-values in the resulting CSR matrix.
Usage of this function greatly decreases memory usage, as the memory for the
array is dynamically allocated on the heap.

The prediction is stored in 'values_size'.

Return value: 
    0 on success.
    HEAP_MEMORY_ERROR if a malloc call fails.
*/
int predict_values_dimension(
    const Matrix* const restrict matrix_a, 
    const Matrix* const restrict matrix_b, 
    uint64_t* const values_size
    );

/*
Initializes an empty result CSR matrix.

Allocates memory for values, colIndices and rowPointers. This function assumes 
that memory for the Matrix struct itself was already allocated.

If the value passed for predict_flag is not 0, the size of the values array is 
mathematically predicted, in contrast to init_empty_csr_matrix_prediction(). 
This prediction is an upper limit to the size of values, which greatly improves 
performance in a very big but relatively sparse matrix (e.g 100x100 but only 1 nnz).

Return values:
    0 on success.
    HEAP_MEMORY_ERROR when memory for one of the subarrays could not be allocated.
*/
int init_empty_csr_matrix(
    const Matrix* const restrict matrix_a, const Matrix* const restrict matrix_b, 
    Matrix* const restrict result, const int predict_flag
    );

/*
Tests if the two given CSR matrices are equal, regardless of the ordering of
values and column indices. The values in the matrices are ascendingly sorted
according to their column indices.

This function may be important for tutors if they wish to check the validity
of the multiplication algorithm.

Return value: 1 if equal, otherwise 0.
*/
int equals(Matrix* const restrict matrix_a, Matrix* const restrict matrix_b);

/*
Sorts (ascending) the values and column indices of the given CSR matrix.

Uses insertion sort as the sorting algorithm.
*/
void sort_matrix(Matrix* const matrix, const uint64_t beg, const uint64_t end);

/*
Prints an ordinary matrix (2D-array of values) on the console.
*/
void print_ordinary_matrix(float** matrix, const uint64_t noRows, const uint64_t noCols);

/*
Prints a matrix in CSR format on the console.
*/
void print_csr_matrix(const Matrix* const csr);

/*
Swaps the column indices on the given indices.

This function is called in sort_matrix().
*/
void _swap_cols(uint64_t* const cols, const uint64_t size, const uint64_t i, const uint64_t j);

/*
Swaps values on the given indices.

This function is called in sort_matrix().
*/
void _swap_values(float* const values, const uint64_t size, const uint64_t i, const uint64_t j);

#endif
