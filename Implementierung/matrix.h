#ifndef MATRIX_H
#define MATRIX_H

#include <stdint.h>
#include <stddef.h>

#include "csrmatrix.h"

// Define constants
#define MATRIX_DIMENSION_ERROR -2  // matrices cannot be multiplied
#define MATRIX_CONVERSION_ERROR -3  // error converting CSR to 2D-array in V4

#define MIN_THREADS 2
#define THREAD_COND_MIN_VALUES 40
#define THREAD_COND_MIN_ROWS 4

/*
Implementation V0, Multithreading (Hauptimplementierung)

The main implementation uses multithreading when it's favourable. Otherwise,
matr_mult_csr_V5() is called, which uses size prediction.

If an error occurs while initializing the subarrays, they are all initialized
to be a NULL pointer. Thus, free() and free_csr_matrix() can be called with no
worries about double freeing a pointer.

Sets errno to:
    MATRIX_DIMENSION_ERROR if the multiplication is mathematically not defined.
    HEAP_MEMORY_ERROR if the matrix cannot be malloc'ed.
    THREAD_START_ERROR if one of the threads cannot be created/started.
*/
void matr_mult_csr(const void* a, const void* b, void* result);

/*
Implementation V1 converts the matrices into a 2D array and then multiplies
them using standard matrix multiplication. The resulting 2D array is then
converted back to a CSR matrix.

errno should be set to 0 before calling this function in order to check if
the matrix was succesfully multiplied.

If an error occurs while initializing the subarrays, they are all initialized
to be a NULL pointer. Thus, free() and free_csr_matrix() can be called with no
worries about double freeing a pointer.

Sets errno to:
    MATRIX_DIMENSION_ERROR if the multiplication is mathematically not defined.
    MATRIX_CONVERSION_ERROR if a conversion between CSR and 2D array form fails.
    HEAP_MEMORY_ERROR if a malloc call fails.
*/
void matr_mult_csr_V1(const void* a, const void* b, void* result);

/*
Implementation V2 uses Gustavson's algorithm to multiply the CSR matrices.
The size of the values array is not predicted, which may lead to a unusually
high memory usage.

If an error occurs while initializing the subarrays, they are all initialized
to be a NULL pointer. Thus, free() and free_csr_matrix() can be called with no
worries about double freeing a pointer.

Sets errno to:
    MATRIX_DIMENSION_ERROR if the multiplication is mathematically not defined.
    HEAP_MEMORY_ERROR if the matrix cannot be malloc'ed.
*/
void matr_mult_csr_V2(const void* a, const void* b, void* result);

/*
V3 implements SIMD using 128 bit SSE registers. The size of the values
array is not predicted.

If an error occurs while initializing the subarrays, they are all initialized
to be a NULL pointer. Thus, free() and free_csr_matrix() can be called with no
worries about double freeing a pointer.

Sets errno to:
    MATRIX_DIMENSION_ERROR if the multiplication is mathematically not defined.
    HEAP_MEMORY_ERROR if the matrix cannot be malloc'ed.
*/
void matr_mult_csr_V3(const void* a, const void* b, void* result);

/*
SIMD with AVX

V4, similar to V3, implements a SIMD approach using 256-bit AVX registers.
The size of the values array is not predicted.

If an error occurs while initializing the subarrays, they are all initialized
to be a NULL pointer. Thus, free() and free_csr_matrix() can be called with no
worries about double freeing a pointer.

Sets errno to:
    MATRIX_DIMENSION_ERROR if the multiplication is mathematically not defined.
    HEAP_MEMORY_ERROR if the matrix cannot be malloc'ed.
*/
void matr_mult_csr_V4(const void* a, const void* b, void* result);

/*
V5 makes use of Gustavson's algorithm but also predicts the size of the values
array in the CSR format. This may lead to a longer runtime but the memory
usage is drastically lower in comparison to other implementations.

errno should be set to 0 before calling this function in order to check if
the matrix was succesfully multiplied.

If an error occurs while initializing the subarrays, they are all initialized
to be a NULL pointer. Thus, free() and free_csr_matrix() can be called with no
worries about double freeing a pointer.

Sets errno to:
    MATRIX_DIMENSION_ERROR if the multiplication is mathematically not defined.
    HEAP_MEMORY_ERROR if the matrix cannot be malloc'ed.
*/
void matr_mult_csr_V5(const void* a, const void* b, void* result);

#endif
