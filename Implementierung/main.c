// We need to define this for CLOCK_MONOTONIC
#define _POSIX_C_SOURCE 200809L

// Default C Library headers
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include <errno.h>

// Our header files
#include "constants.h"
#include "utils.h"
#include "csrmatrix.h"
#include "matrixutils.h"
#include "matrix.h"


// Multiplication algorithm function type
typedef void (*mult_fn)(const void* a, const void* b, void* result);

/*
Function that returns the matrix multiplication algorithm the user wants to use.
*/
mult_fn choose_mult_fn(uint8_t implementation) {
    switch (implementation) {
        case 0:
            return matr_mult_csr;
        case 1:
            return matr_mult_csr_V1;
        case 2:
            return matr_mult_csr_V2;
        case 3:
            return matr_mult_csr_V3;
        case 4:
            return matr_mult_csr_V4;
        case 5:
            return matr_mult_csr_V5;
        default:  // invalid input (this is actually never the case)
            return NULL;
    }
}

int main(int argc, char** argv) {
    // Args that must be provided
    char* filename_matrix_a = NULL;
    char* filename_matrix_b = NULL;
    char* filename_matrix_output = NULL;

    // Optional args
    uint8_t implementation = 0;  // which implementation to use
    int measure_flag = 0;  // flag to measure execution time
    uint64_t number_measures = 1;  // how many times we want to execute the function

    // Used in the main function for printing errors
    char* error_message = NULL;

    int parse_result = parse_arguments(
            argc, argv,
            &filename_matrix_a, &filename_matrix_b, &filename_matrix_output,
            &implementation, &measure_flag, &number_measures,
            &error_message
        );

    switch (parse_result) {
        case ARGPARSE_ERROR:
            main_error:
            // Print error message
            if (error_message == NULL) {  // if there was an error creating the error message
                fprintf(stderr, "%s", HEAP_MEMORY_ERROR_MSG);
            } else {
                fprintf(stderr, "%s", error_message);
                // Free error message memory (heap-allocated string)
                free(error_message);
            }
            // Free filename strings (they are either null ptrs or valid values)
            // Matrices are already freed before reaching here
            free_pointers(3, filename_matrix_a, filename_matrix_b, filename_matrix_output);
            // Print that we are exiting due to an error
            fprintf(stderr, "%s", EXIT_FAIL_MSG);
            // Return failure as specified in stdlib.h
            return EXIT_FAILURE;
        case ARGPARSE_SUCCESS:
            // Read matrix a
            Matrix* matrix_a;
            int op_result = read_matrix_from_file(filename_matrix_a, &matrix_a);
            switch (op_result) {
                case FILE_OPEN_ERROR:
                    set_error_message(&error_message, FILE_OPEN_ERROR_MSG, filename_matrix_a);
                    goto main_error;
                case MATRIX_FILE_FORMAT_ERROR:
                    set_error_message(&error_message, MATRIX_FILE_FORMAT_ERROR_MSG, filename_matrix_a);
                    goto main_error;
                case HEAP_MEMORY_ERROR:
                    set_error_message(&error_message, HEAP_MEMORY_ERROR_MSG);
                    goto main_error;
                default:  // MATRIX_READ_SUCCESS
                    break;
            }

            // Read matrix b
            Matrix* matrix_b;
            op_result = read_matrix_from_file(filename_matrix_b, &matrix_b);
            switch (op_result) {
                case FILE_OPEN_ERROR:
                    // Error reading matrix b, we only have to free matrix a
                    free_csr_matrix(matrix_a);
                    set_error_message(&error_message, FILE_OPEN_ERROR_MSG, filename_matrix_b);
                    goto main_error;
                case MATRIX_FILE_FORMAT_ERROR:
                    free_csr_matrix(matrix_a);
                    set_error_message(&error_message, MATRIX_FILE_FORMAT_ERROR_MSG, filename_matrix_b);
                    goto main_error;
                case HEAP_MEMORY_ERROR:
                    free_csr_matrix(matrix_a);
                    set_error_message(&error_message, HEAP_MEMORY_ERROR_MSG);
                    goto main_error;
                default:  // MATRIX_READ_SUCCESS
                    break;
            }

            // Initialize the result matrix with NULL ptrs (safety for free())
            Matrix* matrix_result = malloc(sizeof(Matrix));
            if (matrix_result == NULL) {
                // Error allocating memory, free a and b
                free_csr_matrices(2, matrix_a, matrix_b);
                set_error_message(&error_message, HEAP_MEMORY_ERROR_MSG);
                goto main_error;
            }
            matrix_result->values = NULL;
            matrix_result->colIndices = NULL;
            matrix_result->rowPointers = NULL;

            // Get implementation/matrix multiplication algorithm the user wants
            mult_fn matr_mult_csr_fn = choose_mult_fn(implementation);
            if (matr_mult_csr_fn == NULL) {
                goto main_error;
            }
            
            if (measure_flag) {
                // Create temporary result matrix here to avoid time measurement
                Matrix* tmp_result = malloc(sizeof(Matrix));
                if (tmp_result == NULL) {
                    // Error allocating memory to result matrix, free a and b
                    free_csr_matrices(2, matrix_a, matrix_b);
                    free(matrix_result);
                    set_error_message(&error_message, HEAP_MEMORY_ERROR_MSG);
                    goto main_error;
                }
                tmp_result->values = NULL;
                tmp_result->colIndices = NULL;
                tmp_result->rowPointers = NULL;

                // Get start time
                struct timespec start;
                clock_gettime(CLOCK_MONOTONIC, &start);

                // Do the first iteration outside of the loop to store results
                errno = 0;
                matr_mult_csr_fn(matrix_a, matrix_b, matrix_result);
                switch (errno) {
                    case MATRIX_CONVERSION_ERROR:
                        set_error_message(&error_message, MATRIX_CONV_ERROR_MSG);
                        free_csr_matrices(3, matrix_a, matrix_b, matrix_result);
                        free(tmp_result);
                        goto main_error;
                    case MATRIX_DIMENSION_ERROR:
                        set_error_message(
                            &error_message, MATRIX_DIM_ERROR_MSG,
                            matrix_a->noRows, matrix_a->noCols,
                            matrix_b->noRows, matrix_b->noCols
                            );
                        free_csr_matrices(3, matrix_a, matrix_b, matrix_result);
                        free(tmp_result);
                        goto main_error;
                    case THREAD_START_ERROR:
                        set_error_message(&error_message, THREAD_START_ERROR_MSG);
                        free_csr_matrices(3, matrix_a, matrix_b, matrix_result);
                        free(tmp_result);
                        goto main_error;
                    case HEAP_MEMORY_ERROR:
                        set_error_message(&error_message, HEAP_MEMORY_ERROR_MSG);
                        free_csr_matrices(3, matrix_a, matrix_b, matrix_result);
                        free(tmp_result);
                        goto main_error;
                    default:  // no error, go on
                        break;
                }

                // Do the matrix multiplication
                for (uint32_t i = 1; i < number_measures; i++) {
                    errno = 0;
                    matr_mult_csr_fn(matrix_a, matrix_b, tmp_result);
                    switch (errno) {
                        case MATRIX_CONVERSION_ERROR:
                            set_error_message(&error_message, MATRIX_CONV_ERROR_MSG);
                            free_csr_matrices(3, matrix_a, matrix_b, matrix_result);
                            free(tmp_result);
                            goto main_error;
                        case MATRIX_DIMENSION_ERROR:
                            set_error_message(
                                &error_message, MATRIX_DIM_ERROR_MSG,
                                matrix_a->noRows, matrix_a->noCols,
                                matrix_b->noRows, matrix_b->noCols
                                );
                            free_csr_matrices(3, matrix_a, matrix_b, matrix_result);
                            free(tmp_result);
                            goto main_error;
                        case THREAD_START_ERROR:
                            set_error_message(&error_message, THREAD_START_ERROR_MSG);
                            free_csr_matrices(3, matrix_a, matrix_b, matrix_result);
                            free(tmp_result);
                            goto main_error;
                        case HEAP_MEMORY_ERROR:
                            set_error_message(&error_message, HEAP_MEMORY_ERROR_MSG);
                            free_csr_matrices(3, matrix_a, matrix_b, matrix_result);
                            free(tmp_result);
                            goto main_error;
                        default:  // no error, go on
                            break;
                    }
                    // Free subarrays of the temporary result matrix
                    free(tmp_result->values);
                    free(tmp_result->colIndices);
                    free(tmp_result->rowPointers);
                }

                // Calculate time it took for the function to execute number_measures times
                struct timespec end;
                clock_gettime(CLOCK_MONOTONIC, &end);
                double time = end.tv_sec - start.tv_sec + 1e-9 * (end.tv_nsec - start.tv_nsec);

                free(tmp_result);  // free the temporary result matrix for good

                // Print time measurement on the console
                printf("Took %g seconds to multiply\n", time);
            } else {
                // Just do the multiplication, no time measurement
                matr_mult_csr_fn(matrix_a, matrix_b, matrix_result);
                // Check if multiplication was successful
                switch (errno) {
                    case MATRIX_CONVERSION_ERROR:
                        set_error_message(&error_message, MATRIX_CONV_ERROR_MSG);
                        free_csr_matrices(2, matrix_a, matrix_b);
                        free(matrix_result);  // the subarrays are already NULL ptrs, no need to free them
                        goto main_error;
                    case MATRIX_DIMENSION_ERROR:
                        set_error_message(
                            &error_message, MATRIX_DIM_ERROR_MSG,
                            matrix_a->noRows, matrix_a->noCols,
                            matrix_b->noRows, matrix_b->noCols
                            );
                        free_csr_matrices(2, matrix_a, matrix_b);
                        free(matrix_result);
                        goto main_error;
                    case THREAD_START_ERROR:
                        set_error_message(&error_message, THREAD_START_ERROR_MSG);
                        free_csr_matrices(2, matrix_a, matrix_b);
                        free(matrix_result);
                        goto main_error;
                    case HEAP_MEMORY_ERROR:
                        set_error_message(&error_message, HEAP_MEMORY_ERROR_MSG);
                        free_csr_matrices(2, matrix_a, matrix_b);
                        free(matrix_result);
                        goto main_error;
                    default:  // no error, go on
                        break;
                }
            }

            // Write result to file
            op_result = write_matrix_to_file(filename_matrix_output, matrix_result);
            switch (op_result) {
                case FILE_OPEN_ERROR:
                    // Error writing result, free all matrices
                    free_csr_matrices(3, matrix_a, matrix_b, matrix_result);
                    set_error_message(&error_message, FILE_OPEN_ERROR_MSG, filename_matrix_output);
                    goto main_error;
                case FILE_WRITE_ERROR:
                    free_csr_matrices(3, matrix_a, matrix_b, matrix_result);
                    set_error_message(&error_message, FILE_WRITE_ERROR_MSG, filename_matrix_output);
                    goto main_error;
                default:  // MATRIX_WRITE_SUCCESS
                    break;
            }

            // Free filename string pointers
            free_pointers(3, filename_matrix_a, filename_matrix_b, filename_matrix_output);
            // Free matrices
            free_csr_matrices(3, matrix_a, matrix_b, matrix_result);
            return EXIT_SUCCESS;
        case ARGPARSE_HELP:
            // Display help message and exit
            free_pointers(3, filename_matrix_a, filename_matrix_b, filename_matrix_output);
            printf("%s", HELP_MSG);
            return EXIT_SUCCESS;
        default:
            goto main_error;
    }

    return EXIT_SUCCESS;
}
