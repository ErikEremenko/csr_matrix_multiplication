#ifndef UTILS_H
#define UTILS_H

// Default C Library headers
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

// Our header files
#include "constants.h"
#include "csrmatrix.h"


// Message constants
#define EXIT_FAIL_MSG "Exiting due to failure...\n"
extern const char* HELP_MSG;  // message to print when -h is passed as an argument
extern const char* HOW_TO_USE_MSG;  //  message to tell the user to use -h or --help
extern const char* ILLEGAL_NUMBER_MEASURES_MSG;  // message to print when a non-positive number of measures was specified
extern const char* ILLEGAL_ARG_MSG;  // message to tell the user that they've passed an illegal argument (like -x)
extern const char* MISSING_ARG_MSG;  // message to print when there is an option missing
extern const char* ILLEGAL_IMPLEMENTATION_MSG;  // message to tell the user the implementation they passed is illegal
extern const char* NON_OPTION_ARGS_MSG;  // message to print when a nonoption argument is given (doesn't have newline at the end)
extern const char* MISSING_FILENAME_A_MSG;  // message to print when filename a is missing
extern const char* MISSING_FILENAME_B_MSG;  // message to print when filename b is missing
extern const char* MISSING_FILENAME_O_MSG;  // message to print when the output filename is missing
extern const char* ALREADY_PARSED_MSG;  // message to print when an argument is given twice (like -a -a)

extern const char* FILE_OPEN_ERROR_MSG;  // message displayed when there is an error opening a given file
extern const char* MATRIX_FILE_FORMAT_ERROR_MSG;  // message to print when the matrix in the file is not correctly formatted
extern const char* FILE_WRITE_ERROR_MSG;  // message displayed when there is an error writing to a given file

extern const char* MATRIX_DIM_ERROR_MSG;  // message displayed when matrix multiplication is undefined
extern const char* MATRIX_CONV_ERROR_MSG;  // message displayed when CSR/2D array conversion fails
extern const char* THREAD_START_ERROR_MSG;  // message displayed when multithreading algorithm fails to start a thread

// Error codes start from -2 because -1 is HEAP_MEMORY_ERROR, which is used throughout the whole program
// Argument parsing error codes
#define ARGPARSE_ERROR -2
#define ARGPARSE_SUCCESS 0
#define ARGPARSE_HELP 1

// Matrix writing error codes
#define ARRAY_WRITE_SUCCESS 0
#define MATRIX_WRITE_SUCCESS 0
#define INCOMPLETE_MATRIX_FILE_ERROR -2
#define MATRIX_FILE_FORMAT_ERROR -3

// File I/O error codes
#define FILE_OPEN_ERROR -2
#define FILE_WRITE_ERROR -3

// Matrix read success code
#define MATRIX_READ_SUCCESS 0
// Passed to matrix reading function when the line to be read is not the last line in the file.
#define NOT_LAST_LINE '\n'


/*
Parses command line arguments and determines which implementation should be used,
from which files the matrices should be read, in which file the resulting matrix
should be stored, whether the time it takes for the program to be executed should
be measured, and also how many times the time should be measured.

Return values:
    ARGPARSE_SUCCESS when all arguments have been parsed correctly.
    ARGPARSE_ERROR when there was an error in parsing, exit prematurely.
    ARGPARSE_HELP when the help message should be displayed, end program.
*/
int parse_arguments(
    int argc,
    char** argv,
    char** filename_matrix_a,
    char** filename_matrix_b,
    char** filename_matrix_output,
    uint8_t* implementation,
    int* measure_flag,
    uint64_t* number_measures,
    char** error_message
);

/*
Reads a matrix in CSR format from the given filename. Allocates memory on the heap for the matrix and
its data (values, column indices, row pointers). This matrix should then be free'd.

On error, the matrix is not allocated.

Return values:
    MATRIX_READ_SUCCESS when the matrix is successfully read.
    FILE_OPEN_ERROR when the file cannot be opened.
    MATRIX_FILE_FORMAT_ERROR if the file is improperly formatted (invalid CSR).
    HEAP_MEMORY_ERROR if there is an error allocating memory to the matrix.
*/
int read_matrix_from_file(const char* filename, Matrix** matrix);

/*
Writes a non-NULL CSR matrix to the given filename. Does NOT free the matrix or its subarrays.

Return values:
    MATRIX_WRITE_SUCCESS when the matrix is successfully written to the file.
    FILE_OPEN_ERROR when the file cannot be opened.
    FILE_WRITE_ERROR when there is an error writing to the file.
*/
int write_matrix_to_file(const char* filename, const Matrix* const matrix);

/*
This function is used to set the error message to be printed on stderr later on.
'error_message' will hold the formatted string at the end of the function call.
'format' is the error message string to be formatted. The function is then variadic, 
it takes in the values that should be inserted in the format.

Error_message is set to NULL when an error occurs allocating memory to it. This can
also be used to check for errors, instead of using the return value of the function.

Return values:
    0 on success. 
    HEAP_MEMORY_ERROR when error_message cannot be malloc'ed.
*/
int set_error_message(char** error_message, const char* format, ...);

/*
Reads a filename string using getopt_long().

This function is called in parse_arguments() and should not be called outside of it.

Return values:
    0 on success.
    HEAP_MEMORY_ERROR if the memory for the string cannot be allocated.
*/
int _parse_filename(char** dest, const char* const src);


/*
Reads a unsigned long long (uint64_t) array from a file, whose values are separated by commas. 
Allocates memory for the array, stores its length in array_size. This memory should then be free'd.
If an error occurs during parsing, the array is freed up before returning.

This function is called in read_matrix_from_file() and should generally not be called outside of it.

eof_flag should be set to EOF when parsing the last line in a file, otherwise NOT_LAST_LINE.

Return values:
    MATRIX_READ_SUCCESS when the matrix was successfully read.
    HEAP_MEMORY_ERROR when memory for the array could not be allocated.
    MATRIX_FILE_FORMAT_ERROR when the file is not correctly formatted.
*/
int _read_uint_64_array(FILE* file, uint64_t** array, uint64_t* array_size, const int eof_flag);

/*
Reads a float array from a file, whose values are separated by commas. Allocates
memory for the array, stores its length in array_size. This memory should then be free'd.
If an error occurs during parsing, the array is freed up before returning.

This function is called in read_matrix_from_file() and should generally not be called outside of it.

eof_flag should be set to EOF when parsing the last line in a file, otherwise NOT_LAST_LINE.

Return values:
    MATRIX_READ_SUCCESS when the matrix was successfully read.
    HEAP_MEMORY_ERROR when memory for the array could not be allocated.
    MATRIX_FILE_FORMAT_ERROR when the file is not correctly formatted.
*/
int _read_float_array(FILE* file, float** array, uint64_t* array_size, const int eof_flag);

/*
Called in read_matrix_to_file() to check if the values parsed from the file actually make up
a valid CSR matrix. The values must be non zero and the size of the parsed array must be less
than or equal to what the matrix can mathematically hold. This check is important because the 
multiplication functions assume that the CSR matrices given have a valid structure.

On error, the value array is freed.

This function is called in read_matrix_from_file() and should not be called outside of it.

Return values:
    0 if the row pointers are properly structured.
    MATRIX_FILE_FORMAT_ERROR if they are not.
*/
int _check_values(
    float* values, const uint64_t values_size,
    const uint64_t noRows, const uint64_t noCols
);

/*
Called in read_matrix_to_file() to check if the column index parsed from the file actually make up
a valid CSR matrix. This is important because the multiplication functions assume that the given
CSR matrices have a valid structure.

On error, the column index array is freed.

This function is called in read_matrix_from_file() and should not be called outside of it.

Return values:
    0 if the row pointers are properly structured.
    MATRIX_FILE_FORMAT_ERROR if they are not.
*/
int _check_col_indices (
    uint64_t* colIndices, const uint64_t col_indices_size, 
    const uint64_t values_size, const uint64_t noCols
    );

/*
Called in read_matrix_to_file() to check if the row pointers parsed from the file actually make up
a valid CSR matrix. This is important because the multiplication functions assume that the given
CSR matrices have a valid structure.

On error, the row pointer array is freed.

This function is called in read_matrix_from_file() and should not be called outside of it.

Return values:
    0 if the row pointers are properly structured.
    MATRIX_FILE_FORMAT_ERROR if they are not.
*/
int _check_row_pointers(
    uint64_t* rowPointers, const uint64_t row_pointers_size, 
    const uint64_t values_size, const uint64_t noRows, const uint64_t noCols
    );

/*
Writes an array of floats to the given file. The values are each separated by commas.

This function is called in write_matrix_to_file() and should not be called outside of it.

Return values:
    ARRAY_WRITE_SUCCESS on success.
    FILE_WRITE_ERROR when the array was not successfully written to the file.
*/
int _write_float_array(FILE* file, float* const array, size_t size);
/*
Writes an array of uint64_t's to the given file. The values are each separated by commas.

This function is called in write_matrix_to_file() and should not be called outside of it.

Return values:
    ARRAY_WRITE_SUCCESS on success.
    FILE_WRITE_ERROR when the array was not successfully written to the file.
*/
int _write_uint_64_array(FILE* file, uint64_t* const array, size_t size);


/*
Frees a heap-allocated CSR matrix and all of its heap-allocated attributes.
*/
void free_csr_matrix(Matrix* matrix);

/*
A function that calls free() on all given pointers.

Used to make the main function clearer and easier to understand.
*/
void free_pointers(const int n, ...);

/*
A function that calls free_csr_matrix() on all given matrices.

Used to make the main function clearer and easier to understand.
*/
void free_csr_matrices(const int n, ...);

/*
Checks for a possible overflow during a memory allocation.

Return values:
    The malloc'ed pointer on success.
    NULL on overflow or if malloc returns NULL.
*/
static inline void* malloc_safe(uint64_t element_size, uint64_t num_elements) {
    // Check for overflow
    uint64_t alloc_size;
    if (__builtin_umull_overflow(element_size, num_elements, &alloc_size)) {
        return NULL;
    }
    return malloc(alloc_size);
}

/*
Checks for a possible overflow during a realloc call.
If there is an overflow or if realloc returns NULL, the original pointer
passed to the function is freed.

Return values:
    The realloc'ed pointer on success.
    NULL on overflow or if realloc returns NULL.
*/
static inline void* realloc_safe(void* const ptr, uint64_t element_size, uint64_t num_elements) {
    // Check for overflow
    uint64_t alloc_size;
    if (__builtin_umull_overflow(element_size, num_elements, &alloc_size)) {
        // overflow, free original pointer
        free(ptr);
        return NULL;
    }
    
    void* tmp = ptr;
    tmp = realloc(tmp, alloc_size);
    if (tmp == NULL) {
        // realloc failed, free original pointer
        free(ptr);
        return NULL;
    }
    return tmp;  // return the new realloc'ed pointer
}

#endif
