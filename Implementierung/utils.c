// C library headers
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <stdarg.h>

// Our headers
#include "csrmatrix.h"
#include "utils.h"

// Constants for argument parsing
const char* matrix_optstring = ":a:b:o:hB::V:";
const struct option matrix_options[] = {
        {"help", no_argument, NULL, 'h'},
        {0, 0, 0, 0}
    };

// Error message constants
const char* HELP_MSG =
"Required arguments:\n"
"  -a <filename>    Input file containing matrix A\n"
"  -b <filename>    Input file containing matrix B\n"
"  -o <filename>    Output file containing resulting matrix\n"
"\n"
"Optional arguments:\n"
"  -h, --help     Display help message and exit program\n"
"  -B<n>    Measure time n times and print (n is optional, default: n = 1)\n"
"  -V <n>    Implementation to use (default: n = 0)\n";

const char* HOW_TO_USE_MSG = "Add -h or --help to learn how to use the program.\n";
const char* ILLEGAL_NUMBER_MEASURES_MSG = "Number of times to measure cannot be \"%s\"\n";
const char* ILLEGAL_ARG_MSG = "Argument '-%c' invalid\n";
const char* MISSING_ARG_MSG = "Missing argument for option -%c\n";
const char* ILLEGAL_IMPLEMENTATION_MSG = "The implementation to use cannot be \"%s\"\n";
const char* NON_OPTION_ARGS_MSG = "Non-option arguments are not allowed: %s\n";
const char* ALREADY_PARSED_MSG = "Argument '-%c' was given twice\n";

const char* MISSING_FILENAME_A_MSG = "Missing filename for matrix A\n";
const char* MISSING_FILENAME_B_MSG = "Missing filename for matrix B\n";
const char* MISSING_FILENAME_O_MSG = "Missing filename for the resulting matrix\n";

const char* FILE_OPEN_ERROR_MSG = "Could not open file \"%s\"\n";
const char* MATRIX_FILE_FORMAT_ERROR_MSG = "File \"%s\" is not correctly formatted\n";
const char* FILE_WRITE_ERROR_MSG = "Error writing to file \"%s\"\n";

const char* MATRIX_DIM_ERROR_MSG = "Incompatible matrix dimensions: %lux%lu and %lux%lu\n";
const char* MATRIX_CONV_ERROR_MSG = "Error converting matrix (CSR-2D array)\n";
const char* THREAD_START_ERROR_MSG = "Error starting a thread during matrix multiplication\n";


// MAIN UTILITY FUNCTIONS BELOW //
// ---------------------------- //


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
) {
    opterr = 0;  // silence error messages from getopt

    //                   a  b  o  B  V
    int flag_array[5] = {0, 0, 0, 0, 0};

    int ch;
    char* endptr;  // used in string to number conversion
    while ((ch = getopt_long(argc, argv, matrix_optstring, matrix_options, NULL)) != -1) {
        switch (ch) {
        case 'h':  // -h or --help
            return ARGPARSE_HELP;
        case 'a':
            if (flag_array[0]) {
                set_error_message(error_message, ALREADY_PARSED_MSG, 'a');
                return ARGPARSE_ERROR;
            }
            flag_array[0] = 1;

            if (_parse_filename(filename_matrix_a, optarg) == HEAP_MEMORY_ERROR) {
                set_error_message(error_message, HEAP_MEMORY_ERROR_MSG);
                return ARGPARSE_ERROR;
            }
            break;
        case 'b':
            if (flag_array[1]) {
                set_error_message(error_message, ALREADY_PARSED_MSG, 'b');
                return ARGPARSE_ERROR;
            }
            flag_array[1] = 1;

            if (_parse_filename(filename_matrix_b, optarg) == HEAP_MEMORY_ERROR) {
                set_error_message(error_message, HEAP_MEMORY_ERROR_MSG);
                return ARGPARSE_ERROR;
            }
            break;
        case 'o':
            if (flag_array[2]) {
                set_error_message(error_message, ALREADY_PARSED_MSG, 'o');
                return ARGPARSE_ERROR;
            }
            flag_array[2] = 1;

            if (_parse_filename(filename_matrix_output, optarg) == HEAP_MEMORY_ERROR) {
                set_error_message(error_message, HEAP_MEMORY_ERROR_MSG);
                return ARGPARSE_ERROR;
            }
            break;
        case 'B':
            if (flag_array[3]) {
                set_error_message(error_message, ALREADY_PARSED_MSG, 'B');
                return ARGPARSE_ERROR;
            }
            flag_array[3] = 1;

                *measure_flag = 1;
                if (optarg == NULL) {
                    // optional argument not given, default to 1
                    *number_measures = 1;
                } else {
                    // Check if a negative number was given
                    if (optarg[0] == '-') {
                        set_error_message(error_message, ILLEGAL_NUMBER_MEASURES_MSG, optarg);
                        return ARGPARSE_ERROR;
                    }

                    // Convert from string to uint64_t
                    errno = 0;
                    *number_measures = strtoull(optarg, &endptr, 10);
                    if (errno || *endptr != '\0' || *number_measures == 0) {  // error converting
                        set_error_message(error_message, ILLEGAL_NUMBER_MEASURES_MSG, optarg);
                        return ARGPARSE_ERROR;
                    }
                }
                break;
            case 'V':
                if (flag_array[4]) {
                    set_error_message(error_message, ALREADY_PARSED_MSG, 'V');
                    return ARGPARSE_ERROR;
                }
                flag_array[4] = 1;

                // Check if a negative number was given
                if (optarg[0] == '-') {
                    set_error_message(error_message, ILLEGAL_IMPLEMENTATION_MSG, optarg);
                    return ARGPARSE_ERROR;
                }

                // Convert from string to uint64_t
                errno = 0;
                uint64_t tmp_implementation = strtoull(optarg, &endptr, 10);
                if (errno || *endptr != '\0') {
                    set_error_message(error_message, ILLEGAL_IMPLEMENTATION_MSG, optarg);
                    return ARGPARSE_ERROR;
                }
                
                // Check if it is a valid implementation number
                if (tmp_implementation > 255 || tmp_implementation >= NUMBER_OF_IMPLEMENTATIONS) {
                    set_error_message(error_message, ILLEGAL_IMPLEMENTATION_MSG, optarg);
                    return ARGPARSE_ERROR;
                }
                *implementation = (uint8_t) tmp_implementation;
                break;
            case '?':  // not a valid argument
                set_error_message(error_message, ILLEGAL_ARG_MSG, optopt);
                return ARGPARSE_ERROR;
            case ':':  // missing option
                set_error_message(error_message, MISSING_ARG_MSG, optopt);
                return ARGPARSE_ERROR;
            default:
                break;
        }
    }

    // Check if -a -b -o were given
    if (*filename_matrix_a == NULL) {
        set_error_message(error_message, MISSING_FILENAME_A_MSG);
        return ARGPARSE_ERROR;
    }

    if (*filename_matrix_b == NULL) {
        set_error_message(error_message, MISSING_FILENAME_B_MSG);
        return ARGPARSE_ERROR;
    }

    if (*filename_matrix_output == NULL) {
        set_error_message(error_message, MISSING_FILENAME_O_MSG);
        return ARGPARSE_ERROR;
    }

    // Check if there are any nonoption args
    if (optind < argc) {
        set_error_message(error_message, NON_OPTION_ARGS_MSG, argv[optind]);
        return ARGPARSE_ERROR;
    }

    return ARGPARSE_SUCCESS;
}

int read_matrix_from_file(const char* filename, Matrix** matrix) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        return FILE_OPEN_ERROR;
    }

    // Read noRows and noCols
    uint64_t* row_col_array;
    uint64_t row_col_array_size = 0;
    int read_result = _read_uint_64_array(file, &row_col_array, &row_col_array_size, NOT_LAST_LINE);
    switch (read_result) {
        case HEAP_MEMORY_ERROR:
        case MATRIX_FILE_FORMAT_ERROR:
            // error parsing no rows and cols, array already freed inside the function
            fclose(file);
            return read_result;
        case MATRIX_READ_SUCCESS:
        default:
            break;
    }
    
    // Check if the first row had exactly two values
    if (row_col_array_size != 2) {
        free(row_col_array);
        fclose(file);
        return MATRIX_FILE_FORMAT_ERROR;
    }

    // Store the arrays values in variables and free up its memory
    uint64_t noRows = row_col_array[0];
    uint64_t noCols = row_col_array[1];
    free(row_col_array);

    // Check if noRows and noCols are valid (> 0)
    if (!noRows || !noCols) {
        fclose(file);
        return MATRIX_FILE_FORMAT_ERROR;
    }

    // Read values
    float* values;
    uint64_t values_size = 0;
    read_result = _read_float_array(file, &values, &values_size, NOT_LAST_LINE);
    switch (read_result) {
        case HEAP_MEMORY_ERROR:
        case MATRIX_FILE_FORMAT_ERROR:
            // error parsing values, array already freed in function
            fclose(file);
            return read_result;
        case MATRIX_READ_SUCCESS:
        default:
            break;
    }

    // Check if values size is more than possible (rows x cols) and check for zero values
    if (_check_values(values, values_size, noRows, noCols) == MATRIX_FILE_FORMAT_ERROR) {
        // array freed on error in function
        fclose(file);
        return MATRIX_FILE_FORMAT_ERROR;
    }

    // Read colIndices
    uint64_t* colIndices;
    uint64_t col_indices_size = 0;
    read_result = _read_uint_64_array(file, &colIndices, &col_indices_size, NOT_LAST_LINE);
    switch (read_result) {
        case HEAP_MEMORY_ERROR:
        case MATRIX_FILE_FORMAT_ERROR:
            // error, col_indices already freed in function call
            free(values);
            fclose(file);
            return read_result;
        case MATRIX_READ_SUCCESS:
        default:
            break;
    }

    // Check if columns indices are valid
    if (_check_col_indices(colIndices, col_indices_size, values_size, noCols) == MATRIX_FILE_FORMAT_ERROR) {
        // col_indices freed in function call
        free(values);
        fclose(file);
        return MATRIX_FILE_FORMAT_ERROR;
    }

    // Read rowPointers and check if they are valid
    uint64_t row_pointers_size = 8;
    uint64_t* rowPointers;
    read_result = _read_uint_64_array(file, &rowPointers, &row_pointers_size, EOF);  // last line in the file!
    switch (read_result) {
        case HEAP_MEMORY_ERROR:
        case MATRIX_FILE_FORMAT_ERROR:
            // error, rowPointers array freed in function call
            free(values);
            free(colIndices);
            fclose(file);
            return read_result;
        case MATRIX_READ_SUCCESS:
        default:
            break;
    }

    // Check if the row pointers are valid
    if (_check_row_pointers(rowPointers, row_pointers_size, values_size, noRows, noCols) == MATRIX_FILE_FORMAT_ERROR) {
        // rowPointers freed in function call
        free(values);
        free(colIndices);
        fclose(file);
        return MATRIX_FILE_FORMAT_ERROR;
    }

    fclose(file);  // we are now done with the file

    // Initialize the CSR matrix
    *matrix = malloc(sizeof(Matrix));
    if (*matrix == NULL) {
        // free all arrays that were read
        free(values);
        free(colIndices);
        free(rowPointers);
        return HEAP_MEMORY_ERROR;
    }

    (*matrix)->noRows = noRows;
    (*matrix)->noCols = noCols;
    (*matrix)->values = values;
    (*matrix)->valuesSize = values_size;
    (*matrix)->colIndices = colIndices;
    (*matrix)->rowPointers = rowPointers;
    (*matrix)->rowPointersSize = row_pointers_size;

    return MATRIX_READ_SUCCESS;
}

int write_matrix_to_file(const char* filename, const Matrix* const matrix) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        return FILE_OPEN_ERROR;
    }

    // Write number of rows and ","
    int length = snprintf(NULL, 0, "%lu,", matrix->noRows);  // used to check if successfully written to file
    if (fprintf(file, "%lu,", matrix->noRows) != length) {
        fclose(file);
        return FILE_WRITE_ERROR;
    }

    // Write number of columns and newline character
    length = snprintf(NULL, 0, "%lu\n", matrix->noCols);
    if (fprintf(file, "%lu\n", matrix->noCols) != length) {
        fclose(file);
        return FILE_WRITE_ERROR;
    }

    // Write values
    if (_write_float_array(file, matrix->values, matrix->valuesSize) != ARRAY_WRITE_SUCCESS) {
        fclose(file);
        return FILE_WRITE_ERROR;
    }
    if (fwrite("\n", 1, 1, file) != 1) {  // write new line
        fclose(file);
        return FILE_WRITE_ERROR;
    }

    // Write col_indices
    if (_write_uint_64_array(file, matrix->colIndices, matrix->valuesSize) != ARRAY_WRITE_SUCCESS) {
        fclose(file);
        return FILE_WRITE_ERROR;
    }
    if (fwrite("\n", 1, 1, file) != 1) {  // write new line
        fclose(file);
        return FILE_WRITE_ERROR;
    }

    // Write row_ptr
    if (_write_uint_64_array(file, matrix->rowPointers, matrix->rowPointersSize) != ARRAY_WRITE_SUCCESS) {
        fclose(file);
        return FILE_WRITE_ERROR;
    }

    fclose(file);

    return MATRIX_WRITE_SUCCESS;
}

int set_error_message(char** error_message, const char* format, ...) {
    // Get variadic argument lists
    va_list ap1, ap2;
    va_start(ap1, format);
    va_copy(ap2, ap1);

    int length = vsnprintf(NULL, 0, format, ap1) + 1;  // +1 for the null byte
    *error_message = malloc(sizeof(char) * length);
    if (*error_message == NULL) {  // error allocating memory for the formatted string
        va_end(ap1);
        va_end(ap2);
        return HEAP_MEMORY_ERROR;
    }
    vsprintf(*error_message, format, ap2);

    // Clean-up
    va_end(ap1);
    va_end(ap2);

    return 0;
}

// HELPER FUNCTIONS BELOW //
// ---------------------- //

int _parse_filename(char** dest, const char* const src) {
    *dest = malloc_safe(sizeof(char), strlen(src) + 1);  // +1 for null byte
    if (*dest == NULL)
        return HEAP_MEMORY_ERROR;
    strcpy(*dest, src);
    return 0;
}

// READING
// -------

int _read_uint_64_array(FILE* file, uint64_t** array, uint64_t* array_size, const int eof_flag) {
    // Define variables used for reading
    int ch;
    char* endptr;  // used for strtoull
    int trailing_comma_flag = 0;

    size_t buffer_size = 128;  // 128 bytes should be more than enough for parsing
    char* buffer = malloc(sizeof(char) * buffer_size);
    if (buffer == NULL) {
        return HEAP_MEMORY_ERROR;
    }
    uint8_t buffer_index = 0;

    // Allocate memory to array
    uint64_t array_index = 0;
    *array_size = 8;
    *array = malloc_safe(sizeof(uint64_t), *array_size);
    if (*array == NULL) {
        free(buffer);
        return HEAP_MEMORY_ERROR;
    }

    uint64_t value;
    while (buffer_index < buffer_size) {
        ch = fgetc(file);

        if (ch >= '0' && ch <= '9') {
            // Read digit into buffer and reset trailing comma flag
            buffer[buffer_index++] = ch;
            trailing_comma_flag = 0;
        } else if (ch == eof_flag) {  // EOF or newline, stop parsing this line
            if (trailing_comma_flag || !array_index) {
                // we have an unmatched comma, e.g 2,5,
                // or a premature EOF/newline
                goto read_uint64_error;
            }

            // Add terminating null byte
            buffer[buffer_index] = '\0';

            // String to integer conversion
            errno = 0;
            value = strtoull(buffer, &endptr, 10);
            if (errno || *endptr != '\0') {
                goto read_uint64_error;
            }

            (*array)[array_index++] = value;

            buffer_index = buffer_size;  // invalidates loop condition, exits loop
            break;
        } else {
            switch (ch) {
                case ',':  // parse buffer and start reading another value
                    if (!buffer_index) {
                        // Comma with no value parsed, error
                        goto read_uint64_error;
                    }

                    // Add terminating null byte
                    buffer[buffer_index] = '\0';

                    // String to integer conversion
                    errno = 0;
                    value = strtoull(buffer, &endptr, 10);
                    if (errno || *endptr != '\0') {
                        goto read_uint64_error;
                    }

                    (*array)[array_index++] = value;
                    if (array_index == *array_size) {
                        // double array size
                        *array_size = (*array_size) * 2;
                        *array = realloc_safe(*array, sizeof(uint64_t), *array_size);
                        if (*array == NULL) {
                            free(buffer);
                            return HEAP_MEMORY_ERROR;
                        }
                    }

                    // Reset buffer index and set trailing comma flag
                    buffer_index = 0;
                    trailing_comma_flag = 1;
                    break;
                default:  // unmatched character
                    read_uint64_error:
                    free(*array);
                    free(buffer);
                    return MATRIX_FILE_FORMAT_ERROR;
            }
        }
    }

    free(buffer);

    // Resize array to just accommodate actual number of elements
    *array_size = array_index;
    *array = realloc_safe(*array, sizeof(uint64_t), *array_size);
    if (*array == NULL) {
        return HEAP_MEMORY_ERROR;
    }

    return MATRIX_READ_SUCCESS;
}

int _read_float_array(FILE* file, float** array, uint64_t* array_size, const int eof_flag) {
    // Define variables used for reading
    int ch;
    int decimal_point_flag = 0;
    char* endptr;  // used for strtof
    int trailing_comma_flag = 0;

    size_t buffer_size = 128;  // 128 bytes should be more than enough for parsing a float
    char* buffer = malloc(sizeof(char) * buffer_size);
    if (buffer == NULL) {
        fclose(file);
        return HEAP_MEMORY_ERROR;
    }
    uint8_t buffer_index = 0;

    // Allocate memory to array
    uint64_t array_index = 0;
    *array_size = 8;
    *array = malloc_safe(sizeof(float), *array_size);
    if (*array == NULL) {
        free(buffer);
        return HEAP_MEMORY_ERROR;
    }

    float value;
    while (buffer_index < buffer_size) {
        ch = fgetc(file);

        if (ch >= '0' && ch <= '9') {
            // Read digit into buffer and reset comma flag
            buffer[buffer_index++] = ch;
            trailing_comma_flag = 0;
        } else if (ch == eof_flag) {  // EOF or newline character
            if (trailing_comma_flag || !array_index) {
                // we have an unmatched comma, e.g 2,5,
                // or a premature EOF/newline
                goto read_float_error;
            }

            // Add terminating null byte at end of buffer
            buffer[buffer_index] = '\0';

            // String to float conversion
            errno = 0;
            value = strtof(buffer, &endptr);
            if (errno == ERANGE || *endptr != '\0') {
                goto read_float_error;
            }
            
            (*array)[array_index++] = value;

            buffer_index = buffer_size;  // exit loop
            break;
        } else {
            switch (ch) {
                case '.':
                    // Check if we already had a decimal point
                    if (decimal_point_flag) {
                        goto read_float_error;
                    } else {
                        buffer[buffer_index++] = '.';
                    }
                    break;
                case ',':  // parse buffer and start reading another value
                    if (!buffer_index) {
                        // Comma with no value parsed, error
                        goto read_float_error;
                    }
                    // Add terminating null byte at end of buffer
                    buffer[buffer_index] = '\0';

                    // String to float conversion
                    errno = 0;
                    value = strtof(buffer, &endptr);
                    if (errno == ERANGE || *endptr != '\0') {
                        goto read_float_error;
                    }

                    (*array)[array_index++] = value;
                    if (array_index == *array_size) {
                        // double array size
                        *array_size = (*array_size) * 2;
                        *array = realloc_safe(*array, sizeof(float), *array_size);
                        if (*array == NULL) {
                            return HEAP_MEMORY_ERROR;
                        }
                    }

                    // Set/reset flags and buffer index to read a new value
                    decimal_point_flag = 0;
                    trailing_comma_flag = 1;
                    buffer_index = 0;

                    break;
                case '-':
                    if (!buffer_index) {  // only the first character can be a minus
                        buffer[buffer_index++] = '-';
                    } else {
                        goto read_float_error;
                    }
                    break;
                default:  // unmatched character (EOF included), error
                    read_float_error:
                    free(*array);
                    free(buffer);
                    return MATRIX_FILE_FORMAT_ERROR;
            }
        }
    }

    // Resize array to accommodate just the actual number of elements in the array
    *array_size = array_index;
    *array = realloc_safe(*array, sizeof(float), *array_size);
    if (*array == NULL) {
        return HEAP_MEMORY_ERROR;
    }

    free(buffer);  // free buffer used for parsing

    return MATRIX_READ_SUCCESS;
}

int _check_values(
    float* values, const uint64_t values_size,
    const uint64_t noRows, const uint64_t noCols
) {
    if (values_size > noRows*noCols) {
        // e.g 3x5 matrix can have a maximum of 15 values but values_size is 20.
        values_check_error:
        free(values);
        return MATRIX_FILE_FORMAT_ERROR;
    }

    // Check for zeroes in array
    for (uint64_t i = 0; i < values_size; i++) {
        if (values[i] == 0.0f) {
            goto values_check_error;
        }
    }

    return 0;
}

int _check_col_indices (
    uint64_t* colIndices, const uint64_t col_indices_size, 
    const uint64_t values_size, const uint64_t noCols
    ) {
    // Check if the size of the column indices matches that of the values array
    if (col_indices_size != values_size) {
        col_indices_check_error: free(colIndices);
        return MATRIX_FILE_FORMAT_ERROR;
    }

    for (uint64_t i = 0; i < col_indices_size; i++) {
        if (colIndices[i] >= noCols) {
            goto col_indices_check_error;
        }
    }

    return 0;
}

int _check_row_pointers(
    uint64_t* rowPointers, const uint64_t row_pointers_size, 
    const uint64_t values_size, const uint64_t noRows, const uint64_t noCols
    ) {
    // At least two row pointers and exactly 1 more than number of rows
    if (row_pointers_size < 2 || row_pointers_size != noRows + 1) {
        row_pointers_check_error: free(rowPointers);
        return MATRIX_FILE_FORMAT_ERROR;
    }

    // The first row pointer must be 0
    uint64_t prev_row_ptr = rowPointers[0];
    if (prev_row_ptr) {
        goto row_pointers_check_error;
    }

    // The difference between two consecutive row pointers must be less than or equal to the number of columns
    uint64_t next_row_ptr;
    for (uint64_t i = 1; i < row_pointers_size; i++) {
        next_row_ptr = rowPointers[i];
        if (next_row_ptr - prev_row_ptr > noCols) {
            goto row_pointers_check_error;
        }
        prev_row_ptr = next_row_ptr;
    }

    // Last row pointer must be equal to number of values
    if (prev_row_ptr != values_size) {
        goto row_pointers_check_error;
    }

    return 0;
}

// WRITING
// -------

int _write_float_array(FILE* file, float* const array, size_t size) {
    if (size--) {  // check if size is greater than 0
        float value;
        int length;  // used to check if fprintf succeeded
        for (size_t i = 0; i < size; i++) {
            value = array[i];
            length = snprintf(NULL, 0, "%g,", value);
            if (fprintf(file, "%g,", value) != length) {
                return FILE_WRITE_ERROR;
            }
        }

        // Write last element without trailing comma
        value = array[size];
        length = snprintf(NULL, 0, "%g", value);
        if (fprintf(file, "%g", value) != length) {
            return FILE_WRITE_ERROR;
        }
    }

    return ARRAY_WRITE_SUCCESS;
}

int _write_uint_64_array(FILE* file, uint64_t* const array, size_t size) {
    if (size--) {  // check if size is greater than 0
        uint64_t value;
        int length;  // used to check if fprintf succeeded
        for (size_t i = 0; i < size; i++) {
            value = array[i];
            length = snprintf(NULL, 0, "%lu,", value);
            if (fprintf(file, "%lu,", value) != length) {
                return FILE_WRITE_ERROR;
            }
        }

        // Write last element without trailing comma
        value = array[size];
        length = snprintf(NULL, 0, "%lu", value);
        if (fprintf(file, "%lu", value) != length) {
            return FILE_WRITE_ERROR;
        }
    }

    return ARRAY_WRITE_SUCCESS;
}

// MEMORY MANAGEMENT FUNCTIONS BELOW //
// --------------------------------- //

void free_csr_matrix(Matrix* matrix) {
    if (matrix) {
        free(matrix->values);
        free(matrix->colIndices);
        free(matrix->rowPointers);
        free(matrix);
    }
}

void free_pointers(const int n, ...) {
    va_list args;
    va_start(args, n);

    void* ptr_to_free;
    for (int i = 0; i < n; i++) {
        ptr_to_free = va_arg(args, void*);
        free(ptr_to_free);
    }
    va_end(args);
}

void free_csr_matrices(const int n, ...) {
    va_list args;
    va_start(args, n);

    Matrix* matrix_to_free;
    for (int i = 0; i < n; i++) {
        matrix_to_free = va_arg(args, Matrix*);
        free_csr_matrix(matrix_to_free);
    }
    va_end(args);
}