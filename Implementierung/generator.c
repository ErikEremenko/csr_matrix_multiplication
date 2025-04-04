/*
This file generates input data for testing. As the tutor, you can execute the commands
below to generate test cases:

gcc -w -O3 -lm -mavx generator.c constants.c utils.c matrix.c matrixutils.c -o generate
./generate -s <seed>

You can use the -s flag to set a seed and generate deterministic test matrices.
Seed is set to time(NULL) if you don't specify the seed or specify an invalid seed.
*/

// We need this to silence VSCode errors for getopt
#define _POSIX_C_SOURCE 200809L

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <errno.h>

#include "csrmatrix.h"
#include "matrixutils.h"
#include "matrix.h"
#include "utils.h"

#define MAX_VALUE 100
#define FILENAME_SIZE 80
int generation_index = 0;

// This case should fail as there is a trailing new line at the end of the file.
const char* ERROR_CASE_0 = "3,4\n"
"1,5,2.,1,1\n"
"10,20\n"
"5,1\n\n\n";
// This case should fail as values contain a zero.
const char* ERROR_CASE_1 = "3,4\n"
"1,0,0,1,1\n"
"10,20\n"
"0,5";
// The test should fail as there is an unmatched trailing comma in a line (two lines).
const char* ERROR_CASE_2 = "10,10\n"
"5,1,2,4,\n"
"0,1,2,3,\n"
"0,4";
// The test should fail as there is a negative value for noRows
const char* ERROR_CASE_3 = "-10,1\n"
"1,1,1\n"
"0,0,0\n"
"0,2";
// The test should fail as there is an unmatched comma in a line.
const char* ERROR_CASE_4 = "10,10\n"
"5,1,2,,,4\n"
"0,0,0,0,1,\n"
"0,4";
// The test should fail as there is a noRows and noCols are 0
const char* ERROR_CASE_5 = "-10,1\n"
"1,1,1\n"
"0,0,0\n"
"0,1,2";
// The test should fail as there is a column index is more than the noCols
const char* ERROR_CASE_6 = "-10,1\n"
"1,1,1\n"
"0,500,0\n"
"0,1,2";
// This case should fail as there is a trailing new line in the values array.
const char* ERROR_CASE_7 = "3,4\n"
"1,3,5,1,1\n"
"0,1,2,1,2\n\n"
"0,2,5";
// This case should fail as there is a space character.
const char* ERROR_CASE_8 = "3,4\n"
"1,6,  7,1,8\n"
"0,1,0,1,2\n"
"0,2,5";
// This case should fail as there are no row pointers.
const char* ERROR_CASE_9 = "3,4\n"
"1,6,7,1,8\n"
"0,1,0,1,2";
// This case should fail as the row pointers don't make sense.
const char* ERROR_CASE_10 = "3,4\n"
"1,6,7,1,8\n"
"0,1,0,1,2\n"
"0,20,10,5,5";

int contains(uint64_t* arr, uint64_t size, uint64_t value) {
    for (uint64_t i = 0; i < size; i++) {
        if (arr[i] == value) { return 1; }
    }
    return 0;
}

void clear_array(uint64_t* arr, uint64_t size) {
    for (uint64_t i = 0; i < size; i++) {
        arr[i] = UINT64_MAX;
    }
}

void generate_error_case(char* filename, const char* error_str, int error_index) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
    error_case_1_error: fprintf(stderr, "Error generating error test case %d\n", error_index);
        return;
    }
    if (fprintf(file, error_str) != strlen(error_str)) {  // use strlen to see if fprintf wrote successfully
        goto error_case_1_error;
    }

    printf("Generation successful for error test case %d\n", error_index);

    fclose(file);
}

Matrix* gen_empty_matrix(uint64_t noRows, uint64_t noCols, uint64_t valuesSize) {
    if (valuesSize < 0)
    {
        fprintf(stderr, "maxNumValues must be greater than or equal 0\n");
    }
    if (valuesSize == 0) { valuesSize = 1; }

    Matrix* matrix = malloc(sizeof(Matrix));
    if (matrix == NULL) {
        return NULL;
    }

    float* values = calloc(sizeof(float), valuesSize);
    if (values == NULL) {
        free(matrix);
        return NULL;
    }

    uint64_t* colIndices = calloc(sizeof(uint64_t), valuesSize);
    if (colIndices == NULL) {
        free(matrix);
        free(values);
        return NULL;
    }

    uint64_t* rowPointers = calloc(sizeof(uint64_t), (noRows + 1));
    if (rowPointers == NULL) {
        free(matrix);
        free(values);
        free(colIndices);
        return NULL;
    }

    matrix->noRows = noRows;
    matrix->noCols = noCols;
    matrix->values = values;
    matrix->valuesSize = valuesSize;
    matrix->colIndices = colIndices;
    matrix->rowPointers = rowPointers;
    matrix->rowPointersSize = noRows + 1;

    return matrix;
}

// for example 10, 10, 15.0
Matrix* gen_rand_matrix(uint64_t noRows, uint64_t noCols, uint64_t valuesSize) {
    Matrix* matrix = gen_empty_matrix(noRows, noCols, valuesSize);
    if (matrix == NULL) {
        fprintf(stderr, "%s", HEAP_MEMORY_ERROR_MSG);
        return NULL;
    }

    if (valuesSize == 0)
    {
        return matrix;
    }


    for (uint64_t i = 0; i < valuesSize; i++) {
        matrix->values[i] = (float)((rand() % MAX_VALUE) + ((float)(rand() % 10) / 10.0));
        while (matrix->values[i] == 0.0)
        {
            matrix->values[i] = (float)(rand() % MAX_VALUE) + ((float)(rand() % 10) / 10.0);
        }


        if (rand() % 4 == 0) {
            matrix->values[i] *= -1;
        }

    }

    uint64_t num_elements_in_row = 0;
    uint64_t random_row = rand() % noRows;
    // Set row pointers
    for (uint64_t i = 0; i < valuesSize; i++) {
        // Returns an index from 1 to noCols-1
        do {
            random_row = 1 + (random_row % noRows);
            num_elements_in_row = matrix->rowPointers[random_row] - matrix->rowPointers[random_row - 1];
        } while (num_elements_in_row >= noCols);

        for (uint64_t j = random_row; j < noRows + 1; j++) {
            matrix->rowPointers[j]++;
        }
    }

    uint64_t row_index = 0;
    uint64_t used_col_counter = 0;
    uint64_t row_end = matrix->rowPointers[row_index + 1];
    uint64_t* used_col_array = malloc(valuesSize * sizeof(uint64_t));
    if (used_col_array == NULL)
    {
        fprintf(stderr, "%s", HEAP_MEMORY_ERROR_MSG);
        return NULL;
    }

    clear_array(used_col_array, valuesSize);
    // Set column indices
    uint64_t random_col;
    for (uint64_t i = 0; i < valuesSize; i++) {

        if (i >= row_end) {
            clear_array(used_col_array, valuesSize);
            used_col_counter = 0;
            row_index++;
            row_end = matrix->rowPointers[row_index + 1];

        }
        do {
            random_col = rand() % (noCols);
            // Set col_index
        } while (contains(used_col_array, used_col_counter, random_col) == 1);

        matrix->colIndices[i] = random_col;
        used_col_array[used_col_counter] = random_col;
        used_col_counter++;
    }
    free(used_col_array);

    return matrix;
}



void generate(uint64_t noRows, uint64_t noCols, uint64_t maxNumValues) {
    char filename[FILENAME_SIZE];
    snprintf(filename, FILENAME_SIZE, "generated/matrix_%d.txt", generation_index);
    Matrix* matrix = gen_rand_matrix(noRows, noCols, maxNumValues);
    if (gen_rand_matrix == NULL) {
        fprintf(stderr, "Generation failed for %d\n", generation_index);
    }
    if (write_matrix_to_file(filename, matrix) != MATRIX_WRITE_SUCCESS) {
        fprintf(stderr, "Error writing matrix to file for %d\n", generation_index);
    }
    free_csr_matrix(matrix);
    printf("Generation successful for %d\n", generation_index);
    generation_index++;
}


int main(int argc, char** argv) {
    // Check if seed was given
    int seed_flag = 0;
    unsigned int seed = 0;
    int ch;
    char* endptr;
    const char* optstring = "s:";
    while ((ch = getopt(argc, argv, optstring)) != -1) {
        switch (ch) {
        case 's':
            errno = 0;
            seed = (unsigned int)strtoul(optarg, &endptr, 10);
            if (errno || *endptr != '\0') {
                // error parsing string, no seed set
                fprintf(stderr, "Invalid seed, defaulting to using time...\n\n");
                break;
            }
            seed_flag = 1;
            break;
        default:
            break;
        }
    }

    // Initialize seed for random number generator here
    if (seed_flag) {
        srand(seed);
    }
    else {  // no seed given, use the current time as the seed
        srand(time(NULL));
    }

    // Standard cases
    generate(10, 10, 10);  // matrix 0, squared
    generate(300, 500, 30);  // matrix 1, can be multipied with 2
    generate(500, 300, 60);  // matrix 2, can be multipied with 1
    generate(5350, 1623, 1457);  // matrix 3, can be multiplied with 4
    generate(1623, 5350, 749);  // matrix 4, can be multiplied with 3
    generate(1, 100, 20);  // matrix 5, can be multiplied with 6
    generate(100, 1, 2);  // matrix 6, can be multiplied with 5


    // Edge cases for correctness tests
    generate(1, 1, 1);  // matrix 7
    generate(5, 5, 0);  // matrix 8

    // Matrix 9, this edge case shows the weakness of V1 (conversion to 2D array)
    generate(16, 16, 2);

    // B.I.G Matrices
    generate(1000, 1000, 10); // matrix 10 can be squared (extra sparse)
    generate(1000, 1000, 1000); // matrix 11 can be squared
    generate(1000, 1000, 50000); // matrix 12 can be squared
    generate(1000, 5000, 50000); // matrix 13, use 10-11-12
    generate(10000, 10000, 1000000);  // matrix 14

    // Error test cases
    int error_index = 0;

    generate_error_case("generated/error_matrix_0.txt", ERROR_CASE_0, error_index++);
    generate_error_case("generated/error_matrix_1.txt", ERROR_CASE_1, error_index++);
    generate_error_case("generated/error_matrix_2.txt", ERROR_CASE_2, error_index++);
    generate_error_case("generated/error_matrix_3.txt", ERROR_CASE_3, error_index++);
    generate_error_case("generated/error_matrix_4.txt", ERROR_CASE_4, error_index++);
    generate_error_case("generated/error_matrix_5.txt", ERROR_CASE_5, error_index++);
    generate_error_case("generated/error_matrix_6.txt", ERROR_CASE_6, error_index++);
    generate_error_case("generated/error_matrix_7.txt", ERROR_CASE_7, error_index++);
    generate_error_case("generated/error_matrix_8.txt", ERROR_CASE_8, error_index++);
    generate_error_case("generated/error_matrix_9.txt", ERROR_CASE_9, error_index++);
    generate_error_case("generated/error_matrix_10.txt", ERROR_CASE_10, error_index++);

    return 0;
}
