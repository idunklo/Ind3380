#include <mpi.h>   
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct
{
    int num_rows;
    int num_cols;
    double *matrix[0];
}
M; 

void allocate_matrix(M *mat, int num_rows, int num_cols)
{
    mat -> num_cols = num_cols;
    mat -> num_rows = num_rows;

    mat -> matrix[0] = (double*)calloc(num_rows*num_cols, sizeof(double));
}

void read_matrix_binaryformat (char* filename, double*** matrix,
        int* num_rows, int* num_cols)
{
    int i;
    FILE* fp = fopen (filename,"rb");
    fread (num_rows, sizeof(int), 1, fp);
    fread (num_cols, sizeof(int), 1, fp);
    /* storage allocation of the matrix */
    *matrix = (double**)malloc((*num_rows)*sizeof(double*));
    (*matrix)[0] = (double*)malloc((*num_rows)*(*num_cols)*sizeof(double));
    for (i=1; i<(*num_rows); i++)
        (*matrix)[i] = (*matrix)[i-1]+(*num_cols);
    /* read in the entire matrix */
    fread ((*matrix)[0], sizeof(double), (*num_rows)*(*num_cols), fp);
    fclose (fp);
}

void write_matrix_binaryformat (char* filename, double** matrix,
        int num_rows, int num_cols)
{
    FILE *fp = fopen (filename,"wb");
    fwrite (&num_rows, sizeof(int), 1, fp);
    fwrite (&num_cols, sizeof(int), 1, fp);
    fwrite (matrix[0], sizeof(double), num_rows*num_cols, fp);
    fclose (fp);
}    


int main(int argc, char *argv[])
{
    //Declaration of variables
    int my_rank, num_procs;  
    int num_rows, num_cols;
    double** matrix_a, matrix_b, matrix_c;   
    M A, B;
    


    /*
    //Declare MPI-suff
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    */
    return 0;  

}

