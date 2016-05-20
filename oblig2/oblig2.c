#include <mpi.h>   
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
   typedef struct
   {
   int num_rows;
   int num_cols;
   double *matrix[0];
   }
   M; 
   */

void print_flush(const char* str, int rank)
{
    printf("%d :: %s\n", rank, str);
    fflush(NULL);
}     

void allocate_matrix(double*** matrix, int num_rows, int num_cols)
{
    int i;
    *matrix = (double**)malloc((num_rows)*sizeof(double*));
    (*matrix)[0] = (double*)calloc((num_rows)*(num_cols),sizeof(double));  
    for (i=1; i<(num_rows); i++)
        (*matrix)[i] = (*matrix)[i-1]+(num_cols);  
    /*mat -> num_cols = num_cols;
      mat -> num_rows = num_rows;*/

}

void deallocate(double ***matrix){
    free((*matrix)[0]);
    free(*matrix);
}


void find_size( int num_rows, int num_cols, int *my_num_rows, int *my_num_cols, int my_rank, int sqrt_p)

{
    int x = my_rank % sqrt_p;;
    int y = (int)floor((double)my_rank / sqrt_p);

    *my_num_cols = floor((x+1) * (double)num_cols / sqrt_p) - floor(x*(double)num_cols/ sqrt_p);
    *my_num_rows = floor((y+1) * (double)num_rows / sqrt_p) - floor(y*(double)num_rows/ sqrt_p); 

    //using pointer to change the value outside the function. 
    //dividing number of columns and row. See example solution exam 2012.  
}    

void MatrixMultiply(int ra, int ca, int rb, int cb, double **a, double **b,
        double ***c)
{
    int i, j, k;

    for (i=0; i<ra; i++)
        for (j=0; j<cb; j++)
            for (k=0; k<rb; k++)
                (*c)[i][j] += a[i][k]*b[k][j];
}

void MatrixMatrixMultiply(double ***a, double ***b, double ***c, int mra, int
        mca, int mrb, int mcb, int *ra, int *ca, int *rb, int *cb, MPI_Comm
        comm)
{
    /*from the teaching book */
    int i, j;
    int num_procs, dims[2], periods[2];
    int myrank, my2drank, mycoords[2];
    int uprank, downrank, leftrank, rightrank, coords[2];
    int shiftsource, shiftdest;
    MPI_Status status; 
    MPI_Comm comm_2d;

    MPI_Comm_size(comm, &num_procs);
    MPI_Comm_rank(comm_2d, &my2drank);
    
    dims[0] = dims[1] = 0;
    MPI_Dims_create(num_procs, 2, dims);
    periods[0]= periods[1] = 1;

    MPI_Cart_create(comm, 2, dims, periods, 1, &comm_2d);
    MPI_Comm_rank(comm_2d, &my2drank);
    MPI_Cart_coords(comm_2d, my2drank, 2, mycoords);

    MPI_Cart_shift(comm_2d, 1, -1, &rightrank, &leftrank);
    MPI_Cart_shift(comm_2d, 0, -1, &downrank, &uprank);

    int ia = my2drank;
    int ib = my2drank;

    MPI_Cart_shift(comm_2d, 1, -mycoords[0], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace((*a)[0], mra*mca, MPI_DOUBLE, shiftdest, 1,
            shiftsource, 1, comm_2d, &status);
    MPI_Sendrecv_replace(&ia, 1, MPI_INT, shiftdest, 1, shiftsource, 1,
            comm_2d, &status);

    MPI_Cart_shift(comm_2d, 0, -mycoords[1], &shiftsource, &shiftdest);
    MPI_Sendrecv_replace((*b)[0], mrb*mcb, MPI_DOUBLE, shiftdest, 1,
            shiftsource, 1, comm_2d, &status);  
    MPI_Sendrecv_replace(&ib, 1, MPI_INT, shiftdest, 1, shiftsource, 1,
            comm_2d, &status);

    for (i=0; i<dims[0]; i++){
        MatrixMultiply(ra[ia], ca[ia], rb[ib], cb[ib], *a, *b, c); /* c=c + a*b */

        MPI_Sendrecv_replace((*a)[0], mra*mca, MPI_DOUBLE, leftrank, 1,
                rightrank, 1, comm_2d, &status);
        MPI_Sendrecv_replace((*b)[0], mrb*mcb, MPI_DOUBLE, uprank, 1, downrank,
                1, comm_2d, &status);

        MPI_Sendrecv_replace(&ia, 1, MPI_INT, leftrank, 1, rightrank, 1,
                comm_2d, &status);
        MPI_Sendrecv_replace(&ib, 1, MPI_INT, uprank, 1, downrank, 1,
                comm_2d, &status);
    }

    MPI_Comm_free(&comm_2d);    
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

void find_max(int *values, int *max, int val_size) {
    // find max value of a general array

    *max = values[0];
    for(int i=1; i<val_size; ++i) {
        if(values[i] > *max) {
            *max = values[i];
        }
    }
}

void fill_matrix(double ***matrix, int nr, int nc) {
    int k = 1;
    for(int i=0; i<nr; ++i) {
        for(int j=0; j<nc; ++j) {
            (*matrix)[i][j] = ((double)(i+j)*0.02) * (double)k;
        }
        k *=-1;
    }
}

void print_matrix(int nr, int nc, double **matrix, int rank, char id) {

    for(int i=0; i<nr; ++i) {
        for(int j=0; j<nc; ++j) {
            printf("%f ", matrix[i][j]);
        }
        printf("row: %i Matrix: %c, rank: %i\n", i, id, rank);
    }
    printf("\n");
}

void test_result(int nrA, int ncA, int nrB, int ncB, double **A, double **B,
        double **result) {

    double **serial_result;
    allocate_matrix(&serial_result, nrA, ncB);
    MatrixMultiply(nrA, ncA, nrB, ncB, A, B, &serial_result);

    int check = 0;
    for(int i=0; i<nrA; ++i) {
        for(int j=0; j<ncB; ++j) {
            if((int)fabs(result[i][j] - serial_result[i][j]) != 0) {
                printf("WRONG! (%i,%i) res: %f, real: %f\n", i,j, result[i][j], serial_result[i][j]);
                check = 1;
            }
        }
    }
    if(!check) {
        printf("Correct\n");
    }
}


int main(int argc, char *argv[])
{
    //Declaration of variables
    int my_rank, num_procs;  
    int num_rows, num_cols;
    double **A, **B, **C;
    int nrA, ncA, nrB, ncB;
    double **a, **b, **c;  
    int nra, nca, nrb, ncb;
    int sqrt_p;
    int *rsA, *csA, *rsB, *csB;
    int mrA, mcA, mrB, mcB;

    //Declare MPI-suff
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    char *Mat_A, *Mat_B, *Mat_C;

    Mat_A = argv[1];
    Mat_B = argv[2];
    Mat_C = argv[3]; 

    sqrt_p = sqrt(num_procs);

    if(my_rank == 0){
        //read_matrix_binaryformat(Mat_A, &A, &nrA, &ncA );
        //read_matrix_binaryformat(Mat_B, &B, &nrB, &ncB );
        nrA = 100;
        ncA = 50;
        nrB = 50;
        ncB = 100;
        allocate_matrix(&A, nrA, ncA);
        allocate_matrix(&B, nrB, ncB);
        fill_matrix(&A, nrA, ncA);
        fill_matrix(&B, nrB, ncB);
        char A_id = 'A', B_id = 'B';
        print_matrix(nrA, ncA, A, my_rank, A_id);
        print_matrix(nrB, ncB, B, my_rank, B_id);
    }

    MPI_Bcast(&nrA, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ncA, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nrB, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ncB, 1, MPI_INT, 0, MPI_COMM_WORLD);  

    find_size(nrA, ncA, &nra, &nca, my_rank, sqrt_p); 
    find_size(nrB, ncB, &nrb, &ncb, my_rank, sqrt_p);

    rsA = (int*)malloc(num_procs*sizeof(int)); 
    csA = (int*)malloc(num_procs*sizeof(int)); 
    rsB = (int*)malloc(num_procs*sizeof(int)); 
    csB = (int*)malloc(num_procs*sizeof(int)); 

    for(int i=0; i<num_procs; ++i) {
        if(i == my_rank) {
            rsA[i] = nra;
            csA[i] = nca;
            rsB[i] = nrb;
            csB[i] = ncb;
        } else {
            MPI_Sendrecv(&nra, 1, MPI_INT, i, 1, &(rsA[i]), 1, MPI_INT, i, 1,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&nca, 1, MPI_INT, i, 1, &(csA[i]), 1, MPI_INT, i, 1,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&nrb, 1, MPI_INT, i, 1, &(rsB[i]), 1, MPI_INT, i, 1,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&ncb, 1, MPI_INT, i, 1, &(csB[i]), 1, MPI_INT, i, 1,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        //printf("nra[%i] = %i, my_rank:%i\n", i, rsA[i], my_rank);
    }

    // root finds large enough sizes
    if(my_rank == 0) {
        find_max(rsA, &mrA, num_procs);
        find_max(csA, &mcA, num_procs);
        find_max(rsB, &mrB, num_procs);
        find_max(csB, &mcB, num_procs);
    }

    // send bigbuffs to all
    MPI_Bcast(&mrA, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mcA, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mrB, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mcB, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // allocate sub matrices
    allocate_matrix(&a, mrA, mcA);
    allocate_matrix(&b, mrB, mcB);
    allocate_matrix(&c, mrA, mcB);

    if(my_rank == 0) {

        // root sets its own matrices
        for(int i=0; i<nra; ++i) {
            for(int j=0; j<nca; ++j) {
                a[i][j] = A[i][j];
            }
        }
        for(int i=0; i<nrb; ++i) {
            for(int j=0; j<ncb; ++j) {
                b[i][j] = B[i][j];
            }
        }

        MPI_Datatype btA;
        MPI_Datatype btB;

        // initialize row and column index variables
        int riA = 0;
        int riB = 0;
        int ciA = csA[0];
        int ciB = csB[0];
        int tmpk = 1;
        for(int k=1; k<num_procs; ++k) {
            //send to slaves 
            // create strided types (blocks/chunks/blarg)
            MPI_Type_vector(rsA[k], csA[k], ncA, MPI_DOUBLE, &btA);
            MPI_Type_create_resized(btA, 0, sizeof(double), &btA);
            MPI_Type_commit(&btA);
            MPI_Type_vector(rsB[k], csB[k], ncB, MPI_DOUBLE, &btB);
            MPI_Type_create_resized(btB, 0, sizeof(double), &btB);
            MPI_Type_commit(&btB);

            // send to slaves
            MPI_Send(&((*A)[ncA*riA + ciA]), 1, btA, k, 1, MPI_COMM_WORLD);
            MPI_Send(&((*B)[ncB*riB + ciB]), 1, btB, k, 2, MPI_COMM_WORLD);

            // free types for next iteration
            MPI_Type_free(&btA);
            MPI_Type_free(&btB);

            ciA += csA[k];
            ciB += csB[k];
            tmpk++;

            if(tmpk == sqrt_p) {
                // jump down to next chunk row
                riA += rsA[k-1];
                riB += rsB[k-1];
                ciA = 0;
                ciB = 0;
                tmpk = 0;
            }
        }

        // deallocate A and B
        //deallocate(&A);
        //deallocate(&B);
    } else {
        // slaves recv from root a chunk/block/somethgin

        // create data types
        MPI_Datatype btA;
        MPI_Datatype btB;
        
        // create strided data types
        MPI_Type_vector(nra, nca, mcA, MPI_DOUBLE, &btA);
        MPI_Type_create_resized(btA, 0, sizeof(double), &btA);
        MPI_Type_commit(&btA);
        MPI_Type_vector(nrb, ncb, mcB, MPI_DOUBLE, &btB);
        MPI_Type_create_resized(btB, 0, sizeof(double), &btB);
        MPI_Type_commit(&btB);

        // recv from root a chunk
        MPI_Recv(&((*a)[0]), 1, btA, 0, 1, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
        MPI_Recv(&((*b)[0]), 1, btB, 0, 2, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);

        // free types
        MPI_Type_free(&btA);
        MPI_Type_free(&btB);
        //char a_id = 'a', b_id = 'b';
        //print_matrix(nra, nca, a, my_rank, a_id);
    }

    char a_id = 'a', b_id = 'b';
    //print_matrix(nra, nca, a, my_rank, a_id);
    //print_matrix(nrb, ncb, b, my_rank, b_id);

    MatrixMatrixMultiply(&a, &b, &c, mrA, mcA, mrB, mcB, rsA, csA, rsB, csB,
            MPI_COMM_WORLD);

    if(my_rank == 0) {

        // allocate result
        allocate_matrix(&C, nrA, ncB);

        // root sets its own matrices
        for(int i=0; i<nra; ++i) {
            for(int j=0; j<ncb; ++j) {
                C[i][j] = c[i][j];
            }
        }

        MPI_Datatype btC;

        int riC = 0;
        int ciC = csB[0];
        int tmpk = 1;
        for(int k=1; k<num_procs; ++k) {
            MPI_Type_vector(rsA[k], csB[k], ncB, MPI_DOUBLE, &btC);
            MPI_Type_create_resized(btC, 0, sizeof(double), &btC);
            MPI_Type_commit(&btC);

            printf("%i\n", ncB*riC+ciC);
            MPI_Recv(&((*C)[ncB*riC+ciC]), 1, btC, k, k, MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE);

            MPI_Type_free(&btC);

            ciC += csB[k];

            if(tmpk == sqrt_p) {
                riC += rsA[k-1];
                ciC = 0;
                tmpk = 0;
            }
        }

        //read_matrix_binaryformat(Mat_A, &A, &nrA, &ncA );
        //read_matrix_binaryformat(Mat_B, &B, &nrB, &ncB );
        test_result(nrA, ncA, nrB, ncB, A, B, C);
        //test_result(nrA, ncA, nrB, ncB, A, B, C);
        deallocate(&A);
        deallocate(&B);
        deallocate(&C);
    } else {
        MPI_Datatype btC;
        MPI_Type_vector(nra, ncb, mcB, MPI_DOUBLE, &btC);
        MPI_Type_create_resized(btC, 0, sizeof(double), &btC);
        MPI_Type_commit(&btC);

        MPI_Send(&((*c)[0]), 1, btC, 0, my_rank, MPI_COMM_WORLD);

        MPI_Type_free(&btC);
    }

    // deallocate sub matrices
    deallocate(&a); 
    deallocate(&b);
    deallocate(&c);

    free(rsA);
    free(csA);
    free(rsB);
    free(csB);
    
    MPI_Finalize();
    return 0;  

}

