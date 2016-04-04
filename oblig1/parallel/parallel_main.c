#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// make use of two functions from the simplejpeg library

void import_JPEG_file(const char *filename, unsigned char **image_chars,
        int *image_height, int *image_width,
        int *num_components);     

void export_JPEG_file(const char *filename, unsigned char *image_chars,
        int image_height, int image_width,
        int num_components, int quality);  

//making data structure
typedef struct
{
    float* data_storage; /* 1D array to store data */
    float** image_data; /* a 2D array of floats */  
    int m;              /* # pixels in y-direction */
    int n;              /* # pixels in x-direction */
}
image; 

//allocate 2D array 
void allocate_image(image *u, int m, int n)
{                                            
    //setting # of pixels in u
    u -> m = m;
    u -> n = n;   

    //allocating matrix for image data 
    u -> data_storage = (float*)malloc(m*n*sizeof(float));
    u -> image_data =  (float**)malloc(m*sizeof(float*));

    int i;
    for (i = 0; i < m; i++){
        u -> image_data[i] = &u -> data_storage[i*n]; } ;
    
    return;
}

//free the storage
void deallocate_image(image *u)
{
    free(u -> data_storage);
    free(u -> image_data); 
}

void convert_jpeg_to_image(const unsigned char* image_chars, image *u)
{
    int i, j;
    for (i = 0; i < u -> m; ++i){
        for (j = 0; j < u -> n; j++){
            u -> image_data[i][j] = (float)image_chars[i*u->n + j]; 
        }
    }
    return;
}

void convert_image_to_jpeg(const image *u, unsigned char* image_chars)  
{
    int i, j;
    for (i = 0; i < u -> m; ++i){
        for (j = 0; j < u -> n; ++j){
            image_chars[(i*u->n) + j] =(unsigned char) u -> image_data[i][j];
        }
    }
    return;
}       

void iso_diffusion_denoising(image *u, image *u_bar, float kappa, int iters)
{ 
    int counter = 0;   
    int i, j;
    while(counter <= iters){
        for (i = 1; i < u->m -1; i++){  //should it be -2 or -1?
            for (j = 1;  j <  u -> n -1; j++){ //unsure about the same as above
                u_bar -> image_data[i][j] = u -> image_data[i][j]
                    + kappa*(u -> image_data[i-1][j] + u -> image_data[i][j-1] 
                    - 4*u -> image_data[i][j] + u -> image_data[i][j+1] + u -> image_data[i+1][j]);
            }
        }
    counter++;
    }
}     

int main(int argc, char *argv[])
{
    int m, n, c, iters;
    int my_m, my_n, my_rank, num_procs;
    float kappa;
    image u, u_bar, whole_image;
    unsigned char *image_chars, *my_image_chars;
    char *input_jpeg_filename, *output_jpeg_filename;
    MPI_Init (&argc, &argv);
    3
        MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &num_procs);

    /* read from command line: kappa, iters, input_jpeg_filename, output_jpeg_filename */
    /* ... */
    if (my_rank==0) {
        import_JPEG_file(input_jpeg_filename, &image_chars, &m, &n, &c);
        allocate_image (&whole_image, m, n);
    }
    MPI_Bcast (&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    /* divide the m x n pixels evenly among the MPI processes */
    my_m = ...;
    my_n = ...;
    allocate_image (&u, my_m, my_n);
    allocate_image (&u_bar, my_m, my_n);
    /* each process asks process 0 for a partitioned region */
    /* of image_chars and copy the values into u */
    /* ... */
    convert_jpeg_to_image (my_image_chars, &u);
    iso_diffusion_denoising_parallel (&u, &u_bar, kappa, iters);
    /* each process sends its resulting content of u_bar to process 0 */
    /* process 0 receives from each process incoming values and */
    /* copy them into the designated region of struct whole_image */
    /* ... */
    if (my_rank==0) {
        convert_image_to_jpeg(&whole_image, image_chars);
        export_JPEG_file(output_jpeg_filename, image_chars, m, n, c, 75);
        deallocate_image (&whole_image);
    }
    deallocate_image (&u);
    deallocate_image (&u_bar);
    MPI_Finalize ();
    return 0;

}
