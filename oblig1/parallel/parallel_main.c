#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// make use of two functions from the simplejpeg library

void print_flush(const char* str, int rank)
{
    printf("%d :: %s\n", rank, str);
    fflush(NULL);
}

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
}

void convert_image_to_jpeg(const image * u, unsigned char* image_chars)  
{
    int i, j;
    for (i = 0; i < u -> m; ++i){
        for (j = 0; j < u -> n; ++j){
            image_chars[i*u->n + j] =(unsigned char) u -> image_data[i][j];
        }
    }
}       

void iso_diffusion_denoising(image *u, image *u_bar, float kappa, int iters)
{ 
    int counter = 0;   
    int i, j;
    while(counter <= iters){
        for (i = 1; i < u->m -1; i++){
            for (j = 1;  j <  u -> n -1; j++){ //unsure about the same as above
                u_bar -> image_data[i][j] = u -> image_data[i][j]
                    + kappa*(u -> image_data[i-1][j] + u -> image_data[i][j-1] 
                            - 4*u -> image_data[i][j] + u -> image_data[i][j+1] + u -> image_data[i+1][j]);
            }
        }
        counter++;
    }

    /*
    tmp = u->image_data;
    u->image_data = u_bar->image_data;
    u_bar->image_data = tmp;  

    Does this work with work with mpi? 
    */ 
}     

int main(int argc, char *argv[])  
{
    int i, j;  //counter variables
    //Declaration of variables
    int m, n, c, iters;
    int my_m, my_rank, num_procs, tmp_m, rows_counter;
    float kappa;
    image u, u_bar, whole_image;
    unsigned char *image_chars, *my_image_chars;
    char *input_jpeg_filename, *output_jpeg_filename;
                 
    //check for enough arguments
    if(argc!=5){
        printf("read from command line: kappa, iters,  input_jpeg_filename, output_jpeg_filename\n");
        return 0;
    }
    
    //Desclare MPI-stuff
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &num_procs); 

    print_flush("MPI init", my_rank);

    MPI_Barrier(MPI_COMM_WORLD);

    /* read from command line: kappa, iters, input_jpeg_filename, output_jpeg_filename */
    kappa = atof(argv[2]);
    iters = atof(argv[1]);
    input_jpeg_filename = argv[3];
    output_jpeg_filename = argv[4];  
    
    //importing the whole image in process 0
    if (my_rank==0) {
        import_JPEG_file(input_jpeg_filename, &image_chars, &m, &n, &c);
        //allocate_image (&whole_image, m, n);
        print_flush("Read & allocated", my_rank);
    }

    //Broadcast so all m an n are the same in all processes
    MPI_Bcast (&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    print_flush("Broadcasted!!!", my_rank);

    /* divide the m x n pixels evenly among the MPI processes */
    my_m = m/num_procs;
    
    //Adding the rest to the last process
    if(my_rank < m%num_procs){
        my_m = my_m + 1;
    }

    //Add padding rows
    if(my_rank == 0 || my_rank == num_procs-1){
        my_m = my_m + 1;
    }
    else{
        my_m = my_m + 2;
    }

    
    // //Allocating image_parts
    allocate_image (&u, my_m, n);
    allocate_image (&u_bar, my_m, n);

    print_flush("Allocated image parts", my_rank);

    //process 0 sends small parts of the picture to the other processes
    if (my_rank == 0){
        rows_counter = my_m-1;

        printf("0 :: %d\n", my_m);

        for(i = 1; i < num_procs; i++){
            MPI_Recv(&tmp_m, 1, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            print_flush("Received size ...", my_rank);
            MPI_Send(&image_chars[(rows_counter-1)*n], tmp_m*n, MPI_UNSIGNED_CHAR, i, 0, MPI_COMM_WORLD);
            print_flush("Sent part ...", my_rank);
            rows_counter += tmp_m - 2; //minus 2 to account for padding
        }
        rows_counter += 1;  //to account for last proc only having on padding row
    }
    else{
        printf("%d :: %d\n", my_rank, my_m);
        MPI_Send(&my_m, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        print_flush("Sent parts size ...", my_rank);
        my_image_chars = (unsigned char*)malloc(my_m*n*sizeof(unsigned char));//allocates space for the image-parts
        print_flush("Allocated 1D uchar storage ...", my_rank);
        MPI_Recv(my_image_chars, my_m*n, MPI_UNSIGNED_CHAR, 0, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        print_flush("Received image data ...", my_rank);
    }

//    print_flush("Distributing data done!", my_rank);
//        
//    if(my_rank != 0) {
//        convert_jpeg_to_image (my_image_chars, &u);
//    }
//    else{
//        convert_jpeg_to_image (image_chars, &u);
//    }
//
//    print_flush("Converted 1D to image!", my_rank);
//
//    //iso_diffusion_denoising(&u, &u_bar, kappa, iters);
//
//    print_flush("ISO diff. done!", my_rank);
//
//    /* each process sends its resulting content of u_bar to process 0 */ 
//    if (my_rank != 0){
//        //remove padding rows
//        tmp_m = my_m - 2;
//        if(my_rank == num_procs-1){
//            tmp_m += 1;
//        }
//        //send u_bar
//        MPI_Send(&tmp_m, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
//        print_flush("Sent part size!", my_rank);
//        MPI_Send(u_bar.image_data[n], tmp_m*n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
//        print_flush("Sent part data!", my_rank);
//    }
//    else{
//        //copy rank 0 u_bar to whole_image
//        for (i = 0; i < my_m-1; ++i){
//            for (j = 0; j < n; ++j){
//                whole_image.image_data[i][j] = u_bar.image_data[i][j];
//            }
//        }
//        print_flush("Copied from u_bar!", my_rank);
//        // recieve data from other procs
//        rows_counter = my_m-1;
//        for(i = 1; i < num_procs; i++){
//            MPI_Recv(&tmp_m, 1, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//            print_flush("Received part size!", my_rank);
//            MPI_Recv(&whole_image.image_data[rows_counter][0], tmp_m*n, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //process 0 recives content og u_bar
//            print_flush("Received part data!", my_rank);
//            rows_counter += tmp_m;
//        }
//    }
//
//    print_flush("Done collecting parts!", my_rank);
//
//    /* process 0 receives from each process incoming values and */
//    /* copy them into the designated region of struct whole_image */
//    /* ... */
//    if (my_rank==0) {
//        convert_image_to_jpeg(&whole_image, image_chars);
//        export_JPEG_file(output_jpeg_filename, image_chars, m, n, c, 75);
//        deallocate_image (&whole_image);
//    }
//    deallocate_image(&u);
//    deallocate_image(&u_bar);
//    if(my_rank != 0) {
//        free(my_image_chars);
//    }

    MPI_Finalize ();
    return 0;

}
