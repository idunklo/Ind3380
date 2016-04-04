//mpicc -cc=gcc-5 -fopenmp 
#include <stdio.h>
#include <stdlib.h>

// make use of two functions from the simplejpeg libraryi
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
    u -> image_data =  (float**)malloc(n*sizeof(float*));
    for (int i = 0; i < n; i++){
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
    for (int i = 0; i < u -> n; ++i){
        for (int j = 0; j < u -> m; ++j){
            u -> image_data[i][j] = (float)image_chars[i*u->n + j]; 
        }
    }
    return;
}

void convert_image_to_jpeg(const image *u, unsigned char* image_chars)  
{
    for (int i = 0; i < u -> n; ++i){
        for (int j = 0; j < u -> m; ++j){
            image_chars[(i*u->n) + j] =(unsigned char) u -> image_data[i][j];
        }
    }
    return;
}       

void iso_diffusion_denoising(image *u, image *u_bar, float kappa, int iters)
{ 
    int counter = 0;
    while(counter <= iters){
        for (int i = 1; i < u->n -1; i++){  //should it be -2 or -1?
            for (int j = 1;  j <  u -> m -1; j++){ //unsure about the same as above
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
    float kappa;
    image u, u_bar;
    unsigned char *image_chars;
    char *input_jpeg_filename, *output_jpeg_filename;

    /* read from command line: kappa, iters, input_jpeg_filename, output_jpeg_filename */
    kappa = atof(argv[1]);
    iters = atof(argv[2]);
    input_jpeg_filename = argv[3];
    output_jpeg_filename = argv[4];

    import_JPEG_file(input_jpeg_filename, &image_chars, &m, &n, &c);
    allocate_image (&u, m, n);
    allocate_image (&u_bar, m, n);  
    convert_jpeg_to_image (image_chars, &u);
    iso_diffusion_denoising (&u, &u_bar, kappa, iters);

    convert_image_to_jpeg (&u_bar, image_chars);
    export_JPEG_file(output_jpeg_filename, image_chars, m, n, c, 75);
    deallocate_image (&u);
    deallocate_image (&u_bar);  
    
    return 0;
}


