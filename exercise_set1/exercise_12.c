#include <stdio.h>
#include <stdlib.h>

void fill1Darray(int n, int *v);

int main()
{  
    int n;
    n = 10;
    int *v = (int*)malloc(n*sizeof(int));

    fill1Darray(n, v);

    for( int i=0; i<n; i++){

         printf("%d ", v[i]);   
    }

    free(v);

    return 0;

}

void fill1Darray(int n, int *v)
{
    /* Filling an 1D array with random numbers */
    srand(8);
    for( int i=0; i<n; i++){
        v[i] = rand() % 50;
    }
}

int find_max()
{
    return 0;
}


