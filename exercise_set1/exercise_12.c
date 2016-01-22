#include <stdio.h>
#include <stdlib.h>

int main()
{  
    int n;
    n = 10;
    int *v = (int*)malloc(n*sizeof(int));

    for( int i=0; i<n; i++){
        v[i] = rand();

         printf("%d ", v[i]);   
    }

    free(v);

    return 0;
}

