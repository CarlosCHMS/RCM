#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<sys/time.h>
#include<omp.h>
#include <stdbool.h>
#include"utils.h"
#include"input.h"
#include"mesh.h"

int main(int argc, char **argv)
{

    struct timeval start, stop;
    
    gettimeofday(&start,NULL);
   
    MESH* mesh = meshInit("./mesh.su2", 1, 0);
          
    meshFree(mesh); 

    gettimeofday(&stop,NULL);
    
    printf("\nDuration %f s\n", duration(start, stop));

    return 0;

}
