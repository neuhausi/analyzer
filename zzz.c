
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#include <malloc.h>
#define EXTERN
#include "utils.h"
#include "files.h"
#include "methods.h"
#include "solution.h"
#include "function.h"
#include "matrix.h"
#include "metrics.h"
#include "var_metrics.h"
#include "stats.h"

int main(void)
{
 
    //static double aaa[20] = {2,5,1,7,19,8,16,20,17,10,14,9,3,15,13,4,6,11,12,18};
    //static double bbb[20] = {20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1};
    //static double ccc[20] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
    static double aaa[19] = {2,5,1,7,19,8,16,17,10,14,9,3,15,13,4,6,11,12,18};
    //static double bbb[19] = {19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1};
    //static double ccc[19] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
    
    int i, j, perm[30];

    sort_by_index(18, aaa, perm);

    for (i = 0; i < 19; i++) {
	j = perm[i];
	printf("i = %2i j = %2i data = %2.1f\n", i, j, aaa[j]);
    }

    return(0);

}
