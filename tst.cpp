#include <stdio.h>

int main(void)
{
#pragma omp parallel for
    for (int i=0; i<100; i++){
        printf("%d\n",i);
    }
}

