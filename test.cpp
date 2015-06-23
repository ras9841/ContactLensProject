#include <stdio.h>

int main(){
    FILE *pFile;
    pFile = fopen("pressure.txt","r");  
    double a[10];
    float f;
    for(int i = 0; i<10; i++){
        fscanf(pFile, "%f", &f);
        a[i] = f;
        printf("F: %lf \n", a[i]);
    }
}
