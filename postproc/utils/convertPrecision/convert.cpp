
#include <stdio.h> 
#include <string.h>

void double2float(FILE *fpIn, FILE *fpOut) {

  float xf;
  double xd; 

  int i = 0;

  while (fread(&xd,sizeof(xd),1,fpIn)) {
    xf = float(xd);
    fwrite(&xf,sizeof(xf),1,fpOut);
    // printf("%i,%f\n", i,x);
      i++;
  }

  printf("i = %d values converted\n", i);

}

void float2double(FILE *fpIn, FILE *fpOut) {

  float xf;
  double xd; 

  int i = 0;

  while (fread(&xf,sizeof(xf),1,fpIn)) {
    xd = double(xf);
    fwrite(&xd,sizeof(xd),1,fpOut);
    // printf("%i,%f\n", i,xf);
      i++;
  }

  printf("i = %d values converted\n", i);

}

int main(int argc, char *argv[]) {

  printf("%s\n%s\n%s\n", argv[1], argv[2], argv[3]);

  FILE *fpIn = fopen(argv[2], "rb");
  FILE *fpOut = fopen(argv[3], "wb");

  if (strcmp(argv[1],"f2d") == 0) {
    printf("convert float to double\n");
    float2double(fpIn, fpOut);
  }
  else if (strcmp(argv[1],"d2f") == 0) {
    printf("convert double to float\n");
    double2float(fpIn, fpOut);
  }
  else {
    printf("option not recognized\n");
  }

  fclose(fpIn);
  fclose(fpOut);

  return 0;
}





