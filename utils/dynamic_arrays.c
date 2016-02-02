#include "dynamic_arrays.h"

int*** allocateArrayI3(int dim1, int dim2, int dim3)
{ 
  int i,j;
  int ***array=(int ***)calloc(dim1,sizeof(int**));
  for (i=0;i<dim1;i++)
    {
      array[i] = (int **)calloc(dim2,sizeof(int*));
      for (j=0;j<dim2;j++)
        array[i][j] = (int *)calloc(dim3,sizeof(int));
    }

  return (int***)array;
}

int** allocateArrayI2(int dim1, int dim2)
{ 
  int i;
  int **array=(int **)calloc(dim1,sizeof(int**));
  for (i=0;i<dim1;i++)
      array[i] = (int *)calloc(dim2,sizeof(int*));

  return (int**)array;
}
