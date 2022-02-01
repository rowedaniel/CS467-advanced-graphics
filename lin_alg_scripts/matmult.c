#include <stdio.h>


int M3d_print_mat (double a[3][3])
{
  int r,c ;
  for (r = 0 ; r < 3 ; r++ ) {
      for (c = 0 ; c < 3 ; c++ ) {
           printf(" %12.4lf ",a[r][c]) ;
      }
      printf("\n") ;
  }

  return 1 ;
} 



int M3d_copy_mat (double a[3][3], double b[3][3])
// a = b
{
  int r,c ;
  for (r = 0 ; r < 3 ; r++ ) {
      for (c = 0 ; c < 3 ; c++ ) {
           a[r][c] = b[r][c] ;
      }
  }

  return 1 ;
} 




int M3d_mat_mult (double res[3][3], double a[3][3], double b[3][3])
// res = a * b
// this is SAFE, i.e. the user can make a call such as 
// M3d_mat_mult(p,  p,q) or M3d_mat_mult(p,  q,p) or  M3d_mat_mult(p, p,p)
{
  double sum ;
  int k ;
  int r,c ;
  double tmp[3][3] ;

  for (r = 0 ; r < 3 ; r++ ) {
      for (c = 0 ; c < 3 ; c++ ) {
           sum = 0.0 ;
           for (k = 0 ; k < 3 ; k++) {
                 sum = sum + a[r][k]*b[k][c] ;
           }
           tmp[r][c] = sum ;
      }
  }


  M3d_copy_mat (res,tmp) ;

  return 1 ;
}
int main() {
	double E1[3][3] = {{0, 1, 0},
		     	   {1, 0, 0},
		     	   {0, 0, 1}};


	double E2[3][3] = {{1, 0, 0},
		     	   {0, 1, 0},
		     	   {-4,0, 1}};
	
	double E3[3][3] = {{1, 0, 0},
		     	   {0, 1, 0},
		     	   {0, 3, 1}};
	
	double E4[3][3] = {{1, 0, 0},
		     	   {0, 1,-1},
		     	   {0, 0, 1}};
	
	double E5[3][3] = {{1, 0, 0},
		     	   {0, 1, 0},
		     	   {0, 0, 1.0/2}};
	
	double E6[3][3] = {{1, 0,-3},
		     	   {0, 1, 0},
		     	   {0, 0, 1}};
	
	double A[3][3]  = {{0, 1, 2},
		     	   {1, 0, 3},
		     	   {4,-3, 8}};
	double C[3][3];
	double D[3][3];


	M3d_print_mat(A);
	printf("\n");
	M3d_mat_mult(C, E1, A);
	M3d_print_mat(C);
	printf("\n");
	M3d_mat_mult(C, E2, C);
	M3d_print_mat(C);
	printf("\n");
	M3d_mat_mult(C, E3, C);
	M3d_print_mat(C);
	printf("\n");
	M3d_mat_mult(C, E4, C);
	M3d_print_mat(C);
	printf("\n");
	M3d_mat_mult(C, E5, C);
	M3d_print_mat(C);
	printf("\n");
	M3d_mat_mult(C, E6, C);
	M3d_print_mat(C);
	printf("\n");
	printf("\n");
	printf("\n");



	M3d_mat_mult(C, E2, E1);
	M3d_mat_mult(C, E3, C);
	M3d_mat_mult(C, E4, C);
	M3d_mat_mult(C, E5, C);
	M3d_mat_mult(C, E6, C);
	printf("transformation matrix:\n");
	M3d_print_mat(C);
	printf("\n");
	printf("\n");
	printf("\n");


	printf("transformation matrix * A:\n");
	M3d_mat_mult(D, C, A);
	M3d_print_mat(D);
	printf("\n");
	
	printf("A * transformation matrix:\n");
	M3d_mat_mult(D, A, C);
	M3d_print_mat(D);
	printf("\n");


	
	/*
	*/

	
}


