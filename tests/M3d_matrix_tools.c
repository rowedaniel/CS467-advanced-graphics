#include <stdio.h>
#include <math.h>


/*

 ( x')          (x)
 ( y')  =   M * (y)  
 ( z')          (z)
 ( 1 )          (1)

instead of (x',y',z',1) = (x,y,z,1) * M  

*/




int M3d_print_mat (double a[4][4])
{
  int r,c ;
  for (r = 0 ; r < 4 ; r++ ) {
      for (c = 0 ; c < 4 ; c++ ) {
           printf(" %12.4lf ",a[r][c]) ;
      }
      printf("\n") ;
  }

  return 1 ;
} 





int M3d_copy_mat (double a[4][4], double b[4][4])
// a = b
{
  int r,c ;
  for (r = 0 ; r < 4 ; r++ ) {
      for (c = 0 ; c < 4 ; c++ ) {
           a[r][c] = b[r][c] ;
      }
  }

  return 1 ;
} 





int M3d_make_identity (double a[4][4])
// a = I
{
  int r,c ;
  for (r = 0 ; r < 4 ; r++ ) {
      for (c = 0 ; c < 4 ; c++ ) {
           if (r == c) a[r][c] = 1.0 ;
               else    a[r][c] = 0.0 ;
      }
  }

  return 1 ;
} 





int M3d_make_translation (double a[4][4], double dx, double dy, double dz)
{
  M3d_make_identity(a) ;
  a[0][3] =  dx ;  a[1][3] = dy ;  a[2][3] = dz ;
  return 1 ;
}





int M3d_make_scaling (double a[4][4], double sx, double sy, double sz)
{
  M3d_make_identity(a) ;
  a[0][0] =  sx ;  a[1][1] = sy ;  a[2][2] = sz ;
  return 1 ;
}












int M3d_make_x_rotation_cs (double a[4][4], double cs, double sn)
// this one assumes cosine and sine are already known
{
  M3d_make_identity(a) ;

  a[1][1] =   cs ;  a[1][2] = -sn ;
  a[2][1] =   sn ;  a[2][2] =  cs ;

  return 1 ;
}



int M3d_make_y_rotation_cs (double a[4][4], double cs, double sn)
// this one assumes cosine and sine are already known
{
  M3d_make_identity(a) ;

  a[0][0] =   cs ;  a[0][2] =  sn ;
  a[2][0] =  -sn ;  a[2][2] =  cs ;

  return 1 ;
}


int M3d_make_z_rotation_cs (double a[4][4], double cs, double sn)
// this one assumes cosine and sine are already known
{
  M3d_make_identity(a) ;

  a[0][0] =   cs ;  a[0][1] = -sn ;
  a[1][0] =   sn ;  a[1][1] =  cs ;

  return 1 ;
}





int M3d_mat_mult (double res[4][4], double a[4][4], double b[4][4])
// res = a * b
// this is SAFE, i.e. the user can make a call such as 
// M3d_mat_mult(p,  p,q) or M3d_mat_mult(p,  q,p) or  M3d_mat_mult(p, p,p)
{
  double sum ;
  int k ;
  int r,c ;
  double tmp[4][4] ;

  for (r = 0 ; r < 4 ; r++ ) {
      for (c = 0 ; c < 4 ; c++ ) {
           sum = 0.0 ;
           for (k = 0 ; k < 4 ; k++) {
                 sum = sum + a[r][k]*b[k][c] ;
           }
           tmp[r][c] = sum ;
      }
  }


  M3d_copy_mat (res,tmp) ;

  return 1 ;
}





int M3d_mat_mult_pt (double P[3],   double m[4][4], double Q[3])
// P = m*Q
// SAFE, user may make a call like M3d_mat_mult_pt (W, m,W) ;
{
  double u,v,t ;

  u = m[0][0]*Q[0] + m[0][1]*Q[1] + m[0][2]*Q[2] + m[0][3] ;
  v = m[1][0]*Q[0] + m[1][1]*Q[1] + m[1][2]*Q[2] + m[1][3] ;
  t = m[2][0]*Q[0] + m[2][1]*Q[1] + m[2][2]*Q[2] + m[2][3] ;  

  P[0] = u ;
  P[1] = v ;
  P[2] = t ;
  
  return 1 ;
}





int M3d_mat_mult_points (double X[], double Y[], double Z[],
                         double m[4][4],
                         double x[], double y[], double z[], int numpoints)
// |X0 X1 X2 ...|       |x0 x1 x2 ...|
// |Y0 Y1 Y2 ...| = m * |y0 y1 y2 ...|
// |Z0 Z1 Z2 ...|       |z0 z1 z2 ...|  
// | 1  1  1 ...|       | 1  1  1 ...|

// SAFE, user may make a call like M3d_mat_mult_points (x,y,z,  m, x,y,z,  n) ;
{
  double u,v,t ;
  int i ;

  for (i = 0 ; i < numpoints ; i++) {
    u = m[0][0]*x[i] + m[0][1]*y[i] + m[0][2]*z[i] + m[0][3] ;
    v = m[1][0]*x[i] + m[1][1]*y[i] + m[1][2]*z[i] + m[1][3] ;
    t = m[2][0]*x[i] + m[2][1]*y[i] + m[2][2]*z[i] + m[2][3] ;    

    X[i] = u ;
    Y[i] = v ;
    Z[i] = t ;
  }

  return 1 ;
}






int M3d_x_product (double res[3], double a[3], double b[3])
// res = a x b  , cross product of two vectors
// SAFE: it is ok to make a call such as
// D3d_x_product (a,  a,b) or
// D3d_x_product (b,  a,b) or
// D3d_x_product (a,  a,a) 
{
    double r[3] ;
    int v ;
    
    r[0] = a[1]*b[2] - b[1]*a[2] ;
    r[1] = b[0]*a[2] - a[0]*b[2] ;
    r[2] = a[0]*b[1] - b[0]*a[1] ;

    res[0] = r[0] ;
    res[1] = r[1] ;
    res[2] = r[2] ;

    if ((res[0] == 0) && (res[1] == 0) && (res[2] == 0)) {
	v = 0 ;
    } else {
	v = 1 ;
    }

    return v ;
}





//===========================================================================
// Custom Daniel stuff
//===========================================================================

int M3d_print_vector(double a[3])
{
	for(int i=0; i<3; i++) {
		printf("%lf ", a[i]);
	}
	printf("\n");
}

double M3d_dot_product(double a[3], double b[3]) {
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double M3d_magnitude(double a[3]) {
	return sqrt(M3d_dot_product(a, a));
}

int M3d_vector_mult_const(double res[3], double a[3], double b) {
	for(int i=0; i<3; ++i) {
		res[i] = a[i]+b;
	}
}

int M3d_vector_add(double res[3], double a[3], double b[3]) {
	for(int i=0; i<3; ++i) {
		res[i] = a[i]+b[i];
	}
}

int M3d_make_unit(double res[3], double a[3]) {
	double magnitude = M3d_magnitude(a);
	M3d_vector_mult_const(res, a, 1/magnitude);
}





//===========================================================================
// For Advanced Graphics :
//===========================================================================




#define SX 0
#define SY 1
#define SZ 2

#define RX 3
#define RY 4
#define RZ 5

#define TX 6
#define TY 7
#define TZ 8

#define NX 9
#define NY 10
#define NZ 11

int M3d_make_movement_sequence_matrix(double v[4][4], double vi[4][4],
                                      int n, int mtype[], double mparam[])
{
	double ftmp[4][4], btmp[4][4];

	M3d_make_identity(v);
	M3d_make_identity(vi);

	for(int i=0; i<n; ++i)
	{
		switch(mtype[i]) {
		case(SX):
			M3d_make_scaling(ftmp, mparam[i],     1.0, 1.0);
			M3d_make_scaling(btmp, 1.0/mparam[i], 1.0, 1.0);
			break;
		case(SY):
			M3d_make_scaling(ftmp, 1.0, mparam[i],     1.0);
			M3d_make_scaling(btmp, 1.0, 1.0/mparam[i], 1.0);
			break;
		case(SZ):
			M3d_make_scaling(ftmp, 1.0, 1.0,  mparam[i]);
			M3d_make_scaling(btmp, 1.0, 1.0,  1.0/mparam[i]);
			break;

		case(RX):
			{
			double angle = mparam[i] * M_PI/180.0;
			M3d_make_x_rotation_cs(ftmp, cos(angle), sin(angle));
			M3d_make_x_rotation_cs(btmp, cos(-angle), sin(-angle));
			break;
			}
		case(RY):
			{
			double angle = mparam[i] * M_PI/180.0;
			M3d_make_y_rotation_cs(ftmp, cos(angle), sin(angle));
			M3d_make_y_rotation_cs(btmp, cos(-angle), sin(-angle));
			break;
			}
		case(RZ):
			{
			double angle = mparam[i] * M_PI/180.0;
			M3d_make_z_rotation_cs(ftmp, cos(angle), sin(angle));
			M3d_make_z_rotation_cs(btmp, cos(-angle), sin(-angle));
			break;
			}

		case(TX):
			M3d_make_translation(ftmp, mparam[i],  0.0, 0.0);
			M3d_make_translation(btmp, -mparam[i], 0.0, 0.0);
			break;
		case(TY):
			M3d_make_translation(ftmp, 0.0, mparam[i],  0.0);
			M3d_make_translation(btmp, 0.0, -mparam[i], 0.0);
			break;
		case(TZ):
			M3d_make_translation(ftmp, 0.0, 0.0, mparam[i]);
			M3d_make_translation(btmp, 0.0, 0.0, -mparam[i]);
			break;

		case(NX):
			M3d_make_scaling(ftmp, -1.0, 1.0, 1.0);
			M3d_make_scaling(btmp, -1.0, 1.0, 1.0);
			break;
		case(NY):
			M3d_make_scaling(ftmp, 1.0, -1.0, 1.0);
			M3d_make_scaling(btmp, 1.0, -1.0, 1.0);
			break;
		case(NZ):
			M3d_make_scaling(ftmp, 1.0, 1.0,  -1.0);
			M3d_make_scaling(btmp, 1.0, 1.0,  -1.0);
			break;
		}
		M3d_mat_mult(v, ftmp, v);
		M3d_mat_mult(vi, vi, btmp);
	}
}
  


int M3d_view(double v[4][4], double vi[4][4],  double eyeA[3], double coiA[3], double upA[3])
{

 int i = 0;
 int mtype[100] ;
 double mparam[100] ;

 double new_coiA[3], new_upA[3];

 // translate eye to origin
 mtype[i] = TX ;  mparam[i] =  -eyeA[0]                         ; i++ ;
 mtype[i] = TY ;  mparam[i] =  -eyeA[1]                         ; i++ ;
 mtype[i] = TZ ;  mparam[i] =  -eyeA[2]                         ; i++ ;

 


 // rotate center of interest to go along z-axis
 
 // update new_coiA to reflect its new position
 M3d_make_movement_sequence_matrix(v,vi,  i,mtype,mparam) ;
 M3d_mat_mult_pt(new_coiA, v, coiA);

 // add new transformation (make coiA[1]=0)
 mtype[i] = RX ;  mparam[i] =  atan2(new_coiA[1], new_coiA[2])*180/M_PI  ; i++ ;

 // update new_coiA to reflect its new position
 M3d_make_movement_sequence_matrix(v,vi,  i,mtype,mparam) ;
 M3d_mat_mult_pt(new_coiA, v, coiA);
 // add new transformation (make coiA[0]=0)
 mtype[i] = RY ;  mparam[i] =  -atan2(new_coiA[0], new_coiA[2])*180/M_PI ; i++ ;






 // rotate up vector to point along y-axis
 
 // transform new_upA to reflect its new position
 M3d_make_movement_sequence_matrix(v,vi,  i,mtype,mparam) ;
 M3d_mat_mult_pt(new_upA, v, upA);

 // add new transformation to put up vector where it's supposed to be 
 mtype[i] = RZ ;  mparam[i] =  atan2(new_upA[0], new_upA[1])*180/M_PI    ; i++ ;
 M3d_make_movement_sequence_matrix(v,vi,  i,mtype,mparam) ;
 /*
 */

}
