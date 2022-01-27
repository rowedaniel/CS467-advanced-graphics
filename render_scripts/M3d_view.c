

#include "FPToolkit.c"
#include "M3d_matrix_tools.c"

int M3d_view(double v[4][4], double vi[4][4],  double eyeA[3], double coiA[3], double upA[3])
{

 int i ;
 int mtype[100] ;
 double mparam[100] ;

 // rescale stick-figure to match the screensize
 i = 0 ;
 mtype[i] = TX ;  mparam[i] =  -eyeA[0]                ; i++ ;
 mtype[i] = TY ;  mparam[i] =  -eyeA[1]                ; i++ ;
 mtype[i] = TY ;  mparam[i] =  -eyeA[2]                ; i++ ;
 mtype[i] = RZ ;  mparam[i] =  atan2(coiA[1], coiA[0]) ; i++ ;
 mtype[i] = RY ;  mparam[i] =  atan2(coiA[0], coiA[2]) ; i++ ;
 mtype[i] = RZ ;  mparam[i] =  atan2(coiA[0], coiA[1]) ; i++ ;

 M3d_make_movement_sequence_matrix(v,vi,  i,mtype,mparam) ;
}
