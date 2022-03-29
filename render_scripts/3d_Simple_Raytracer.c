#include "FPToolkit.c"
#include "M3d_matrix_tools.c"

#define M_YON 10000
#define M 100
#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 800

double obmat[M][4][4] ;
double obinv[M][4][4] ;
double color[M][3] ;
double reflectivity[M] ;
int (*grad[M])(double gradient[3], int onum, double intersection[3]);
void (*draw[M])(int onum);
double (*intersection[M])(double start[3], double change[3]);
int    num_objects ;



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// To support the light model :
double light_in_eye_space[3] ;
double AMBIENT      = 0.2 ;
double MAX_DIFFUSE  = 0.5 ;
double SPECPOW      = 50 ;



int Light_Model (double irgb[3],
                 double s[3],
                 double p[3],
                 double n[3],
                 double argb[3])
// s,p,n in eyespace

// irgb == inherent color of object (input to this function)
// s = location of start of ray (probably the eye)
// p = point on object (input to this function)
// n = normal to the object at p (input to this function)
// argb == actual color of object (output of this function)
// globals : AMBIENT, MAX_DIFFUSE, SPECPOW, light_in_eye_space[3]

// return 1 if successful, 0 if error
{

  double len ;
  double N[3] ; 
  len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]) ;
  if (len == 0) return 0 ;
  N[0] = n[0]/len ;  N[1] = n[1]/len ;  N[2] = n[2]/len ;

  double E[3] ;
  E[0] = s[0] - p[0] ; 
  E[1] = s[1] - p[1] ; 
  E[2] = s[2] - p[2] ; 
  len = sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]) ;
  if (len == 0) return 0 ;
  E[0] /= len ;  E[1] /= len ;  E[2] /= len ;
  double NdotE = N[0]*E[0] + N[1]*E[1] + N[2]*E[2] ;

  double L[3] ;
  L[0] = light_in_eye_space[0] - p[0] ; 
  L[1] = light_in_eye_space[1] - p[1] ; 
  L[2] = light_in_eye_space[2] - p[2] ; 
  len = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]) ;
  if (len == 0) return 0 ;
  L[0] /= len ;  L[1] /= len ;  L[2] /= len ;
  double NdotL = N[0]*L[0] + N[1]*L[1] + N[2]*L[2] ;





  double max_ambient_and_diffuse = AMBIENT + MAX_DIFFUSE ;
     // this needs to occur BEFORE you possibly jump to LLL below




  double intensity ;
  if (NdotL*NdotE < 0) {
    // eye and light are on opposite sides of polygon
    intensity = AMBIENT ; 
    goto LLL ;
  } else if ((NdotL < 0) && (NdotE < 0)) {
    // eye and light on same side but normal pointing "wrong" way
    N[0] *= (-1.0) ;    N[1] *= (-1.0) ;    N[2] *= (-1.0) ; 
    NdotE *= (-1.0) ;   // don't use NdotE below, probably should eliminate this
  }


  // ignore Blinn's variant
  double R[3] ; // Reflection vector of incoming light
  R[0] = 2*NdotL*N[0] - L[0] ;
  R[1] = 2*NdotL*N[1] - L[1] ;
  R[2] = 2*NdotL*N[2] - L[2] ;

  double EdotR = E[0]*R[0] + E[1]*R[1] + E[2]*R[2] ;

  double diffuse ;
  if (NdotL <= 0.0) { diffuse = 0.0 ; }
  else { diffuse = MAX_DIFFUSE*NdotL ; }

  double specular ;
  if (EdotR <= 0.0) { specular = 0.0 ; }
  else { specular = (1.0 - max_ambient_and_diffuse)*pow(EdotR,SPECPOW) ;}

  // printf("%lf %lf\n",diffuse,specular) ;
  intensity = AMBIENT + diffuse + specular ;



 LLL : ;

  double f,g ;
  if (intensity <= max_ambient_and_diffuse) {
    f = intensity / max_ambient_and_diffuse ;
    argb[0] = f * irgb[0] ;
    argb[1] = f * irgb[1] ;
    argb[2] = f * irgb[2] ;
  } else {
    f = (intensity - max_ambient_and_diffuse) / 
                           (1.0 - max_ambient_and_diffuse) ;
    g = 1.0 - f ;
    argb[0] = g * irgb[0] + f ;
    argb[1] = g * irgb[1] + f ;
    argb[2] = g * irgb[2] + f ;
  }

  return 1 ;
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int solve_quadratic(double a, double b, double c, double t[2])
{
  const double discriminant = b*b-4*a*c;
  if(discriminant< 0) {
	  return 0;
  } else if(discriminant== 0) {
	  t[0] = -b + sqrt(discriminant) / (2 * a);
	  return 1;
  }
  t[0] = (-b + sqrt(discriminant)) / (2 * a);
  t[1] = (-b - sqrt(discriminant)) / (2 * a);
  return 2;
  
}



// ================ Sphere stuff ===============
int sphere_grad(double gradient[3], int onum, double intersection[3]) {
  // for spheres, gradient is <2x, 2y, 2z>
  M3d_mat_mult_pt(gradient, obinv[onum], intersection);
  M3d_vector_mult_const(gradient, gradient, 2);
}

double sphere_intersection(double start[3], double change[3])
{
  // solve quadratic (only works for objects which started as unit sphere)
  double a = M3d_dot_product(change, change);
  double b = 2*M3d_dot_product(start, change);
  double c = M3d_dot_product(start, start) -1;

  double t[2];
  double t_closest = M_YON + 1;
  int num_solutions = solve_quadratic(a, b, c,  t);

  // find closest solution
  for(int i = 0; i < num_solutions; ++i) {
    if(t[i] < 0) {
      continue;
    }
    if(t[i] < t_closest) {
      t_closest = t[i];
    }
  }

  return t_closest;
}

// ==========================================================================


// ================ plane stuff ===============
int plane_grad(double gradient[3], int onum, double intersection[3]) {
  gradient[0] = 0;
  gradient[1] = 0;
  gradient[2] = 1;
}
double plane_intersection(double start[3], double change[3]) {
  if(change[2] == 0) {
    return M_YON + 1;
  }
  double t = -start[2] / change[2];
  if(t < 0) {
    return M_YON + 1;
  }

  double x = start[0] + t * change[0];
  double y = start[1] + t * change[1];
  if(-1 > x || x > 1 ||
     -1 > y || y > 1) {
    return M_YON + 1;
  }
  return t;
}




/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



int get_normal(double normal[3], int onum, double intersection[3]) {
  // get gradient
  double gradient[3];
  grad[onum](gradient, onum, intersection);

  // finally, transform back into object space
  double inv_T[4][4];
  M3d_transpose(inv_T, obinv[onum]);
  M3d_mat_mult_pt(normal, inv_T, gradient);

  M3d_normalize(normal, normal);
}

int ray_recursive(double Rsource[3], double Rtip[3], double argb[3], int n)
{

  double t_closest = M_YON + 1;
  int saved_onum = -1;


  // line = start + t*change
  for(int onum = 0; onum < num_objects; ++onum) {
    double tip[3];
    double start[3];
  
    // calculate and start
    M3d_vector_copy(start, Rsource);
    M3d_vector_copy(tip, Rtip);

    // transform line to object space
    M3d_mat_mult_pt(start, obinv[onum], start);
    M3d_mat_mult_pt(tip, obinv[onum], tip);

    // subtract start from change
    double change[3];
    M3d_vector_mult_const(change, start, -1);
    M3d_vector_add(change, change, tip);

    double t = intersection[onum](start, change);
    if(t < t_closest) {
      // new closest point! Save this one.
      t_closest = t;
      saved_onum = onum;
    }


  }

  if(t_closest > M_YON) {
    for(int j=0; j<3; ++j) {
      // didn't hit anything, so make color 0
      argb[j] = 0;
    }
	  return 0;
  }


  // get the point of intersection
  // first, calculate change
  double change[3];
  M3d_vector_mult_const(change, Rsource, -1);
  M3d_vector_add(change, change, Rtip);

  // then get point = Rsource + t(change)
  double point[3] = {0,0,0};
  // reduce t_closest slightly to avoid z-fighting
  t_closest *= 0.999;
  M3d_vector_mult_const(point, change, t_closest);
  M3d_vector_add(point, Rsource, point);

  // calculate normal
  double normal[3];
  get_normal(normal, saved_onum, point);

  // set color to this object's color:
  double o_color[3];
  Light_Model(color[saved_onum],
              Rtip, point, normal,
              o_color);


  // recurse!
  if (n > 0) {


    // calculate reflection angle
    // get transformation matrix to reflect across plane defined by normal vector
    double look[3], reflection[3];
    // look = intersection - source
    M3d_vector_mult_const(look, Rsource, -1);
    M3d_vector_add(look, look, point);
    M3d_normalize(look, look);

    // reflection = look - 2(look * normal)normal
    // get the point of intersection
    // first, calculate change
    M3d_vector_mult_const(reflection, normal, -2*M3d_dot_product(look, normal));
    M3d_vector_add(reflection, reflection, look);


    // new tip = reflection + intersection
    double new_tip[3];
    M3d_vector_add(new_tip, reflection, point);

    ray_recursive(point, new_tip, argb, n-1);

    // mix the colors
    // (for now, just have hardcoded 70% reflection)
    // TODO: add better support for reflection
  }

  for(int j=0; j<3; ++j) {
    const double r = reflectivity[saved_onum];
    argb[j] = (1-r) * o_color[j] + r * argb[j];
  }

}







//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

int ray(double Rsource[3], double Rtip[3], double argb[3])
{
  ray_recursive(Rsource, Rtip, argb, 4);
}


int screen_to_coords(double point[3], double screen_pos[2]) {
  point[0] = screen_pos[0] * 2 / SCREEN_WIDTH - 1;
  point[1] = screen_pos[1] * 2 / SCREEN_HEIGHT - 1;
  point[2] = 1;
}




int Draw_all(double light[3], double eye[3], double coi[3], double up[3])
{
  double view_mat[4][4], view_mat_inv[4][4];
  double argb[3];

  double Rsource[3];
  M3d_view(view_mat, view_mat_inv, eye, coi, up);

  // get objects into eye space
  for(int onum=0; onum < num_objects; ++onum) {
    M3d_mat_mult(obmat[onum], view_mat, obmat[onum]) ;
    M3d_mat_mult(obinv[onum], obinv[onum], view_mat_inv) ;
  }
  // also get light in eye space
  M3d_mat_mult_pt(light_in_eye_space, view_mat, light);

  // testing
  double testPoint[3] = {0,0,0};
  M3d_mat_mult_pt(testPoint, obmat[0], testPoint);
  printf("test point:\n");
  M3d_print_vector(testPoint);


  const double step = 1;

  for(int i=0; i<SCREEN_WIDTH; i+=step) {
    for(int j=0; j<SCREEN_HEIGHT; j+=step) {
      double screen_pos[2], Rtip[3];
      screen_pos[0] = i * 1.0;
      screen_pos[1] = j * 1.0;
      screen_to_coords(Rtip, screen_pos);


      ray (Rsource, Rtip, argb) ; 
      G_rgb(argb[0], argb[1], argb[2]);
      G_point(screen_pos[0], screen_pos[1]);

      if(argb[0] != 0) {
        /*
        printf("Rsource: \n");
        M3d_print_vector(Rsource);
        printf("Rtip: \n");
        M3d_print_vector(Rtip);
        */
      }
    }
  }

  // multiply object transformation matricies by inverse view matrix to reset
  for(int onum=0; onum < num_objects; ++onum) {
    M3d_mat_mult(obmat[onum], view_mat_inv, obmat[onum]) ;
    M3d_mat_mult(obinv[onum], obinv[onum], view_mat) ;
  }

}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


int test01()
{
  // ======================== Build Objects ==========================
  double Tvlist[100];
  int Tn, Ttypelist[100];

  num_objects = 0 ;

  //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////
  color[num_objects][0] = 1.0 ;
  color[num_objects][1] = 1.0 ; 
  color[num_objects][2] = 1.0 ;
  reflectivity[num_objects] = 0.9;

  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  2   ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  2   ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  2   ; Tn++ ;

  M3d_make_movement_sequence_matrix(obmat[num_objects], obinv[num_objects], Tn, Ttypelist, Tvlist);

  grad[num_objects] = plane_grad;
  intersection[num_objects] = plane_intersection;
  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////
  color[num_objects][0] = 1.0 ;
  color[num_objects][1] = 1.0 ; 
  color[num_objects][2] = 1.0 ;
  reflectivity[num_objects] = 0.9;

  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  2   ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  2   ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  -2  ; Tn++ ;
  Ttypelist[Tn] = RY ; Tvlist[Tn] =  60  ; Tn++ ;

  M3d_make_movement_sequence_matrix(obmat[num_objects], obinv[num_objects], Tn, Ttypelist, Tvlist);

  grad[num_objects] = plane_grad;
  intersection[num_objects] = plane_intersection;
  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////
  color[num_objects][0] = 1.0 ;
  color[num_objects][1] = 1.0 ; 
  color[num_objects][2] = 1.0 ;
  reflectivity[num_objects] = 0.9;

  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  2   ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  2   ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  -2  ; Tn++ ;
  Ttypelist[Tn] = RY ; Tvlist[Tn] = 120  ; Tn++ ;

  M3d_make_movement_sequence_matrix(obmat[num_objects], obinv[num_objects], Tn, Ttypelist, Tvlist);

  grad[num_objects] = plane_grad;
  intersection[num_objects] = plane_intersection;
  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////
  color[num_objects][0] = 1.0 ;
  color[num_objects][1] = 0.0 ; 
  color[num_objects][2] = 0.0 ;
  reflectivity[num_objects] = 0;

  Tn = 0 ;

  M3d_make_movement_sequence_matrix(obmat[num_objects], obinv[num_objects], Tn, Ttypelist, Tvlist);

  grad[num_objects] = sphere_grad;
  intersection[num_objects] = sphere_intersection;
  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////

  // =======================================================================

    
  // ====================== Set up eye & space =============================
    double eye[3], coi[3], up[3];
    double light_in_world_space[3];
    double t;

    coi[0] = 0;
    coi[1] = 0;
    coi[2] = 0;

    up[0] = 0;
    up[1] = 1;
    up[2] = 0;

    light_in_world_space[0] = 0;
    light_in_world_space[1] = 10;
    light_in_world_space[2] = 10;
  // =======================================================================

    t = 0;
    do
    {
      G_rgb(0,0,0) ;
      G_clear() ;

      // rotate view around
      eye[0] = 5*sin(t);
      eye[1] = 0; //4*sin(t/100);
      eye[2] = -5*cos(t);
      printf("eye: \n");
      M3d_print_vector(eye);

      up[0] = eye[0];
      up[1] = 1;
      up[2] = eye[2];
      t += 0.1;

      Draw_all(light_in_world_space, eye, coi, up);

    } while(G_wait_key() != 'q');

    G_save_image_to_file("3d_Simple_Raytracer.xwd") ;
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




int main()
{
  G_init_graphics(SCREEN_WIDTH,SCREEN_HEIGHT);
  test01() ;
}
