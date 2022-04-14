#include "FPToolkit.c"
#include "M3d_matrix_tools.c"
#include "xwd_tools_03.c"

#define M_HITHER 0.01
#define M_YON 10000
#define M 100
#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 800
#define SIMPLE_COLOR 0
#define TEXTURE_COLOR 1

// color
double color_type[M] ;
double color[M][3] ;
int image_IDs[M];
int image_size[M][2];
double reflectivity[M] ;

// object info
double obmat[M][4][4] ;
double obinv[M][4][4] ;
int (*grad[M])(double gradient[3], int onum, double intersection[3]);
double (*intersection[M])(double start[3], double change[3]);
int (*to_parametric[M])(double point[3], double P[2]);

int    num_objects ;


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int get_color(int onum, double point[3], double rgb[3])
{
  if(color_type[onum] == SIMPLE_COLOR) {
    for(int i=0; i<3; ++i) {
      rgb[i] = color[onum][i];
    }
  } else if(color_type[onum] == TEXTURE_COLOR){
    // convert to object space
    double o_point[3];
    M3d_mat_mult_pt(o_point, obinv[onum], point);



    double P[2];
    to_parametric[onum](o_point, P);
    int x = floor( image_size[onum][0] * P[0] );
    int y = floor( image_size[onum][1] * P[1] );
    int e = get_xwd_map_color(image_IDs[onum], x,y,rgb);
  }
}



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// returns object number and t value for closest object intersected by ray
int ray_get_closest_obj(double Rsource[3], double Rtip[3], double point[3], double * ret_t)
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

  if(t_closest == M_YON + 1) {
    // no intersection, so give an onum of -1
    return -1;
  }
  (*ret_t) = t_closest;
  return saved_onum;
}

int cast_ray(double Rsource[3], double Rtip[3], double point[3])
{
  double t;
  int saved_onum = ray_get_closest_obj(Rsource, Rtip, point, &t);
  if(saved_onum == -1) {
    return -1;
  }


  // get the point of intersection
  // first, calculate change
  double change[3];
  M3d_vector_mult_const(change, Rsource, -1);
  M3d_vector_add(change, change, Rtip);

  // then get point = Rsource + t(change)
  M3d_vector_mult_const(point, change, t);
  M3d_vector_add(point, Rsource, point);
  return saved_onum;
}


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


  // flip normals if they're pointing the wrong direction
  if(NdotL < 0) {
    NdotL = -NdotL;
    NdotE = -NdotE;
    M3d_vector_mult_const(N, N, -1);
  }






  double max_ambient_and_diffuse = AMBIENT + MAX_DIFFUSE ;
     // this needs to occur BEFORE you possibly jump to LLL below




  double intensity ;

  // raycast to light, stop if you hit anything
  // first, move point out slightly from start, to prevent colliding with the same object
  double point[3];
  M3d_vector_mult_const(point, L, M_HITHER);
  M3d_vector_add(point, point, p);
  // then get the intersection
  double blocking_t;
  int blocking_onum = ray_get_closest_obj(point, light_in_eye_space, point, &blocking_t);
  // finally subtract p from point to get the vector

  if(blocking_onum != -1  && \
    blocking_t < 1)
  {
    intensity = AMBIENT;
    goto LLL ;
  }

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
int sphere_to_parametric(double point[3], double P[2])
{
  // parameterization = (cos(u)cos(v), sin(u)cos(v), sin(v))
  // so v=asin(z), u=atan2(y,x)
  // finally rescale so that 0<=u,v<=1
  P[1] = (point[1] + 1) / 2;
  if(point[1] <= -1 || point[1] >= 1) {
    P[0] = 0;
    return 0;
  }
  P[0] = (atan2(point[0]/sqrt(1-P[1]*P[1]), point[2]/sqrt(1-P[1]*P[1])) + M_PI) / M_PI / 2;
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
int plane_to_parametric(double point[3], double P[2])
{
  // plane parameterization is (u, v, 0)
  // rescale so that 0<=u,v<=1
  P[0] = (point[0] + 1) / 2;
  P[1] = (point[1] + 1) / 2;
}
// ==========================================================================



// ================ triangle stuff ===============
int triangle_grad(double gradient[3], int onum, double intersection[3]) {
  gradient[0] = 0;
  gradient[1] = 0;
  gradient[2] = 1;
}
double triangle_intersection(double start[3], double change[3]) {
  // triangle goes from points (0,0,0), (1,0,0), and (0,1,0)
  if(change[2] == 0) {
    return M_YON + 1;
  }
  double t = -start[2] / change[2];
  if(t < 0) {
    return M_YON + 1;
  }

  double x = start[0] + t * change[0];
  double y = start[1] + t * change[1];
  if(0 > x || x > 1 ||
     0 > y || y > x) {
    return M_YON + 1;
  }
  return t;
}
int triangle_to_parametric(double point[3], double P[2])
{
  // triangle parameterization is (u, v, 0)
  // rescale so that 0<=u,v<=1
  P[0] = point[0];
  P[1] = point[1];
}
// ==========================================================================

// ================ hyperboloid stuff ===============
int hyperboloid_grad(double gradient[3], int onum, double intersection[3]) {
  // for hyperboloids, gradient is <2x, -2y, 2z>
  M3d_mat_mult_pt(gradient, obinv[onum], intersection);
  M3d_vector_mult_const(gradient, gradient, 2);
  gradient[1] *= -1;
}

double hyperboloid_intersection(double start[3], double change[3])
{
  double change_negY[3] = {change[0], -change[1], change[2]};
  double start_negY[3] = {start[0], -start[1], start[2]};

  double a = M3d_dot_product(change, change_negY);
  double b = 2*M3d_dot_product(start, change_negY);
  double c = M3d_dot_product(start, start_negY) -1;

  double t[2];
  double t_closest = M_YON + 1;
  int num_solutions = solve_quadratic(a, b, c,  t);

  // find closest solution
  for(int i = 0; i < num_solutions; ++i) {
    if(t[i] < 0) {
      continue;
    }

    // make sure it's in the hyperbola bounds -1 <= y <= 1
    double y = start[1] + change[1] * t[i];
    if(-1 > y || y > 1) {
      continue;
    }

    if(t[i] < t_closest) {
      t_closest = t[i];
    }
  }

  return t_closest;
}
int hyperboloid_to_parametric(double point[3], double P[2])
{
  // hyperboloid parameterization is (cos(u)sec(v), sin(u)sec(v), tan(v))
  // inverse is u=atan2(z,x), v=atan(y)
  // rescale so that 0<=u,v<=1
  P[0] = (atan2(point[2], point[0]) + M_PI) / M_PI / 2;
  P[1] = (-atan(point[1]) + M_PI / 2) / M_PI;
}
// ==========================================================================



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
  // default color to black
  for(int j=0; j<3; ++j) {
    argb[j] = 0;
  }

  // get point of intersection
  double point[3];
  int saved_onum = cast_ray(Rsource, Rtip, point);
  if(saved_onum == -1) {
    // no intersection
    return 0;
  }


  // calculate normal
  double normal[3];
  get_normal(normal, saved_onum, point);

  // set color to this object's color:
  double o_color[3];
  get_color(saved_onum, point, o_color);
  Light_Model(o_color,
              Rtip,
              point,
              normal,
              o_color);

  double new_color[3] = {0, 0, 0};


  // recurse!
  if (n > 0 && reflectivity[saved_onum] > 0) {


    // calculate reflection angle
    // get transformation matrix to reflect across plane defined by normal vector
    double look[3], reflection[3];
    // look = intersection - source
    M3d_vector_mult_const(look, Rsource, -1);
    M3d_vector_add(look, look, point);
    M3d_normalize(look, look);

    // reflection = look - 2(look dot normal)normal
    M3d_vector_mult_const(reflection, normal, -2*M3d_dot_product(look, normal));
    M3d_vector_add(reflection, reflection, look);

    // move point out slightly, to avoid colliding with the same object again
    M3d_vector_mult_const(reflection, reflection, M_HITHER);
    M3d_vector_add(point, point, reflection);

    // new tip = reflection + intersection
    double new_tip[3];
    M3d_vector_add(new_tip, reflection, point);

    ray_recursive(point, new_tip, new_color, n-1);

  }

  for(int j=0; j<3; ++j) {
    const double r = reflectivity[saved_onum];
    argb[j] = (1-r) * o_color[j] + r * new_color[j];
  }

}







//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

int ray(double Rsource[3], double Rtip[3], double argb[3])
{
  ray_recursive(Rsource, Rtip, argb, 1);
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
  int error;

  num_objects = 0 ;

  //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////
  color_type[num_objects] = TEXTURE_COLOR;
  image_IDs[num_objects] = init_xwd_map_from_file("clock.xwd");
  if(image_IDs[num_objects] == -1) { printf("File load failure!\n"); exit(0); }
  error = get_xwd_map_dimensions(image_IDs[num_objects],image_size[num_objects]); 
  if(error == -1) { printf("File load failure!\n"); exit(0); }
  /*
  color_type[num_objects] = SIMPLE_COLOR;
  color[num_objects][0] = 1;
  color[num_objects][1] = 1;
  color[num_objects][2] = 1;
  */

  reflectivity[num_objects] = 0.0;

  Tn = 0 ;
  Ttypelist[Tn] = RX ; Tvlist[Tn] =  90   ; Tn++ ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  6   ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  6   ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  6   ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  -4   ; Tn++ ;

  M3d_make_movement_sequence_matrix(obmat[num_objects], obinv[num_objects], Tn, Ttypelist, Tvlist);

  grad[num_objects] = plane_grad;
  intersection[num_objects] = plane_intersection;
  to_parametric[num_objects] = plane_to_parametric;
  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////
  color_type[num_objects] = TEXTURE_COLOR;
  image_IDs[num_objects] = init_xwd_map_from_file("clock.xwd");
  if(image_IDs[num_objects] == -1) { printf("File load failure!\n"); exit(0); }
  error = get_xwd_map_dimensions(image_IDs[num_objects],image_size[num_objects]); 
  if(error == -1) { printf("File load failure!\n"); exit(0); }
  /*
  color_type[num_objects] = SIMPLE_COLOR;
  color[num_objects][0] = 1;
  color[num_objects][1] = 1;
  color[num_objects][2] = 1;
  */

  reflectivity[num_objects] = 1.0;

  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  6   ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  6   ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  6   ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  0   ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  -4   ; Tn++ ;

  M3d_make_movement_sequence_matrix(obmat[num_objects], obinv[num_objects], Tn, Ttypelist, Tvlist);

  grad[num_objects] = plane_grad;
  intersection[num_objects] = plane_intersection;
  to_parametric[num_objects] = plane_to_parametric;
  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////
  /*
  color_type[num_objects] = TEXTURE_COLOR;
  image_IDs[num_objects] = init_xwd_map_from_file("earth_jeff.xwd");
  if(image_IDs[num_objects] == -1) { printf("File load failure!\n"); exit(0); }
  error = get_xwd_map_dimensions(image_IDs[num_objects],image_size[num_objects]); 
  if(error == -1) { printf("File load failure!\n"); exit(0); }

  reflectivity[num_objects] = 0.0;

  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  2   ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  2   ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  2   ; Tn++ ;

  M3d_make_movement_sequence_matrix(obmat[num_objects], obinv[num_objects], Tn, Ttypelist, Tvlist);

  grad[num_objects] = sphere_grad;
  intersection[num_objects] = sphere_intersection;
  to_parametric[num_objects] = sphere_to_parametric;
  num_objects++ ; // don't forget to do this
  */
  //////////////////////////////////////////////////////////////
  for(int i=0; i<5; ++i) {
    color_type[num_objects] = TEXTURE_COLOR;
    to_parametric[num_objects] = plane_to_parametric;
    image_IDs[num_objects] = init_xwd_map_from_file("clock.xwd");
    if(image_IDs[num_objects] == -1) { printf("File load failure!\n"); exit(0); }
    error = get_xwd_map_dimensions(image_IDs[num_objects],image_size[num_objects]); 
    if(error == -1) { printf("File load failure!\n"); exit(0); }
    /*
    color_type[num_objects] = SIMPLE_COLOR;
    color[num_objects][0] = 1;
    color[num_objects][1] = 1;
    color[num_objects][2] = 1;
    */

    reflectivity[num_objects] = i/5.0;

    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  0.5   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  0.5   ; Tn++ ;
    Ttypelist[Tn] = SZ ; Tvlist[Tn] =  0.5   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  4     ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  -2    ; Tn++ ;
    Ttypelist[Tn] = RY ; Tvlist[Tn] =  i/5.0 * 360   ; Tn++ ;

    M3d_make_movement_sequence_matrix(obmat[num_objects], obinv[num_objects], Tn, Ttypelist, Tvlist);

    grad[num_objects] = hyperboloid_grad;
    intersection[num_objects] = hyperboloid_intersection;
  to_parametric[num_objects] = hyperboloid_to_parametric;
    num_objects++ ; // don't forget to do this

  }
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
    up[1] = 0;
    up[2] = 0;

    light_in_world_space[0] = 0;
    light_in_world_space[1] = 10;
    light_in_world_space[2] = 0;
  // =======================================================================

    t = 1;
    do
    {
      G_rgb(0,0,0) ;
      G_clear() ;

      // rotate view around
      eye[0] = 10*sin(t);
      eye[1] = 0; //4*sin(t/100);
      eye[2] = 10*cos(t);

      up[0] = eye[0];
      up[1] = eye[1]+1;
      up[2] = eye[2];
      t -= 0.1;

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
