
#include "FPToolkit.c"
#include "M3d_matrix_tools.c"
#include "xwd_tools_03.c"

#define M_HITHER 0.001
#define M_YON 2000
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
int (*to_parametric[M])(double point[3], double P[2]);
double (*SDF[M])(double point[3]);

// debugging
int debug = 0;

int    num_objects ;
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

double phi(double point[3]) {
  //return 1;

  // for spheres, phi(x,y,z) = 2/(1 + x^2 + y^2 + z^2)
  return 2 / (1 + M3d_dot_product(point, point));
}

double phi_grad(double point[3], double grad[3]) {
  //M3d_vector_mult_const(point, point, 0);
  //return 0;

  // note: this is for 2d sphereical geometry only, z is ignored.

  // (taken from Paul Allen)
  // d(phi)/dx = -4x / (1+x^2+y^2)^2
  // d(phi)/dy = -4y / (1+x^2+y^2)^2
  double denom = (1+M3d_dot_product(point, point));
  denom = denom*denom;
  M3d_vector_mult_const(grad, point, -4/denom);
}

double V_second_deriv(double point[3], double V[3],   double second_deriv[3]) {
  //M3d_vector_mult_const(second_deriv, second_deriv, 0);
  //return 0;

  double grad[3];
  phi_grad(point, grad);

  // TODO: ask Paul why I have a -1 here, but he has a +1 here
  double denom = -phi(point);

  second_deriv[0] = (grad[0]*V[0]*V[0]
                     + 2*grad[1]*V[0]*V[1] - grad[0]*V[1]*V[1]
                     + 2*grad[2]*V[0]*V[2] - grad[0]*V[2]*V[2]
                    ) / denom;
  second_deriv[1] = (grad[1]*V[1]*V[1]
                     + 2*grad[0]*V[1]*V[0] - grad[1]*V[0]*V[0]
                     + 2*grad[2]*V[1]*V[2] - grad[1]*V[2]*V[2]
                    ) / denom;
  second_deriv[2] = (grad[2]*V[2]*V[2]
                     + 0*grad[0]*V[2]*V[0] - grad[2]*V[0]*V[0]
                     + 0*grad[1]*V[2]*V[1] - grad[2]*V[1]*V[1]
                    ) / denom;
}



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



const double scaling = 3;
const double D = 0.01;

const double tan_half = 1; // for now, hard-code the half angle to be 45 deg
const double H = tan_half;
const double HalfWinWidth  = 0.5*SCREEN_WIDTH;
const double HalfWinHeight = 0.5*SCREEN_HEIGHT;

// global view matrix
double view_mat[4][4], view_inv[4][4];

int screen_to_camera(double point[3], double screen_pos[2])
{
  point[0] = H*(screen_pos[0] - HalfWinWidth)/HalfWinWidth ;
  point[1] = H*(screen_pos[1] - HalfWinHeight)/HalfWinHeight ;
  point[2] = 0;
}

int screen_to_ray(double point[3], double screen_pos[2])
{
  screen_to_camera(point, screen_pos);
  point[2] = 1;
  M3d_vector_mult_const(point,    point, D);
  M3d_mat_mult_pt(point,    view_inv, point);
}

int screen_to_coords(double point[3], double screen_pos[2]) {
  screen_to_camera(point, screen_pos);
  M3d_mat_mult_pt(point,    view_inv, point);
}

int coords_to_screen(double screen_pos[2], double point[3]) {
  double tmp[3];
  M3d_mat_mult_pt(tmp,    view_mat, point);
  //printf("after view mat  : "); M3d_print_vector(tmp);
  screen_pos[0] = tmp[0] * HalfWinWidth   / H + HalfWinWidth;
  screen_pos[1] = tmp[1] * HalfWinHeight  / H + HalfWinHeight;

  /*
  screen_pos[0] = (tmp[0]/scaling + 1) * SCREEN_WIDTH / 2;
  screen_pos[1] = (tmp[1]/scaling + 1) * SCREEN_WIDTH / 2;
  */
}


int debug_draw_point(double point[3])
{
  if(!debug) {
    return 0;
  }
  double s[2];
  //printf("before view mat : "); M3d_print_vector(point);
  coords_to_screen(s, point);
  G_fill_circle(s[0], s[1], 0.5);
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


// gets an object's base color (prior to light model)
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

// casts a ray out from Rsource to Rtip, and checks if it intersects with any object.
// Returns:
//   - point[3] - point of intersection with object (typically slightly moved away from object)
//   - V[3]     - direction the ray came from
// (though because the geometry is non-Euclidean, will probably not reach Rtip)
// 
int cast_ray(double Rsource[3], double Rtip[3], double point[3], double V[3])
{
  const double delta_t = 0.0001;

  double dV_dt[3], next_point[3], obj_point[3];
  int onum;

  // V = Rtip - Rsource
  M3d_vector_mult_const(V, Rsource, -1);
  M3d_vector_add(V, Rtip, V);
  // currently, normalize V. TODO: should I keep doing this?
  M3d_normalize(V, V);


  M3d_vector_copy(next_point, Rsource);

  for(int n=0; n < 78000; ++n) {

    if(M3d_magnitude(Rsource) > M_YON) {
      return -1;
    }

    // get next point in the ray
    //
    // update point to be the old next_point, then
    // calculate new next_point and update V
    M3d_vector_copy(point, next_point);
    // next_point = point + (V * delta_t)
    M3d_vector_mult_const(next_point, V, delta_t);
    M3d_vector_add(next_point, next_point, point);
    // V += dV/dt
    V_second_deriv(next_point, V,   dV_dt);
    M3d_vector_mult_const(dV_dt,   dV_dt, delta_t);
    M3d_vector_add(V,   V, dV_dt);



    {
      // DEBUG
      //G_rgb(1,0,0);
      debug_draw_point(point);
      //G_wait_key();
    }

    // check for collision with any object
    for(onum=0; onum < num_objects; ++onum) {
      // transform point to object space
      M3d_mat_mult_pt(obj_point, obinv[onum], next_point);

      double sdf = SDF[onum](obj_point);
      if (sdf <= 0) {
        return onum;
      }
    }

  }

  //printf("no collision!\n");
  return -1;
}



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// To support the light model :
double light_in_eye_space[3] ;
double AMBIENT      = 1.0 ;
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


  /*
  // raycast to light, stop if you hit anything
  // first, move point out slightly from start, to prevent colliding with the same object
  double point[3];
  M3d_vector_mult_const(point, L, M_HITHER);
  M3d_vector_add(point, point, p);
  // then get the intersection
  double blocking_t;
  int blocking_onum = ray_get_closest_obj(point, light_in_eye_space, point, &blocking_t);

  if(blocking_onum != -1  && \
    blocking_t < 1)
  {
    intensity = AMBIENT;
    goto LLL ;
  }
  */

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
double sphere_SDF(double point[3]) {
  // f(x,y,z) = x^2 + y^2 + z^2 - 1
  return M3d_dot_product(point, point) - 1;
}

double inv_sphere_SDF(double point[3]) {
  // inverted sphere where in is out and out is in
  return -sphere_SDF(point);
}

int sphere_grad(double gradient[3], int onum, double intersection[3]) {
  // for spheres, gradient is <2x, 2y, 2z>
  M3d_mat_mult_pt(gradient, obinv[onum], intersection);
  M3d_vector_mult_const(gradient, gradient, 2);
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
double plane_SDF(double point[3]) {
  if(-1 > point[0] || point[0] > 1 ||
     -1 > point[1] || point[1] > 1) {
    return M_HITHER;
  }
  return point[2]*point[2];
}

int plane_grad(double gradient[3], int onum, double intersection[3]) {
  gradient[0] = 0;
  gradient[1] = 0;
  gradient[2] = 1;
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


int ray_to_rgb_recursive(double Rsource[3], double Rtip[3], double argb[3], int n)
{
  // default color to black
  for(int j=0; j<3; ++j) {
    argb[j] = 0;
  }

  // get point of intersection and look vector
  double point[3], look[3];
  int saved_onum = cast_ray(Rsource, Rtip, point, look);
  if(saved_onum == -1) {
    // no intersection
    return 0;
  }



  double normal[3];
  // for now, no light model
  // instead, just set color to this object's color:
  double o_color[3];
  get_color(saved_onum, point, o_color);



  /*
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
  */

  G_rgb(o_color[0], o_color[1], o_color[2]);
  debug_draw_point(point);
  //G_line(Rsource[0], Rsource[1], Rtip[0], Rtip[1]);

  double new_color[3] = {0, 0, 0};


  // recurse!
  if (n > 0 && reflectivity[saved_onum] > 0) {

    // calculate reflection angle
    // get transformation matrix to reflect across plane defined by normal vector
    /*
    // look = intersection - source
    M3d_vector_mult_const(look, Rsource, -1);
    M3d_vector_add(look, look, point);
    M3d_normalize(look, look);
    */

    // reflection = look - 2(look * normal)normal
    double reflection[3];
    M3d_vector_mult_const(reflection, normal, -2*M3d_dot_product(look, normal));
    M3d_vector_add(reflection, reflection, look);

    // move point out slightly, to avoid colliding with the same object again
    M3d_vector_mult_const(reflection, reflection, M_HITHER);
    M3d_vector_add(point, point, reflection);

    // new tip = reflection + intersection
    double new_tip[3];
    M3d_vector_add(new_tip, reflection, point);

    {
      // DEBUG
      G_rgb(o_color[0], o_color[1], o_color[2]);
      double a[2],b[2], debug_tip[3];
      M3d_vector_mult_const(debug_tip, reflection, 1);
      M3d_vector_add(debug_tip, debug_tip, point);

      coords_to_screen(a, point);
      coords_to_screen(b, debug_tip);

      G_line(a[0], a[1], b[0], b[1]);
      G_fill_circle(a[0], a[1], 4);
    }

    ray_to_rgb_recursive(point, new_tip, new_color, n-1);

  }

  for(int j=0; j<3; ++j) {
    const double r = reflectivity[saved_onum];
    argb[j] = (1-r) * o_color[j] + r * new_color[j];
  }

}







//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


int ray_to_rgb(double Rsource[3], double Rtip[3], double argb[3])
{
  ray_to_rgb_recursive(Rsource, Rtip, argb, 0);
}


/*
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


      ray_to_rgb (Rsource, Rtip, argb) ; 
      G_rgb(argb[0], argb[1], argb[2]);
      G_point(screen_pos[0], screen_pos[1]);

      if(argb[0] != 0) {
        // printf("Rsource: \n");
        // M3d_print_vector(Rsource);
        // printf("Rtip: \n");
        // M3d_print_vector(Rtip);
      }
    }
  }

  // multiply object transformation matricies by inverse view matrix to reset
  for(int onum=0; onum < num_objects; ++onum) {
    M3d_mat_mult(obmat[onum], view_mat_inv, obmat[onum]) ;
    M3d_mat_mult(obinv[onum], obinv[onum], view_mat) ;
  }

}
*/




//========================================================================================
//==================================== drawing shapes ====================================
//========================================================================================
void (*draw[M])(int onum);
void Draw_ellipsoid (int onum)
{
  int n,i ;
  double t, xyz[3] ;
  double x,y ;

  G_rgb (color[onum][0],color[onum][1],color[onum][2]) ;
  
  n = 10000 ;
  for (i = 0 ; i < n ; i++) {
    t = i*2*M_PI/n ;
    xyz[0] = cos(t) ;
    xyz[1] = sin(t) ;
    xyz[2] = 0 ;
    //printf("start           : "); M3d_print_vector(xyz);
    M3d_mat_mult_pt(xyz, obmat[onum], xyz) ;

    double p[2];
    //printf("before view mat : "); M3d_print_vector(xyz);
    //M3d_mat_mult_pt(xyz, view_mat, xyz) ;
    coords_to_screen(p, xyz);
    G_point(p[0], p[1]);
    //printf("final           : %lf, %lf\n", p[0], p[1]);
    /*
    x = xyz[0] ;
    y = xyz[1] ;
    G_point(x,y) ;
    */
  }
}


void Draw_plane (int onum)
{
  int n,i ;
  double t, xyz[3] ;
  double x,y ;

  G_rgb (color[onum][0],color[onum][1],color[onum][2]) ;
  
  n = 1000 ;
  for (i = 0 ; i < n ; i++) {
    t = (i * 1.0 /  n  - 0.5) * 2;
    xyz[0] = t ;
    xyz[1] = 0 ;
    xyz[2] = 0 ;
    M3d_mat_mult_pt(xyz, obmat[onum], xyz) ;
    x = xyz[0] ;
    y = xyz[1] ;
    G_point(x,y) ;
  }
}

void Draw_hyperbola(int onum)
{
  int n,i ;
  double t, xyz[3] ;
  double x,y ;

  G_rgb (color[onum][0],color[onum][1],color[onum][2]) ;
  
  n = 1000 ;
  for (i = 0 ; i < n ; i++) {
    t = i * 2 * M_PI / n - M_PI;
    xyz[0] = 1/cos(t) ;
    xyz[1] = tan(t)   ;
    xyz[2] = 0 ;
    if(-1 > xyz[1] || xyz[1] > 1) {
      continue;
    }
    M3d_mat_mult_pt(xyz, obmat[onum], xyz) ;
    x = xyz[0] ;
    y = xyz[1] ;
    G_point(x,y) ;
  }
}
// ==========================================================================





void Draw_the_scene()
{
  int onum ;
  for (onum = 0 ; onum < num_objects ; onum++) {
    draw[onum](onum) ;
  }
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////





int test01()
{
    double Tvlist[100];
    int Tn, Ttypelist[100];
    double m[4][4], mi[4][4];
    double Rsource[3];
    double Rtip[3];
    double argb[3] ;

    num_objects = 0 ;

    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 1.0 ;
    color[num_objects][1] = 1.0 ; 
    color[num_objects][2] = 1.0 ;
    color_type[num_objects] = SIMPLE_COLOR;
    reflectivity[num_objects] = 0.0;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  0.10   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  0.10   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_copy_mat(obmat[num_objects], m);
    M3d_copy_mat(obinv[num_objects], mi) ;

    grad[num_objects] = sphere_grad;
    draw[num_objects] = Draw_ellipsoid;
    SDF[num_objects] = sphere_SDF;
    num_objects++ ; // don't forget to do this
    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 0.3 ;
    color[num_objects][1] = 0.3 ; 
    color[num_objects][2] = 1.0 ;
    color_type[num_objects] = SIMPLE_COLOR;
    reflectivity[num_objects] = 0.0;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  2.10   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  2.10   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_copy_mat(obmat[num_objects], m) ;
    M3d_copy_mat(obinv[num_objects], mi) ;

    grad[num_objects] = sphere_grad;
    draw[num_objects] = Draw_ellipsoid;
    SDF[num_objects] = inv_sphere_SDF;
    num_objects++ ; // don't forget to do this        
    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 0.1 ;
    color[num_objects][1] = 0.5 ; 
    color[num_objects][2] = 0.0 ;
    color_type[num_objects] = SIMPLE_COLOR;
    reflectivity[num_objects] = 0.0;
	
    Tn = 0 ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  0.400   ; Tn++ ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  0.200   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] = -0.800   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  0.800   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_copy_mat(obmat[num_objects], m) ;
    M3d_copy_mat(obinv[num_objects], mi) ;

    grad[num_objects] = sphere_grad;
    draw[num_objects] = Draw_ellipsoid;
    SDF[num_objects] = sphere_SDF;
    num_objects++ ; // don't forget to do this        
    //////////////////////////////////////////////////////////////
    /*
    color[num_objects][0] = 0.5 ;
    color[num_objects][1] = 1.0 ; 
    color[num_objects][2] = 1.0 ;
    color_type[num_objects] = SIMPLE_COLOR;
    reflectivity[num_objects] = 0.0;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  0.130   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  0.30   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] = -0.15   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  0.100   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  2.600   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_copy_mat(obmat[num_objects], m) ;
    M3d_copy_mat(obinv[num_objects], mi) ;

    grad[num_objects] = sphere_grad;
    draw[num_objects] = Draw_ellipsoid;
    SDF[num_objects] = sphere_SDF;
    num_objects++ ; // don't forget to do this        
    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 0.5 ;
    color[num_objects][1] = 0.5 ; 
    color[num_objects][2] = 0.5 ;
    color_type[num_objects] = SIMPLE_COLOR;
    reflectivity[num_objects] = 0.0;
	
    Tn = 0 ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  0.800   ; Tn++ ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  0.800   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  0.900   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] = -1.000   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_copy_mat(obmat[num_objects], m) ;
    M3d_copy_mat(obinv[num_objects], mi) ;

    grad[num_objects] = sphere_grad;
    draw[num_objects] = Draw_ellipsoid;
    SDF[num_objects] = sphere_SDF;
    num_objects++ ; // don't forget to do this        
    //////////////////////////////////////////////////////////////
    */

    

    light_in_eye_space[0] = 500;
    light_in_eye_space[1] = 500;
    light_in_eye_space[2] = 0;

    // view matrix
    double eye[3], coi[3], up[3];
    double origin[3] = {0,0,0};

    eye[0] = 1;
    eye[1] = 0;
    eye[2] = 0;

    coi[0] = 0;
    coi[1] = 0;
    coi[2] = 0;

    up[0] = eye[0];
    up[1] = eye[1]+1;
    up[2] = eye[2];

    M3d_make_translation(view_mat,  eye[0],  eye[1],  eye[2]);
    M3d_make_translation(view_inv, -eye[0], -eye[1], -eye[2]);

    debug = 1;


    int draw_screen() {
      G_rgb(0,0,0) ;
      G_clear() ;

      double origin[3] = {0,0,0};
      double s[2], p[2], o[2];
      coords_to_screen(s, Rsource);
      coords_to_screen(p, Rtip);
      coords_to_screen(o, origin);
      G_rgb(1,0,1) ; G_fill_circle(s[0], s[1], 3) ;
      G_rgb(1,0,1) ; G_fill_circle(p[0], p[1], 3) ;
      G_rgb(1,0.5,0.2) ; G_fill_circle(o[0], o[1], 3) ;

      G_rgb(1,1,0) ; G_line(s[0], s[1], p[0], p[1]) ;
      //G_rgb(1,0,1) ; G_line(s[0]+50,s[1]-50,  s[0]+50,s[1]+50) ;



      Draw_the_scene() ;
    }

    int update_ray_points() {
    }



    draw_screen();
    
    // first frame, raytrace the whole line
    const int res = 10;
    double colors[SCREEN_HEIGHT][3];

    double p[2];
    M3d_mat_mult_pt(Rsource, view_inv, origin);
    for(int y_pix=0; y_pix<SCREEN_HEIGHT; y_pix += res) {
      int x_pix = SCREEN_WIDTH;
      p[0] = x_pix;
      p[1] = y_pix;

      screen_to_coords(Rtip, p);


      argb[0] = 1;
      argb[1] = 0;
      argb[2] = 1;
      ray_to_rgb (Rsource, Rtip, argb) ; 
      for(int j=0; j<3; ++j) {
        colors[y_pix][j] = argb[j];
      }
    }

    double p1[2], p2[2];
    G_rgb(1,0,1);
    G_fill_rectangle(HalfWinWidth-5, 0, 10, SCREEN_HEIGHT);


    for(int y_pix=0; y_pix<SCREEN_HEIGHT; y_pix += res) {
      G_rgb(colors[y_pix][0], colors[y_pix][1], colors[y_pix][2]);
      G_line(HalfWinWidth-2, y_pix, HalfWinWidth+2, y_pix);
    }
    G_save_image_to_file("2d_Simple_Raymarcher.xwd") ;

    while(1)
    {



      double p[2];
      G_wait_click(p);
      if(p[1] < 50) { break; }

      screen_to_ray(Rtip, p);
      Rtip[2] = 0;


      draw_screen();

      // draw ray
      ray_to_rgb (Rsource, Rtip, argb) ; 
      G_rgb(argb[0], argb[1], argb[2]);
      ray_to_rgb (Rsource, Rtip, argb) ; 
    }

    G_save_image_to_file("2d_Simple_Raymarcher.xwd") ;
}











int test02()
{
    double Tvlist[100];
    int Tn, Ttypelist[100];
    double m[4][4], mi[4][4];
    double Rsource[3];
    double Rtip[3];
    double argb[3] ;

    num_objects = 0 ;

    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 1.0 ;
    color[num_objects][1] = 1.0 ; 
    color[num_objects][2] = 1.0 ;
    color_type[num_objects] = SIMPLE_COLOR;
    reflectivity[num_objects] = 0.0;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  0.10   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  0.10   ; Tn++ ;
    Ttypelist[Tn] = SZ ; Tvlist[Tn] =  0.10   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_copy_mat(obmat[num_objects], m);
    M3d_copy_mat(obinv[num_objects], mi) ;

    SDF[num_objects] = sphere_SDF;
    num_objects++ ; // don't forget to do this
    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 0.0 ;
    color[num_objects][1] = 0.0 ; 
    color[num_objects][2] = 1.0 ;
    color_type[num_objects] = SIMPLE_COLOR;
    reflectivity[num_objects] = 0.0;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  0.20   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  0.20   ; Tn++ ;
    Ttypelist[Tn] = SZ ; Tvlist[Tn] =  0.10   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  0.50   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_copy_mat(obmat[num_objects], m);
    M3d_copy_mat(obinv[num_objects], mi) ;

    SDF[num_objects] = sphere_SDF;
    num_objects++ ; // don't forget to do this
    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 0.0 ;
    color[num_objects][1] = 0.0 ; 
    color[num_objects][2] = 1.0 ;
    color_type[num_objects] = SIMPLE_COLOR;
    reflectivity[num_objects] = 0.0;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  0.05   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  0.05   ; Tn++ ;
    Ttypelist[Tn] = SZ ; Tvlist[Tn] =  0.10   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  -0.30   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_copy_mat(obmat[num_objects], m);
    M3d_copy_mat(obinv[num_objects], mi) ;

    SDF[num_objects] = sphere_SDF;
    num_objects++ ; // don't forget to do this
    //////////////////////////////////////////////////////////////


    
    light_in_eye_space[0] = 0;
    light_in_eye_space[1] = 0;
    light_in_eye_space[2] = 10;

    // view matrix
    double eye[3], coi[3], up[3];
    double origin[3] = {0,0,0};

    eye[0] = 0;
    eye[1] = 0;
    eye[2] = 2;

    coi[0] = 0;
    coi[1] = 0;
    coi[2] = 0;

    up[0] = eye[0];
    up[1] = eye[1]+1;
    up[2] = eye[2];

    M3d_view(view_mat, view_inv, eye, coi, up);

    // first frame, raytrace the whole line
    const int res = 20;
    double colors[SCREEN_HEIGHT][3];

    G_rgb(0, 0.1, 0);
    G_clear();

    double p[2];
    M3d_mat_mult_pt(Rsource, view_inv, origin);
    for(int x_pix=0; x_pix<SCREEN_WIDTH; x_pix += res) {
      for(int y_pix=0; y_pix<SCREEN_HEIGHT; y_pix += res) {
    //for(int x_pix = 300; x_pix<500; x_pix += res) {
    //  for(int y_pix = 300; y_pix<500; y_pix += res) {
        p[0] = x_pix;
        p[1] = y_pix;

        screen_to_ray(Rtip, p);

        ray_to_rgb (Rsource, Rtip, argb) ; 

        G_rgb(argb[0], argb[1], argb[2]);
        G_point(x_pix, y_pix);
        G_display_image();
      }
    }

    printf("finished render\n");
    G_save_image_to_file("2d_Simple_Raymarcher.xwd") ;
    G_wait_key();
}




//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




int main()
{
  G_init_graphics(800,800);
  test02() ;
}
