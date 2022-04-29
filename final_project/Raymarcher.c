
//#include "FPToolkit.c"
#include "M3d_matrix_tools.c"
#include "xwd_tools_03.c"

#define M_HITHER 0.001
#define M_YON 2000
#define M 100
#define SCREEN_WIDTH 800
#define SCREEN_HEIGHT 800
#define SIMPLE_COLOR 0
#define TEXTURE_COLOR 1
#define MAX_LIGHT_OBJ_RES  200

// color
double color_type[M] ;
double color[M][3] ;
int image_IDs[M];
int image_size[M][2];
double reflectivity[M] ;

// object info
double obmat[M][4][4] ; // object space -> world space
double obinv[M][4][4] ; // world space -> object space
int (*grad[M])(double gradient[3], int onum, double intersection[3]); // used in normal calc
int (*to_parametric[M])(double point[3], double P[2]); // used in texture mapping
double (*SDF[M])(double point[3]); // used in raymarch collision detection
double baked_lights[M][MAX_LIGHT_OBJ_RES][MAX_LIGHT_OBJ_RES][3];

double (*phi)(double point[3]);
double (*phi_grad)(double point[3], double grad[3]);

// debugging
int debug = 0;
int do_lightmodel = 0;
int id;

int    num_objects ;
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

double euclid_phi(double point[3]) {
  return 1;
}

double euclid_phi_grad(double point[3], double grad[3]) {
  for(int i=0; i<3; ++i) {
    grad[i] = 0;
  }
  return 0;
}



double sphere_phi(double point[3]) {
  // for spheres, phi(x,y,z) = 2/(1 + x^2 + y^2 + z^2)
  return 2 / (1 + M3d_dot_product(point, point));
}

double sphere_phi_grad(double point[3], double grad[3]) {
  // note: this is for 2d sphereical geometry only, z is ignored.

  // (taken from Paul Allen)
  // d(phi)/dx = -4x / (1+x^2+y^2)^2
  // d(phi)/dy = -4y / (1+x^2+y^2)^2
  double denom = (1+M3d_dot_product(point, point));
  denom = denom*denom;
  M3d_vector_mult_const(grad, point, -4/denom);
}



double V_second_deriv(double point[3], double V[3],   double second_deriv[3]) {
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
                     + 2*grad[0]*V[2]*V[0] - grad[2]*V[0]*V[0]
                     + 2*grad[1]*V[2]*V[1] - grad[2]*V[1]*V[1]
                    ) / denom;
}



//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////



const double scaling = 3;
const double D = 1; //0.01;

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
  point[2] = 0.0;
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
  //G_fill_circle(s[0], s[1], 0.5);
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


// gets an object's base color (prior to light model)
int get_color(int onum, double P[2], double rgb[3])
{
  if(color_type[onum] == SIMPLE_COLOR) {
    for(int i=0; i<3; ++i) {
      rgb[i] = color[onum][i];
    }
  } else if(color_type[onum] == TEXTURE_COLOR){
    // convert to object space
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

  for(int n=0; n < 100000; ++n) {

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
      ////G_rgb(1,0,0);
      debug_draw_point(point);
      ////G_wait_key();
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
double light_in_world_space[3] ;
int light_resolution             =  50 ;
int light_obj_res                =  50 ;
int light_interpolation_distance =   5 ; 

int bake_light() {

  if(!do_lightmodel) {
    return 0;
  }

  int onum, u, v, i;
  double tmp_baked_lights[num_objects][light_obj_res][light_obj_res][3];

  // init to 0
  for( onum=0; onum<num_objects; ++onum) {
    for( u=0; u<light_obj_res; ++u) {
      for( v=0; v<light_obj_res; ++v) {
        for( i=0; i<3; ++i) {
          tmp_baked_lights[onum][u][v][i] = 0;
        }
      }
    }
  }


  // cast out a bunch of rays from the light source


  double u_d, v_d, Rtip[3];
  double point[3], V[3];

  int param_u, param_v;
  double param_point[2];

  // cast out rays in a sphere around light source
  for( u=0; u<light_resolution; ++u) {
  //for( u=0; u<50; ++u) {
    //for( v=0; v<1; ++v) {
    printf("starting u-cycle %d of %d\n", u, light_resolution);
    for( v=0; v<light_resolution; ++v) {

      u_d = u*2.0*M_PI / light_resolution ;
      v_d = 0; //v*2*M_PI / light_resolution - M_PI;
      Rtip[0] = cos(u_d)*cos(v_d) ;
      Rtip[1] = sin(u_d)*cos(v_d) ;
      Rtip[2] = sin(v_d) ;

      M3d_vector_add(Rtip, light_in_world_space, Rtip);

      onum = cast_ray(light_in_world_space, Rtip, point, V);
      if(onum == -1) {
        continue;
      }

      // move point to object space
      M3d_mat_mult_pt(point, obinv[onum], point);

      // get parametric coords of point
      to_parametric[onum](point, param_point);
      param_u = param_point[0] * (light_obj_res-1);
      param_v = param_point[1] * (light_obj_res-1);

      M3d_vector_mult_const(tmp_baked_lights[onum][param_u][param_v], V, -1);

      {
        if(debug) {
          // DEBUG
          double p[3],s[2];
          M3d_mat_mult_pt(point, obmat[onum], point);
          M3d_vector_mult_const(p, tmp_baked_lights[onum][param_u][param_v], 1);
          M3d_vector_add(p, point, p);
          coords_to_screen(s, p);
          //G_fill_circle(s[0], s[1], 2);
        }
      }



    }
  }


  // now interpolate
  // for each object, go through the entire parametrized space
  // and interpolate each value
  int interp_u, interp_v;
  double tmp[3];
  for( onum=0; onum<num_objects; ++onum) {
    for( u=0; u<light_obj_res; ++u) {
      for( v=0; v<light_obj_res; ++v) {

        // copy value over
        M3d_vector_copy(baked_lights[onum][u][v], tmp_baked_lights[onum][u][v]);

        // if light already hit this point, skip
        if(M3d_magnitude(baked_lights[onum][u][v]) > 0) {
          continue;
        }

        for(interp_u = u-light_interpolation_distance;
            interp_u < u+light_interpolation_distance;
            ++interp_u) {

          if(interp_u < 0 || interp_u >= light_obj_res) { continue; }

          for(interp_v = v-light_interpolation_distance;
              interp_v < v+light_interpolation_distance;
              ++interp_v) {

            if(interp_v < 0 || interp_v >= light_obj_res) { continue; }

            double dist = (interp_u - u)*(interp_u - u) + (interp_v - v)*(interp_v - v);
            dist /= light_interpolation_distance*light_interpolation_distance;
            M3d_vector_mult_const(tmp,
                                  tmp_baked_lights[onum][interp_u][interp_v],
                                  1.0 / (1.0 + dist));
            M3d_vector_add(baked_lights[onum][u][v], baked_lights[onum][u][v], tmp);

          }
        }

        const double mag = M3d_magnitude(baked_lights[onum][u][v]);
        if(mag > 0) {
          M3d_vector_mult_const(baked_lights[onum][u][v], baked_lights[onum][u][v], 1.0/mag);
        }
      }
    }
  }
}


// light model constants
double AMBIENT      = 0.2 ;
double MAX_DIFFUSE  = 0.5 ;
double SPECPOW      = 50 ;


int Light_Model (double irgb[3],
                 double E[3],
                 double N[3],
                 double L[3],
                 double argb[3])
// s,p,n in eyespace

// irgb == inherent color of object (input to this function)
// E = vector from this poin to eye (input to this function)
// n = normal to the object at p (input to this function)
// L = vector from this point to light (input to this function)
// argb == actual color of object (output of this function)
// globals : AMBIENT, MAX_DIFFUSE, SPECPOW, light_in_eye_space[3]

// return 1 if successful, 0 if error
{



  // this needs to occur BEFORE you possibly jump to LLL below
  double max_ambient_and_diffuse = AMBIENT + MAX_DIFFUSE ;
  double intensity ;

  double len ;

  // normalize all vectors
  if(!M3d_normalize(E, E)) { return 0; }
  if(!M3d_normalize(N, N)) { return 0; }
  if(!M3d_normalize(L, L)) {
    // if L has magnitude 0, then that means it's in shadow, so jump
    intensity = AMBIENT;
    goto LLL;
  }


  /*
  printf("E: "); M3d_print_vector(E);
  printf("N: "); M3d_print_vector(N);
  printf("L: "); M3d_print_vector(L);
  */


  double NdotE = N[0]*E[0] + N[1]*E[1] + N[2]*E[2] ;
  double NdotL = N[0]*L[0] + N[1]*L[1] + N[2]*L[2] ;


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
  else { specular = (1.0 - max_ambient_and_diffuse)*pow(EdotR,SPECPOW) ; }

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

  P[1] = point[1];
  if(point[1] <= -1 || point[1] >= 1) {
    P[1] = (P[1] + 1) / 2;
    P[0] = 0;
    return 0;
  }
  P[0] = (atan2(point[2]/sqrt(1-P[1]*P[1]), point[0]/sqrt(1-P[1]*P[1])) + M_PI) / M_PI / 2;
  P[1] = (P[1] + 1) / 2;
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


  // set color to this object's color:
  double o_color[3];
  get_color(saved_onum, point, o_color);


  // calculate normal
  double normal[3];
  get_normal(normal, saved_onum, point);

  if(do_lightmodel) {

    // flip look vector so it points TO eye FROM point
    M3d_vector_mult_const(look, look, -1);

    // get parametrized coordinates 
    double o_point[3], P[2];
    M3d_mat_mult_pt(o_point, obinv[saved_onum], point);
    to_parametric[saved_onum](o_point, P);

    // DEBUG: draw point
    const int param_x = P[0] * (light_obj_res-1);
    const int param_y = P[1] * (light_obj_res-1);




    if(debug) {
      double p[3], s[2], s2[2];
      const int x = param_x;
      const int y = param_y;
      const double len = 1;

      // normal vector
      M3d_vector_mult_const(p, normal, len);
      M3d_vector_add(p, point, p);
      coords_to_screen(s, p);
      coords_to_screen(s2, point);
      //G_line(s[0], s[1], s2[0], s2[1]);

      // light vector
      M3d_vector_mult_const(p, baked_lights[saved_onum][x][y], len);
      M3d_vector_add(p, point, p);
      coords_to_screen(s, p);
      coords_to_screen(s2, point);
      //G_rgb(1,1,0);
      //G_line(s[0], s[1], s2[0], s2[1]);

      // look vector
      M3d_vector_mult_const(p, look, len);
      M3d_vector_add(p, point, p);
      coords_to_screen(s, p);
      coords_to_screen(s2, point);
      //G_rgb(0,1,0);
      //G_line(s[0], s[1], s2[0], s2[1]);
    }




    Light_Model(o_color,
                look,
                normal,
                baked_lights[saved_onum][param_x][param_y],
                o_color);
  }

  //G_rgb(o_color[0], o_color[1], o_color[2]);
  debug_draw_point(point);
  ////G_line(Rsource[0], Rsource[1], Rtip[0], Rtip[1]);

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
      //G_rgb(o_color[0], o_color[1], o_color[2]);
      double a[2],b[2], debug_tip[3];
      M3d_vector_mult_const(debug_tip, reflection, 1);
      M3d_vector_add(debug_tip, debug_tip, point);

      coords_to_screen(a, point);
      coords_to_screen(b, debug_tip);

      //G_line(a[0], a[1], b[0], b[1]);
      //G_fill_circle(a[0], a[1], 4);
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
      //G_rgb(argb[0], argb[1], argb[2]);
      //G_point(screen_pos[0], screen_pos[1]);

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

  //G_rgb (color[onum][0],color[onum][1],color[onum][2]) ;
  
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
    //G_point(p[0], p[1]);
    //printf("final           : %lf, %lf\n", p[0], p[1]);
    /*
    x = xyz[0] ;
    y = xyz[1] ;
    //G_point(x,y) ;
    */
  }
}


void Draw_plane (int onum)
{
  int n,i ;
  double t, xyz[3] ;
  double x,y ;

  //G_rgb (color[onum][0],color[onum][1],color[onum][2]) ;
  
  n = 1000 ;
  for (i = 0 ; i < n ; i++) {
    t = (i * 1.0 /  n  - 0.5) * 2;
    xyz[0] = t ;
    xyz[1] = 0 ;
    xyz[2] = 0 ;
    M3d_mat_mult_pt(xyz, obmat[onum], xyz) ;
    x = xyz[0] ;
    y = xyz[1] ;
    //G_point(x,y) ;
  }
}

void Draw_hyperbola(int onum)
{
  int n,i ;
  double t, xyz[3] ;
  double x,y ;

  //G_rgb (color[onum][0],color[onum][1],color[onum][2]) ;
  
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
    //G_point(x,y) ;
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


int test01_scene(double eye[3], double coi[3], double up[3]) {
    double Tvlist[100];
    int Tn, Ttypelist[100];
    double m[4][4], mi[4][4];

    num_objects = 0 ;

    //////////////////////////////////////////////////////////////
    const double dist = 4;
    for(int x=0; x<3; ++x) {
      for(int y=1; y<2; ++y) {
        color[num_objects][0] = 1.0 ;
        color[num_objects][1] = 1.0 ; 
        color[num_objects][2] = 1.0 ;

        color[num_objects][x] = 0 ;
        color[num_objects][y] = 0.5 ;

        color_type[num_objects] = SIMPLE_COLOR;
        reflectivity[num_objects] = 0.0;
      
        Tn = 0 ;
        Ttypelist[Tn] = SX ; Tvlist[Tn] =  1        ; Tn++ ;
        Ttypelist[Tn] = SY ; Tvlist[Tn] =  1        ; Tn++ ;
        Ttypelist[Tn] = SZ ; Tvlist[Tn] =  1        ; Tn++ ;
        Ttypelist[Tn] = TX ; Tvlist[Tn] = 1*dist   ; Tn++ ;
        Ttypelist[Tn] = TY ; Tvlist[Tn] =  dist*(x-1)   ; Tn++ ;
        Ttypelist[Tn] = TZ ; Tvlist[Tn] =  dist*(y-1)   ; Tn++ ;
      
        M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
        M3d_copy_mat(obmat[num_objects], m);
        M3d_copy_mat(obinv[num_objects], mi) ;

        SDF[num_objects] = sphere_SDF;
        grad[num_objects] = sphere_grad;
        to_parametric[num_objects] = sphere_to_parametric;

        draw[num_objects] = Draw_ellipsoid; // for 2d
        num_objects++ ; // don't forget to do this
      }
    }
    //////////////////////////////////////////////////////////////

    // place camera
    eye[0] = -1;
    eye[1] = 0;
    eye[2] = 0;

    coi[0] = 0;
    coi[1] = 0;
    coi[2] = 0;

    up[0] = eye[0];
    up[1] = eye[1]+1;
    up[2] = eye[2];

    light_in_world_space[0] = -1;
    light_in_world_space[1] = 0;
    light_in_world_space[2] = 0;


    M3d_view(view_mat, view_inv, eye, coi, up);


}




int test01a()
{

    // view matrix
    double eye[3], coi[3], up[3];
    double origin[3] = {0,0,0};

    test01_scene(eye, coi, up);


    double Rsource[3];
    double Rtip[3];
    double argb[3] ;

    // =====================================================
    // manually make view matrix
    double Tvlist[100];
    int Tn, Ttypelist[100];
    const double scaling = 0.10;

    Tn = 0 ;
    //Ttypelist[Tn] = RY ; Tvlist[Tn] =   90      ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] = -eye[0]   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] = -eye[1]   ; Tn++ ;
    Ttypelist[Tn] = TZ ; Tvlist[Tn] = -eye[2]   ; Tn++ ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  scaling  ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  scaling  ; Tn++ ;
    Ttypelist[Tn] = SZ ; Tvlist[Tn] =  scaling  ; Tn++ ;
  
    M3d_make_movement_sequence_matrix(view_mat, view_inv, Tn, Ttypelist, Tvlist);
    // =====================================================

    debug = 1;
    do_lightmodel = 1;



    int draw_screen() {
      //G_rgb(0,0,0) ;
      //G_clear() ;

      double origin[3] = {0,0,0};
      double s[2], p[2], o[2];
      coords_to_screen(s, Rsource);
      coords_to_screen(p, Rtip);
      coords_to_screen(o, origin);
      //G_rgb(1,0,1) ; G_fill_circle(s[0], s[1], 3) ;
      //G_rgb(1,0,1) ; G_fill_circle(p[0], p[1], 3) ;
      //G_rgb(1,0.5,0.2) ; G_fill_circle(o[0], o[1], 3) ;

      //G_rgb(1,1,0) ; G_line(s[0], s[1], p[0], p[1]) ;
      ////G_rgb(1,0,1) ; G_line(s[0]+50,s[1]-50,  s[0]+50,s[1]+50) ;



      Draw_the_scene() ;
    }



    draw_screen();
    bake_light();
    //G_wait_key();
    
    // first frame, raytrace the whole line
    const int res = 20;
    double colors[SCREEN_HEIGHT][3];

    double p[2];
    M3d_mat_mult_pt(Rsource, view_inv, origin);
    int x_pix = SCREEN_WIDTH;
    if(eye[0] > 0) {
      x_pix = 0;
    }
    for(int y_pix=0; y_pix<SCREEN_HEIGHT; y_pix += res) {

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
    //G_rgb(1,0,1);
    //G_fill_rectangle(x_pix-5, 0, 10, SCREEN_HEIGHT);


    for(int y_pix=0; y_pix<SCREEN_HEIGHT; y_pix += res) {
      //G_rgb(colors[y_pix][0], colors[y_pix][1], colors[y_pix][2]);
      //G_line(x_pix-2, y_pix, x_pix+2, y_pix);
    }
    //G_save_image_to_file("2d_Simple_Raymarcher.xwd") ;

    while(1)
    {



      double p[2];
      //G_wait_click(p);
      if(p[1] < 50) { break; }

      screen_to_coords(Rtip, p);

      //printf("Rsource: "); M3d_print_vector(Rsource);
      //printf("Rtip: "); M3d_print_vector(Rtip);

      draw_screen();

      // draw ray
      ray_to_rgb (Rsource, Rtip, argb) ; 
      //G_rgb(argb[0], argb[1], argb[2]);
      ray_to_rgb (Rsource, Rtip, argb) ; 
    }

    //G_save_image_to_file("2d_Simple_Raymarcher.xwd") ;
}











int test01b()
{
    
    double eye[3], coi[3], up[3];
    test01_scene(eye, coi, up);

    bake_light();

    double Rsource[3];
    double Rtip[3];
    double argb[3] ;
    
    // view matrix
    double origin[3] = {0,0,0};

    ////G_rgb(0, 0.1, 0);
    //G_rgb(0, 0, 0);
    //G_clear();

    // clear the screen
    for(int x=0; x<SCREEN_WIDTH; ++x) {
    	for(int y=0; y<SCREEN_HEIGHT; ++y) {
        set_xwd_map_color(id, x,y, 0.0, 0.0, 0.0);
      }
    }

    double p[2];

    int render_point() {
      screen_to_ray(Rtip, p);

      ray_to_rgb (Rsource, Rtip, argb) ; 
		  int e = set_xwd_map_color(id, p[0], p[1],  argb[0], argb[1], argb[2]);

      //G_rgb(argb[0], argb[1], argb[2]);
      //G_point(p[0], p[1]);
      //G_display_image();
    }

    // first frame, raytrace the whole plane
    const int res = 5;
    double colors[SCREEN_HEIGHT][3];

    M3d_mat_mult_pt(Rsource, view_inv, origin);
    for(int x_pix=300; x_pix<500; x_pix += res) {
      for(int y_pix=250; y_pix<650; y_pix += res) {
    //for(int x_pix = 300; x_pix<500; x_pix += res) {
    //  for(int y_pix = 300; y_pix<500; y_pix += res) {
        p[0] = x_pix;
        p[1] = y_pix;

        render_point();
      }
    }

    printf("finished render\n");
    xwd_map_to_named_xwd_file(id, "9sphere_spherical_lightmodel.xwd");
    //G_save_image_to_file("2d_Simple_Raymarcher.xwd") ;

    while(1) {
      //G_wait_click(p);
      if(p[1] < 50) { break; }

      render_point();
      //printf("Rsource: "); M3d_print_vector(Rsource);
      //printf("Rtip: "); M3d_print_vector(Rtip);


    }

}




//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




int main()
{

  phi = euclid_phi;
  phi_grad = euclid_phi_grad;
  /*
  */

  /*
  phi = sphere_phi;
  phi_grad = sphere_phi_grad;
  */

  //G_init_graphics(800,800);
  id = create_new_xwd_map(SCREEN_WIDTH, SCREEN_HEIGHT);
  if(id == -1) { printf("failure!\n"); exit(0); }
  test01a() ;
}
