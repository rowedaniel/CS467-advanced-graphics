#include "Raymarcher.c"

double dist = 6;
int res = 5;
int gridSize = 40;
int y_grid = -1;
int id;

int fnum = 0;
double move_per_frame = 1;

int build_scene(double eye[3], double coi[3], double up[3]) {
    double Tvlist[100];
    int Tn, Ttypelist[100];
    double m[4][4], mi[4][4];

    num_objects = 0 ;

    //////////////////////////////////////////////////////////////
    if(y_grid == -1) {
      y_grid = gridSize / 2 - 1;
    }
    for(int x=0; x<gridSize; ++x) {
      for(int y=y_grid+1; y<gridSize-y_grid; ++y) {
        color[num_objects][0] =       x * 1.0 / gridSize;
        color[num_objects][1] = 1.0 - x * 1.0 / gridSize;
        color[num_objects][2] =       y * 1.0 / gridSize;

        color_type[num_objects] = SIMPLE_COLOR;
        reflectivity[num_objects] = 0.0;
      
        Tn = 0 ;
        Ttypelist[Tn] = SX ; Tvlist[Tn] =  1        ; Tn++ ;
        Ttypelist[Tn] = SY ; Tvlist[Tn] =  1        ; Tn++ ;
        Ttypelist[Tn] = SZ ; Tvlist[Tn] =  1        ; Tn++ ;
        Ttypelist[Tn] = TX ; Tvlist[Tn] = 1*dist   ; Tn++ ;
        Ttypelist[Tn] = TY ; Tvlist[Tn] = 2.5 * (x-gridSize/2)   ; Tn++ ;
        Ttypelist[Tn] = TZ ; Tvlist[Tn] = 2.5 * (y-gridSize/2)   ; Tn++ ;
      
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
    // build black hole itself (speeds up rendering)
    color[num_objects][0] = 1;
    color[num_objects][1] = 1;
    color[num_objects][2] = 1;

    color_type[num_objects] = SIMPLE_COLOR;
    reflectivity[num_objects] = 0.0;
  
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  0.98       ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  0.98       ; Tn++ ;
    Ttypelist[Tn] = SZ ; Tvlist[Tn] =  0.98       ; Tn++ ;
  
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_copy_mat(obmat[num_objects], m);
    M3d_copy_mat(obinv[num_objects], mi) ;

    SDF[num_objects] = sphere_SDF;
    grad[num_objects] = sphere_grad;
    to_parametric[num_objects] = sphere_to_parametric;

    draw[num_objects] = Draw_ellipsoid; // for 2d
    num_objects++ ; // don't forget to do this
    //////////////////////////////////////////////////////////////

    // place camera
    eye[0] = -20;
    eye[1] = 0;
    eye[2] = move_per_frame * fnum;

    coi[0] = 0;
    coi[1] = 0;
    coi[2] = 0;

    up[0] = eye[0];
    up[1] = eye[1]+1;
    up[2] = eye[2];


    M3d_view(view_mat, view_inv, eye, coi, up);


}



int do_2d()
{

    // view matrix
    double eye[3], coi[3], up[3];
    double origin[3] = {0,0,0};

    build_scene(eye, coi, up);


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



    int draw_screen() {
      //rgb(0,0,0) ;
      //clear() ;

      double origin[3] = {0,0,0};
      double s[2], p[2], o[2];
      coords_to_screen(s, Rsource);
      coords_to_screen(p, Rtip);
      coords_to_screen(o, origin);
      //rgb(1,0,1) ; G_fill_circle(s[0], s[1], 3) ;
      //rgb(1,0,1) ; G_fill_circle(p[0], p[1], 3) ;
      //rgb(1,0.5,0.2) ; G_fill_circle(o[0], o[1], 3) ;

      //rgb(1,1,0) ; G_line(s[0], s[1], p[0], p[1]) ;
      ////rgb(1,0,1) ; G_line(s[0]+50,s[1]-50,  s[0]+50,s[1]+50) ;



      Draw_the_scene() ;
    }



    draw_screen();
    bake_light();
    //wait_key();
    
    // first frame, raytrace the whole line
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
    //rgb(1,0,1);
    //fill_rectangle(x_pix-5, 0, 10, SCREEN_HEIGHT);


    for(int y_pix=0; y_pix<SCREEN_HEIGHT; y_pix += res) {
      //rgb(colors[y_pix][0], colors[y_pix][1], colors[y_pix][2]);
      //line(x_pix-2, y_pix, x_pix+2, y_pix);
    }
    //save_image_to_file("Raymarcher.xwd") ;

    xwd_map_to_named_xwd_file(id, "blackhole_render.xwd");

}











int do_3d()
{
    
    double eye[3], coi[3], up[3];
    build_scene(eye, coi, up);

    bake_light();

    double Rsource[3];
    double Rtip[3];
    double argb[3] ;
    
    // view matrix
    double origin[3] = {0,0,0};

    ////rgb(0, 0.1, 0);
    //rgb(0, 0, 0);
    //clear();

    double p[2];

    int render_point() {
      screen_to_ray(Rtip, p);

      ray_to_rgb (Rsource, Rtip, argb) ; 

      int e = set_xwd_map_color(id, p[0], p[1],  argb[0], argb[1], argb[2]);

    }

    // first frame, raytrace the whole plane
    double colors[SCREEN_HEIGHT][3];

    M3d_mat_mult_pt(Rsource, view_inv, origin);
    for(int x_pix=0; x_pix<SCREEN_WIDTH; x_pix += res) {
      printf("starting x-cycle %d\n", x_pix);
      for(int y_pix=0; y_pix<SCREEN_HEIGHT; y_pix += res) {
    //for(int x_pix = 300; x_pix<500; x_pix += res) {
    //  for(int y_pix = 300; y_pix<500; y_pix += res) {
        p[0] = x_pix;
        p[1] = y_pix;

        render_point();
      }
    }

    printf("finished render\n");
    //save_image_to_file("Raymarcher.xwd") ;

    char filename[100];
    sprintf(filename, "blackhole_render%04d.xwd", fnum);
    xwd_map_to_named_xwd_file(id, filename);


}




//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




