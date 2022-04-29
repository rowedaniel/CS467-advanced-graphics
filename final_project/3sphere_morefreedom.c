#include "Raymarcher.c"

int res = 5;

int build_scene(double eye[3], double coi[3], double up[3]) {
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



    draw_screen();
    bake_light();
    //G_wait_key();
    
    // first frame, raytrace the whole line
    double colors[SCREEN_HEIGHT][3];

    double p[2], s[2];

    while(1)
    {


      G_wait_click(s);
      if(s[1] < 50) { break; }
      G_wait_click(p);
      if(p[1] < 50) { break; }

      screen_to_coords(Rtip, p);
      screen_to_coords(Rsource, s);

      //printf("Rsource: "); M3d_print_vector(Rsource);
      //printf("Rtip: "); M3d_print_vector(Rtip);

      draw_screen();

      // draw ray
      ray_to_rgb (Rsource, Rtip, argb) ; 
      if(argb[0] == 0 && argb[1] == 0 && argb[2] == 0) {
        G_rgb(1, 0, 1);
      } else {
        G_rgb(argb[0], argb[1], argb[2]);
      }
      ray_to_rgb (Rsource, Rtip, argb) ; 
    }

    G_save_image_to_file("Raymarcher.xwd") ;
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

    //G_rgb(0, 0.1, 0);
    G_rgb(0, 0, 0);
    G_clear();

    double p[2];

    int render_point() {
      screen_to_ray(Rtip, p);

      ray_to_rgb (Rsource, Rtip, argb) ; 

      G_rgb(argb[0], argb[1], argb[2]);
      G_point(p[0], p[1]);
      G_display_image();
    }

    // first frame, raytrace the whole plane
    double colors[SCREEN_HEIGHT][3];

    M3d_mat_mult_pt(Rsource, view_inv, origin);
    for(int x_pix=300; x_pix<500; x_pix += res) {
      for(int y_pix=0; y_pix<SCREEN_HEIGHT; y_pix += res) {
    //for(int x_pix = 300; x_pix<500; x_pix += res) {
    //  for(int y_pix = 300; y_pix<500; y_pix += res) {
        p[0] = x_pix;
        p[1] = y_pix;

        render_point();
      }
    }

    printf("finished render\n");
    G_save_image_to_file("Raymarcher.xwd") ;

    while(1) {
      G_wait_click(p);
      if(p[1] < 50) { break; }

      render_point();
      //printf("Rsource: "); M3d_print_vector(Rsource);
      //printf("Rtip: "); M3d_print_vector(Rtip);


    }

}




//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




