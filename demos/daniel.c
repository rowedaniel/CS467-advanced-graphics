/*
 In this version, we "think" in world space instead
 of in eyespace.  This necessitates only a few changes :

 1.  Instead of passing the light_in_world_space through
     the view matrix to create the light_in_eye_space, 
     instead we merely copy it, leaving the light essentially
     in world space.

 2.  The object matrices, obmat[] and obinv[] no longer
     involve the view matrix and its inverse.

 3.  The endpoints of the initial ray in eyespace, origin 
     and screen_pt, are now passed through the view inverse
     to place them in world space, and the first ray is 
     launched from there.

*/


    origin[0] = 0;
    origin[1] = 0;
    origin[2] = 0;

    // if you wan't screen points that are close to the observer :
    double D,H,HalfWinSize ;
    D = 0.01 ;   // make this even smaller if you want
    H = D*tan_half ;
    HalfWinSize = 0.5*WinSize ;
      // the only real reason for doing this is as a prelude
      // to Daniel's project of making a ray creeper in a distorted
      // space

    
    for(x_pix = 0; x_pix < WinSize; x_pix++){
      for(y_pix = 0; y_pix < WinSize; y_pix++){
	
	
        screen_pt[0] = H*(x_pix - HalfWinSize)/HalfWinSize ;
        screen_pt[1] = H*(y_pix - HalfWinSize)/HalfWinSize ;
        screen_pt[2] = D ;


        // this is what you would have done in eyespace
        // ray_to_rgb(origin, screen_pt, argb,6) ; // 6 levels of reflection	


        // this is what you do if you want to follow rays in worldspace :
        double rstart[3], rend[3] ;
        M3d_mat_mult_pt(rstart,   vi, origin) ; // vi is view matrix inverse
        M3d_mat_mult_pt(rend  ,   vi, screen_pt) ;
        ray_to_rgb(rstart, rend, argb,6) ; // 6 levels of reflection	

	
        G_rgb(argb[0], argb[1], argb[2]) ;
        G_point(x_pix, y_pix) ;

      } // end for y_pix

    } // end for x_pix


