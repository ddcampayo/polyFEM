#include"main.h"
#include"periodic.h"
#include"sim_pars.h"
//extern Triangulation T;

//#include <random>

 

extern sim_pars simu;

#include"fields.h"

//extern const FT LL;

//typedef FT (*field)(const FT x,const FT y, bool deriv=false); // pointer to function returning field


FT field_r(const FT x,const FT y, const FT r, bool deriv=false) ;
FT field_Zalesak(const FT x,const FT y) ;
FT field_cos(const FT x, bool deriv=false) ;
FT field_sin(const FT x, bool deriv=false) ;
FT field_linear(const FT x,const FT y, bool deriv=false) ;
FT field_quad(const FT x,const FT y, bool deriv=false) ;
Vector_2 field_rotation(const FT x,const FT y, bool deriv=false) ;
FT field_sin_cos(const FT x,const FT y, bool deriv=false) ;
void alpha_set_mean( Triangulation& T , const FT& mean ) ;



Vector_2 Gresho_v( const FT x, const FT y) {

  const FT tiny = 1e-10;

  const FT rc1 = 0.2;
  const FT rc2 = 0.4;

  FT r2 =  x*x + y*y ;
  FT r  =  std::sqrt( r2 ) ;
  Vector_2  u_theta( Vector_2( -y , x ) / ( r + tiny) ) ;

  FT amp = 0;

  if( r <= rc1 ) amp = 5*r ;
  else if( r <= rc2 )  amp = 2 - 5*r;

  return amp * u_theta ;

}


void set_vels_Gresho(Triangulation& T) {
  
  for(F_v_it vit=T.finite_vertices_begin();
      vit != T.finite_vertices_end();
      vit++) {

    FT x=vit->point().x();
    FT y=vit->point().y();

    vit->U.set( Gresho_v( x , y) ) ;

  }

  return;
}





FT L2_vel_Gresho( Triangulation& T) {

  FT L2=0;
  int nn=0;
  for(F_v_it vit=T.finite_vertices_begin();
      vit != T.finite_vertices_end();
      vit++) {

    FT x=vit->point().x();
    FT y=vit->point().y();

    Vector_2 U0 = Gresho_v( x , y) ;
    Vector_2 U  = vit->U.val();
    L2 += ( U - U0 ).squared_length();
    ++nn;
    
  }

  return L2 / nn;
}




void set_forces_Kolmo(Triangulation& T) {

  for(F_v_it vit=T.finite_vertices_begin();
      vit != T.finite_vertices_end();
      vit++) {

    FT f0 = 4 * M_PI * M_PI ;
    FT y = vit->point().y();

    vit->force.set( f0 * Vector_2(field_cos(y) , 0) );
  }

  return;
}

void set_forces_Kolmo(Triangulation& T, const int nn ) {

  FT f0 = 4 * M_PI * M_PI * simu.v0() ;

  for(F_v_it vit=T.finite_vertices_begin();
      vit != T.finite_vertices_end();
      vit++) {

    FT y = vit->point().y();

    vit->force.set( f0 * Vector_2(field_cos( nn * y ) , 0) );
  }

  return;
}



void set_fields_cos(Triangulation& T) {

  for(F_v_it vit=T.finite_vertices_begin();
      vit != T.finite_vertices_end();
      vit++) {

    FT x=vit->point().x();
    FT y=vit->point().y();

    vit->rold.set( vit->point() );

//    vit->p.set( field_quad(x,y) );
//    vit->p.set( (field_cos(2*x) + field_cos(2*y))/4.0 ) ;
//    vit->alpha.set( field_r(x,y) ) ;
    vit->alpha.set( field_sin(x) * field_sin(y) ) ;
    vit->U.set( Vector_2( field_cos(x) , 0 ));
//    vit->U.set( Vector_2( field_sin(x) , 0 ));
    vit->Uold.set( vit->U.val() );

  }

  return;

}




void set_fields_TG(Triangulation& T) {

  for(F_v_it vit=T.finite_vertices_begin();
      vit != T.finite_vertices_end();
      vit++) {

    FT x=vit->point().x();
    FT y=vit->point().y();

    vit->rold.set( vit->point() );

//    vit->p.set( field_quad(x,y) );
//    vit->p.set( (field_cos(2*x) + field_cos(2*y))/4.0 ) ;
//    vit->alpha.set( field_r(x,y) ) ;
    vit->alpha.set( field_sin(x) * field_sin(y) ) ;
    vit->U.set( Vector_2( field_sin_cos(x,y) , - field_sin_cos(y,x) ));
//    vit->U.set( Vector_2( field_sin(x) , 0 ));
    vit->Uold.set( vit->U.val() );

  }

  return;
}


void set_fields_Zalesak(Triangulation& T) {

  static bool first=true;

  if(first) {
    set_vels_rotating(T);
    first=false;
  }

  for(F_v_it vit=T.finite_vertices_begin();
      vit != T.finite_vertices_end();
      vit++) {

    FT x=vit->point().x();
    FT y=vit->point().y();

    vit->rold.set( vit->point() );
    vit->alpha.set( field_Zalesak(x,y) ) ;
    vit->Uold.set( vit->U.val() );
  }


  return;
}



void fidelity(Triangulation& T, std::ofstream& log_file ) {

  FT dd = 0, dd2 = 0 ;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)    {
    FT al = fv->alpha() ;

    FT x=fv->point().x();
    FT y=fv->point().y();

    FT al0= field_Zalesak(x,y) ;

    dd += (al-al0) * (al-al0) ;
    dd2+= al0 * al0 ;
  }

//   cout << "dd = " << dd << endl ;

  log_file << "  " << std::sqrt( dd / dd2 );

  return;

}




void set_vels_rotating(Triangulation& T) {

  for(F_v_it vit=T.finite_vertices_begin();
      vit != T.finite_vertices_end();
      vit++) {

    FT x=vit->point().x();
    FT y=vit->point().y();

    vit->U.set( field_rotation(x,y) );

  }

  return;
}



Vector_2 field_rotation(const FT x,const FT y, bool deriv) {

  const FT cut   = 1.4;
  const FT width = 0.1;
  
  FT r2= x*x + y*y;

  FT r = std::sqrt(r2);
  if( r > cut )
    return CGAL::NULL_VECTOR ;
  else {
    FT cross = 1; //(1 - std::tanh( ( r - cut) / width  )) / 2;

    return cross * 2 * M_PI * Vector_2( -y , x ) ;

  }

}

// Vector_2 field_rotation(const FT x,const FT y, bool deriv) {

//   const FT radius=0.9;
//   FT dt2 = simu.dt() / 2.0;

//   FT vx= - 2 * M_PI *  y;
//   FT vy =  2 * M_PI *  x;

//   FT x2 = x + dt2 * vx;
//   FT y2 = y + dt2 * vy;

//   if( ( x2 * x2 + y2 * y2 ) > radius * radius )  return CGAL::NULL_VECTOR ;
//   else return 2 * M_PI * Vector_2( -y2 , x2 ) ;

// }


FT field_r(const FT x,const FT y, const FT radius, bool deriv) {

  //  const FT radius=0.5;

  if( (x*x + y*y) < radius*radius )  return 1;
  else return 0;

}


void set_alpha_under_cos( Triangulation& T ) {

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)    {

    FT x=fv->point().x();
    FT y=fv->point().y();

//   FT dt2 = simu.dt() / 2.0;
    FT h= 0.0 * LL * field_cos(x);

    FT val = 0.1;

    if (y < h)
      val *= -1;
    // else
    //   val =  1;

    fv->alpha.set( val );
	
  }

  return;

}

void set_alpha_cos( Triangulation& T  ) {

  FT amp = simu.v0() ;
  
  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++) {
    FT x=fv->point().x();
    FT y=fv->point().y();


    fv->alpha.set( amp * field_cos(x) );
  }
  
  alpha_set_mean( T, 0 );

  return;

}


#include <ctime>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>


void set_alpha_random( Triangulation& T  ) {

  boost::mt19937 randomNumbergenerator( time( 0 ) );

  //const FT limit = 0.1 ;

  FT limit = simu.v0() ;
  
  typedef boost::random::uniform_real_distribution< FT > uniform;

  uniform distribution( -limit , limit );

  boost::variate_generator< boost::mt19937&, uniform >  gen( randomNumbergenerator, distribution );
  
  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++) 
    fv->alpha.set(  gen() );
  //    fv->chempot.set(  gen() );


  alpha_set_mean( T, 0 );

  return;

}


void alpha_set_mean( Triangulation& T , const FT& mean ) {

  FT prev_mean=0;
  int NN=0;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)    {

    prev_mean += fv->alpha.val();
    ++NN ;
  }

  prev_mean /= NN;

  FT offset = prev_mean - mean ;
  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)   
    fv->alpha.set( fv->alpha.val() - offset  );

  return;

}



void zero_mean_v( Triangulation& T ,  const kind::f vectorf ) {

  Vector_2 mean = CGAL::NULL_VECTOR  ;
  int NN=0;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)    {

    mean = mean + fv->vf( vectorf ).val();
    ++NN ;
  }

  mean = (1.0/ FT(NN)) * mean;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)   {
    Vector_2 corr = fv->vf( vectorf ).val() - mean ;
    fv->vf( vectorf ).set( corr );
  }


  // check.-

  mean = CGAL::NULL_VECTOR  ;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)    {

    mean = mean + fv->vf( vectorf ).val();
    ++NN ;
  }

  mean = (1.0/ FT(NN)) * mean;

  cout << "mean = " << mean << endl;
  
  return;

}


void set_alpha_circle( Triangulation& T , const FT& rr ) {

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)    {

    FT x=fv->point().x();
    FT y=fv->point().y();

    FT val = 2 * field_r( x , y , rr ) - 1 ;
    fv->alpha.set( val );

  }


  alpha_set_mean( T ,0 );

  return;

}



FT field_Zalesak(const FT x,const FT y) {

  const FT radius=0.5;
  const FT w= (1.0/6.0) / 2.0 + 0.001;
  const FT top= (1.0/3.0);

  bool in_circle=false;
  bool in_strip=false;

  if( (x*x + y*y) < radius*radius )  in_circle = true ;
  if( (x*x < w*w ) && (y<top))  in_strip = true ;

  if( in_circle && (!in_strip)) return 1;
  else return 0;

}




FT field_cos(const FT x, bool deriv) {
  // if(!deriv)
  //   return x;
  // else
  //   return 0;
 
    return std::cos(2*M_PI*x/LL);

}


FT field_sin(const FT x , bool deriv) {
  return std::sin(2* M_PI * x/LL  );
}


FT field_sin_cos(const FT x,const FT y, bool deriv) {
  return std::sin( 2*M_PI * x/LL  ) * std::cos(2* M_PI * y/LL  );
}

FT field_quad(const FT x,const FT y, bool deriv) {

  if(!deriv)
    return x*x;
  else
    return 2;
}

FT field_linear(const FT x,const FT y, bool deriv) {

  if(!deriv)
    return x+1;
  else
    return 0;
}

