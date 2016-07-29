#include"main.h"
#include"sim_pars.h"
#include"periodic.h"

extern sim_pars simu;

struct data_kept {
  int idx;
  Point pos;
  Point rold;
  Vector_2 Uold;
  Vector_2 laplU;
  Vector_2 U;
  FT divU;
  Vector_2 Ustar;
  FT alpha;
  FT p;
  Vector_2 gradp;
  FT vol;

  data_kept(const F_v_it fv) {
    idx = fv->idx();
    U = fv->U();
    divU = fv->divU();
    rold = fv->rold();
    Uold = fv->Uold();
    Ustar = fv->Ustar();
    laplU = fv->laplU();
    alpha=fv->alpha();
    p= fv->p();
    gradp= fv->gradp();
    vol= fv->vol();
  }

  void restore(Vertex_handle fv) {
    fv->idx.set( idx );
    fv->U.set( U );
    fv->divU.set( divU );
    fv->rold.set( rold );
    fv->Uold.set( Uold );
    fv->Ustar.set( Ustar );
    fv->laplU.set( laplU );
    fv->alpha.set( alpha );
    fv->p.set( p );
    fv->gradp.set( gradp );
    fv->vol.set(vol);
  }

};

FT move(FT dt) {

  vector<data_kept> prev;

  FT dd2=0;

  for(F_v_it fv=Tp.finite_vertices_begin();
      fv!=Tp.finite_vertices_end();
      fv++) {
    data_kept data(fv);

    Vector_2  vel = fv->U();

    Vector_2 disp = dt * vel;

    Periodic_point rr=Tp.periodic_point(fv);

    Point rnow=Tp.point(rr); // current point

    Point r0=fv->rold(); // starting point

    Point rnew= r0 + disp;

    Vector_2 disp2 = per_vect(rnew,rnow);
    dd2 += sqrt(disp2.squared_length())/simu.h();

    //    cout << "New position: " << r0 ;

    data.pos = per_point( rnew );

    //    cout << " ---> " << data.pos  << endl ;

    prev.push_back (data);

  }

//  cout << "relative displacement " << sqrt(dd2)/simu.no_of_points()/simu.h()  << endl ;
  dd2 /= simu.no_of_particles();

//  cout << "relative displacement " << dd2 << endl ;

  Tp.clear(); // clears the triangulation !!

  for(vector<data_kept>::iterator data=prev.begin();
      data!=prev.end();
      data++) {

    // cout << "Inserting back at " << data->pos << endl ;

    Vertex_handle fv=Tp.insert(data->pos);

    data->restore(fv);

    // return info to vertices


  }

  Tp.convert_to_1_sheeted_covering();

  return dd2;
}


void move_info(void) {

  for(F_v_it fv=Tp.finite_vertices_begin();
      fv!=Tp.finite_vertices_end();
      fv++) {

    Periodic_point rr=Tp.periodic_point(fv);

    Point r=Tp.point(rr); // current point

    fv->rold.set(r);
    fv->Uold.set(fv->U());
    fv->Ustar.set(fv->U());
    //    fv->p.set(fv->p() +  fv->pstar() );
  }


  for(F_v_it fv=Tm.finite_vertices_begin();
      fv!=Tm.finite_vertices_end();
      fv++) {
    fv->Uold.set(fv->U());
  }

  return;
}


void u_new(FT dt) {

  for(F_v_it fv=Tm.finite_vertices_begin();
      fv!=Tm.finite_vertices_end();
      fv++) {

//  for(F_v_it fv=Tp.finite_vertices_begin();
//      fv!=Tp.finite_vertices_end();
//      fv++) {

    Vector_2 Ustar = fv->Ustar.val() ;
    Vector_2 gradp = fv->gradp.val() ;
    Vector_2 U = Ustar - dt * gradp;

    // relaxation mixing .-
    FT alpha=simu.alpha();
    Vector_2 U0=fv->U() ;
    fv->U.set( alpha*U0+ (1-alpha)*U );

  }

  return;

}


void u_star(FT dt , bool semi ) { 

  for(F_v_it fv=Tm.finite_vertices_begin();
      fv!=Tm.finite_vertices_end();
      fv++) {

    Vector_2 U0 = fv->Uold() ;

    Vector_2 f;

    if(semi)
      f = simu.mu() * fv->laplU() - fv->gradp();
    else
      f= simu.mu() * fv->laplU() ; // Include all forces here, but for the pressure grad!

    Vector_2 Ustar = U0 + dt * f;

    fv->Ustar.set(Ustar) ;

  }

  return;

}
