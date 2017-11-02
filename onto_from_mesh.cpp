#include"main.h"
#include"sim_pars.h"
#include"periodic.h"
#include"linear.h"

extern sim_pars simu;

FT values_at(const Triangulation& T, const Point& p0, const bool FEM,
               const kind::f scalarf ) ;

Vector_2 values_at_v(const Triangulation& T, const Point& p0, const bool FEM,
               const kind::f vectorf ) ;


FT vect_prod(const Vector_2& v1 , const Vector_2& v2 ) {
  return v1.x() * v2.y() - v1.y() * v2.x() ;
}

#include"onto_from_mesh.h"

// local util functions.-

void lumped_full_common(Triangulation& Tfrom, Triangulation& Tto, const kind::f scalarf);
void lumped_full_common_v(Triangulation& Tfrom, Triangulation& Tto, const kind::f vectorf);

void lumped_full_common_vertices(Triangulation& Tfrom, Triangulation& Tto, const kind::f scalarf);
void lumped_full_common_vertices_v(Triangulation& Tfrom, Triangulation& Tto, const kind::f vectorf);

//void lumped_full_common_interp(void) ;

void lumped_full_common_centroid(const kind::f scalarf);
void lumped_full_common_centroid_v(const kind::f vectorf);

void FEM_hs(const Face_handle& fc, const Point& p0,
	    std::vector<Vertex_handle>& v,
	    std::vector<FT>& hh    );



void reset_v(Triangulation& T, const kind::f vectorf) {

  for(F_v_it vit=T.vertices_begin();
      vit != T.vertices_end();
      vit++) {
    vit->vf(vectorf).reset();
  }

  return;
}

void incr_v(Triangulation& T, const kind::f vectorfrom,const kind::f vectorincr, const kind::f vectorto) {

  for(F_v_it vit=T.vertices_begin();
      vit != T.vertices_end();
      vit++) {
    Vector_2  v_from    = vit->vf(vectorfrom).val();
    Vector_2  Delta_v = vit->vf(vectorincr).val();
    
    vit->vf(vectorto).set(  v_from + Delta_v  );
  }

  return;
}




void from_mesh_full( Triangulation& Tfrom, Triangulation& Tto, linear& algebra_p, const kind::f scalarf) {

  lumped_full_common( Tfrom , Tto , scalarf);

  algebra_p.mass_s(scalarf);

  return;
}


void from_mesh_full_v( Triangulation& Tfrom, Triangulation& Tto, linear& algebra_p, const kind::f vectorf) {

  lumped_full_common_v( Tfrom, Tto , vectorf);

  algebra_p.mass_v(vectorf);

  return;
}


void from_mesh_lumped( Triangulation& Tfrom, Triangulation& Tto, const kind::f scalarf) {

  lumped_full_common(  Tfrom, Tto , scalarf);

  for(F_v_it vit=Tto.vertices_begin();
      vit != Tto.vertices_end();
      vit++) {
    FT V=vit->vol();
    vit->sf(scalarf) /= V ;
  }

  return;
}


void from_mesh_lumped_v(  Triangulation& Tfrom, Triangulation& Tto, const kind::f vectorf) {

  lumped_full_common_v( Tfrom, Tto , vectorf);

  for(F_v_it vit=Tto.vertices_begin();
      vit != Tto.vertices_end();
      vit++) {
    FT V=vit->vol();
    vit->vf(vectorf) /= V ;
  }

  return;
}



void onto_mesh_full( Triangulation& Tfrom, Triangulation& Tto, linear& algebra, const kind::f scalarf) {

  lumped_full_common( Tfrom, Tto , scalarf);

  algebra.mass_s(scalarf);

  return;
}



void onto_mesh_full_v( Triangulation& Tfrom, Triangulation& Tto, linear& algebra, const kind::f vectorf) {

  lumped_full_common_v( Tfrom, Tto , vectorf);

  algebra.mass_v(vectorf);

  return;
}



void lumped_full_common(Triangulation& Tfrom, Triangulation& Tto, const kind::f scalarf) {
  lumped_full_common_vertices( Tfrom, Tto, scalarf );
  //lumped_full_common_centroid(scalarf);
   //   lumped_full_common_interp();
   return;
}



void lumped_full_common_v(Triangulation& Tfrom, Triangulation& Tto,const kind::f vectorf) {
  lumped_full_common_vertices_v( Tfrom, Tto, vectorf );
  //   lumped_full_common_centroid_v(vectorf);
   //   lumped_full_common_interp();
   return;
}


void onto_mesh_delta( Triangulation& Tfrom, Triangulation& Tto, const kind::f scalarf) {

  for(F_v_it vit=Tto.vertices_begin();
      vit != Tto.vertices_end();
      vit++) {

    Point pm =  vit->point();
//    cout << pm << endl ;

    FT ss=values_at( Tfrom , pm , simu.FEMp(), scalarf );

    vit->sf( scalarf).set( ss );
  }

  return;

}

void onto_mesh_delta_v(Triangulation& Tfrom, Triangulation& Tto, const kind::f vectorf) {

  for(F_v_it vit=Tto.vertices_begin();
      vit != Tto.vertices_end();
      vit++) {

    Point pm =  vit->point();

    Vector_2 vv=values_at_v(Tfrom , pm , simu.FEMp(), vectorf );

    vit->vf( vectorf ).set( vv );
  }

  return;

}


void FEM_hs(const Face_handle& fc, const Point& p0,
	    std::vector<Vertex_handle>& v,
	    std::vector<FT>& hh    ) {  

  for(int i=0; i < 3 ; ++i)    v[i] = fc->vertex(i);

  std::vector<Point> p(3);

  for(int i=0; i < 3 ; ++i)    p[i] = v[i]->point();

  Vector_2 v01 = per_vect( p[0] , p[1] );
  Vector_2 v02 = per_vect( p[0] , p[2] );

  FT A = fc->area;

  Vector_2 rr = per_vect( p[0] , p0 );

  hh[1] = vect_prod( rr , v02 ) / (2*A);

  hh[2] = vect_prod( v01 , rr ) / (2*A);

  hh[0] = 1 - hh[1] - hh[2];

  return;
}

		       
FT values_at(const Triangulation& T, const Point& p0, const bool FEM,
	       const kind::f scalarf) {

  std::vector<Vertex_handle> v(3);
  std::vector<FT> hh(3);

  //  cout << "Locating p0 = " << p0 << endl ;
  Face_handle fc=T.locate( p0 );
  //  cout << "   located" << endl;

  FEM_hs( fc , p0, v, hh);

  FT ss=0;

  for(int i0 = 0; i0 < 3 ; ++i0) {

    ss += v[i0]->sf( scalarf ).val() * hh[i0];

    if(FEM) continue;
      
    typedef std::vector<Vertex_handle> v_vrtcs;
    typedef std::vector<FT> v_A;

    //      cout << i0  << endl;

    v_vrtcs  cvs=fc->connects[i0]; // cv: connected vertices
    v_A coeffs=fc->coeffs[i0];

    int i1= (i0+1)%3;
    int i2= (i0+2)%3;

    FT overlap = hh[i1] * hh[i2] ; // * A/12 ;

    for(int n=0;n<cvs.size();n++) {

      FT mm;
	
      FT AA=coeffs[n];
      Vertex_handle cv=cvs[n];

//	cout << n << " " << mm << endl;

      mm=  AA * overlap;

      ss += cv->sf( scalarf ).val() * mm;
    }
  }

  return ss;
}


Vector_2 values_at_v(const Triangulation& T, const Point& p0, const bool FEM,
	       const kind::f vectorf) {

  std::vector<Vertex_handle> v(3);
  std::vector<FT> hh(3);

  //  cout << "Locating p0 = " << p0 ;
  Face_handle fc=T.locate( p0 );
  //  cout << "   located" << endl;

  FEM_hs( fc , p0, v, hh);

  Vector_2 vv(0,0);

  for(int i0 = 0; i0 < 3 ; ++i0) {

    vv = vv + hh[i0] * v[i0]->vf( vectorf ).val() ;

    if(FEM) continue;
      
    typedef std::vector<Vertex_handle> v_vrtcs;
    typedef std::vector<FT> v_A;

    //      cout << i0  << endl;

    v_vrtcs  cvs=fc->connects[i0]; // cv: connected vertices
    v_A coeffs=fc->coeffs[i0];

    int i1= (i0+1)%3;
    int i2= (i0+2)%3;

    FT overlap = hh[i1] * hh[i2] ; // * A/12 ;

    for(int n=0;n<cvs.size();n++) {

      FT mm;
	
      FT AA=coeffs[n];
      Vertex_handle cv=cvs[n];

//	cout << n << " " << mm << endl;

      mm=  AA * overlap;

      vv = vv + mm * cv->vf( vectorf ).val() ;
    }
  }

  return vv;
}



void onto_mesh_lumped(Triangulation& Tfrom, Triangulation& Tto, const kind::f scalarf) {

  lumped_full_common( Tfrom , Tto ,  scalarf );

  for(F_v_it vit=Tto.vertices_begin();
      vit != Tto.vertices_end();
      vit++) {
    FT V=vit->vol();

    vit->sf(scalarf) /= V ;

    //    cout << vit->p() << endl;
  }

  return;
}







void lumped_full_common_vertices( Triangulation& Tfrom, Triangulation& Tto, const kind::f scalarf ) {

  // Linear approximation for the integrals,
  // as opposed to quadratic
  // May be tied to simu.FEMp() or not!
  
  bool linear = false ; // = simu.FEMp();
  
  for(F_v_it vit=Tto.vertices_begin();
      vit != Tto.vertices_end();
      vit++) 
    vit->sf(scalarf).reset();

  for(F_f_it fit=Tto.faces_begin();
      fit != Tto.faces_end();
      fit++) {

    vector<FT> field_at_v( 3 );

    for(int i0 = 0; i0 < 3 ; ++i0) {
      Vertex_handle vm = fit->vertex(i0);

      Point pm=per_point( vm->point() );

      //      cout << "Finding value at " << pm << endl;
      field_at_v[i0] =  values_at(Tfrom , pm , simu.FEMp(), scalarf );

    }

    // Field opposite vertex i:
    vector<FT> field_op_v( 3 );

    //        if( !simu.FEMp() ) {
    if( !linear ) {
      for(int i0 = 0; i0 < 3 ; ++i0) {
	int i1 = (i0+1) % 3;
	int i2 = (i1+1) % 3;

	Vertex_handle v1 = fit->vertex(i1);
	Vertex_handle v2 = fit->vertex(i2);

	Point p1=per_point( v1->point() );
	Point p2=per_point( v2->point() );

	Vector_2 v12=per_vect(p1,p2);

	Point pm=per_point( p1 + v12 / 2.0 );

	field_op_v[i0] = values_at( Tfrom , pm , simu.FEMp(), scalarf );

      }

    }


    // Do the area integrals of the field values:

    FT a=fit->area;
    //    FT vol = a / 3.0;

    std::vector<FT>         pm(3);

    //    if( !simu.FEMp() ) 
    if( !linear )
      for(int i0 = 0; i0 < 3 ; ++i0) {
	int i1 = (i0+1) % 3;
	int i2 = (i1+1) % 3;

	pm[i0] =	 4 * field_op_v[i0] - 2 * field_at_v[i1]  - 2 * field_at_v[i2] ;

      }

    for(int i0 = 0; i0 < 3 ; ++i0) {

      Vertex_handle vm = fit->vertex(i0);

      int i1 = (i0+1) % 3;
      int i2 = (i1+1) % 3;

      vm->sf(scalarf) += ( 2 * field_at_v[i0] + field_at_v[i1] + field_at_v[i2] ) * a / 12.0 ;

      FT p0m, p1m, p2m;
      Vector_2 U0m, U1m, U2m;
      Vector_2 Uold0m, Uold1m, Uold2m;

      //      if( !simu.FEMp() )
      if( !linear ) 
	vm->sf(scalarf) += ( pm[i0] + 2 * pm[i1] + 2 * pm[i2] ) * a / 60.0 ;

      if(simu.FEMm()) continue;

      typedef std::vector<Vertex_handle> v_vrtcs;
      typedef std::vector<FT> v_A;

      v_vrtcs  cvs=fit->connects[i0]; // cv: connected vertices
      v_A coeffs=fit->coeffs[i0];

      //      FT vol2 = a / 12.0;

      for(int n=0;n<cvs.size();n++) {
	
	FT AA=coeffs[n];
	Vertex_handle cv=cvs[n];

	cv->sf(scalarf) += ( field_at_v[i0] + 2 * field_at_v[i1] + 2* field_at_v[i2] ) * AA * a / 60.0 ;

	//	if( !simu.FEMp() )
	  if( !linear ) 
	    cv->sf(scalarf) += ( 2 * pm[i0] +  pm[i1] +  pm[i2] ) * AA * a / 180.0 ;

      }

    }
  }

  return;
}




void lumped_full_common_vertices_v( Triangulation& Tfrom, Triangulation& Tto, const kind::f vectorf ) {

    
  bool linear = false ; // = simu.FEMp();

  for(F_v_it vit=Tto.vertices_begin();
      vit != Tto.vertices_end();
      vit++) 
    vit->vf(vectorf).reset();

  for(F_f_it fit=Tto.faces_begin();
      fit != Tto.faces_end();
      fit++) {

    vector<Vector_2> field_at_v( 3 );

    for(int i0 = 0; i0 < 3 ; ++i0) {
      Vertex_handle vm = fit->vertex(i0);

      Point pm=per_point( vm->point() );

      field_at_v[i0] =  values_at_v(Tfrom , pm , simu.FEMp(), vectorf );

    }

    vector<Vector_2> field_op_v( 3 );

    //        if( !simu.FEMp() ) {
    if( !linear ) {
      for(int i0 = 0; i0 < 3 ; ++i0) {
	int i1 = (i0+1) % 3;
	int i2 = (i1+1) % 3;

	Vertex_handle v1 = fit->vertex(i1);
	Vertex_handle v2 = fit->vertex(i2);

	Point p1=per_point( v1->point() );
	Point p2=per_point( v2->point() );

	Vector_2 v12=per_vect(p1,p2);

	Point pm=per_point( p1 + v12 / 2.0 );

	field_op_v[i0] = values_at_v( Tfrom , pm , simu.FEMp(), vectorf );

      }

    }


    // Do the area integrals of the field values:

    FT a=fit->area;
    //    FT vol = a / 3.0;

    std::vector<Vector_2>         pm(3);


    //    if( !simu.FEMp() ) 
    if( !linear )
      for(int i0 = 0; i0 < 3 ; ++i0) {
	int i1 = (i0+1) % 3;
	int i2 = (i1+1) % 3;

	pm[i0] =	 4 * field_op_v[i0] - 2 * field_at_v[i1]  - 2 * field_at_v[i2] ;

      }

    for(int i0 = 0; i0 < 3 ; ++i0) {

      Vertex_handle vm = fit->vertex(i0);

      int i1 = (i0+1) % 3;
      int i2 = (i1+1) % 3;

      vm->vf(vectorf) += ( 2 * field_at_v[i0] + field_at_v[i1] + field_at_v[i2] ) * a / 12.0 ;

      FT p0m, p1m, p2m;
      Vector_2 U0m, U1m, U2m;
      Vector_2 Uold0m, Uold1m, Uold2m;

      
    //    if( !simu.FEMp() ) 
      if( !linear )
	vm->vf(vectorf) += ( pm[i0] + 2 * pm[i1] + 2 * pm[i2] ) * a / 60.0 ;

      if(simu.FEMm()) continue;

      typedef std::vector<Vertex_handle> v_vrtcs;
      typedef std::vector<FT> v_A;

      v_vrtcs  cvs=fit->connects[i0]; // cv: connected vertices
      v_A coeffs=fit->coeffs[i0];

      //      FT vol2 = a / 12.0;

      for(int n=0;n<cvs.size();n++) {
	
	FT AA=coeffs[n];
	Vertex_handle cv=cvs[n];

	cv->vf(vectorf) += ( field_at_v[i0] + 2 * field_at_v[i1] + 2* field_at_v[i2] ) * AA * a / 60.0 ;

	if( !simu.FEMp() )
	  cv->vf(vectorf) += ( 2 * pm[i0] +  pm[i1] +  pm[i2] ) * AA * a / 180.0 ;

      }

    }
  }

  return;
}



void from_mesh(Triangulation& Tfrom, Triangulation& Tto, const kind::f scalarf) {

  for(F_v_it vit=Tto.vertices_begin();
      vit != Tto.vertices_end();
      vit++) {

    Point pm =  vit->point();

    vit->sf( scalarf).set( values_at(Tfrom , pm , simu.FEMm(), scalarf) );

  }

  return;

}




void from_mesh_v(Triangulation& Tfrom, Triangulation& Tto, const kind::f vectorf) {

  for(F_v_it vit=Tto.vertices_begin();
      vit != Tto.vertices_end();
      vit++) {

    Point pm =  vit->point();

    vit->vf( vectorf ).set( values_at_v(Tfrom , pm , simu.FEMm(), vectorf) );

  }

  return;

}




void flip_volumes(Triangulation& Tpart, Triangulation& Tmesh, bool FEM) {

  for(F_v_it fv=Tmesh.finite_vertices_begin();
      fv!=Tmesh.finite_vertices_end();
      fv++) {
    //    fv->fvol.reset();

    Point pm = fv->point();

    Face_handle fc=Tpart.locate( pm );

    std::vector<Vertex_handle> v(3);
    std::vector<FT> hh(3);

    FEM_hs( fc , pm, v, hh);

    FT fvol=0;

    for(int i0=0; i0< 3 ;i0++)
      fvol += hh[i0] * v[i0]->vol();

    fv->fvol.set( fvol ) ;

  }

  cout << "FLIP volumes computed" << endl;
}




void onto_mesh_flip(Triangulation& Tpart, Triangulation& Tmesh, bool FEM, const kind::f scalarf ) {


  for(F_v_it fv=Tmesh.finite_vertices_begin();
      fv!=Tmesh.finite_vertices_end();
      fv++) {
    fv->sf(scalarf).reset();

    Point pm = fv->point();

    Face_handle fc=Tpart.locate( pm );

    std::vector<Vertex_handle> v(3);
    std::vector<FT> hh(3);

    FEM_hs( fc , pm, v, hh);

    for(int i0=0; i0< 3 ;i0++)
      fv->sf(scalarf) += hh[i0] * v[i0]->sf(scalarf).val() * v[i0]->vol();

    fv->sf(scalarf) /= fv->fvol.val();


  }

}


void onto_mesh_flip_v(Triangulation& Tpart, Triangulation& Tmesh, bool FEM, const kind::f vectorf ) {

  for(F_v_it fv=Tmesh.finite_vertices_begin();
      fv!=Tmesh.finite_vertices_end();
      fv++) {
    fv->vf(vectorf).reset();

    Point pm = fv->point();

    Face_handle fc=Tpart.locate( pm );

    std::vector<Vertex_handle> v(3);
    std::vector<FT> hh(3);

    FEM_hs( fc , pm, v, hh);

    for(int i0=0; i0< 3 ;i0++)
      fv->vf(vectorf) += hh[i0] * v[i0]->vf(vectorf).val() * v[i0]->vol();

    fv->vf(vectorf) /= fv->fvol.val();

  }

}

