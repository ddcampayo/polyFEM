#include"main.h"
//#include"sim_pars.h"
//#include"linear.h"
//#include"periodic.h"

//extern sim_pars simu;




void areas(Triangulation& T) {

  for(F_f_it fc=T.finite_faces_begin();
      fc!=T.finite_faces_end();
      fc++) {

    Periodic_triangle pt=T.periodic_triangle(fc);
  
    Triangle t=T.triangle(pt);

    fc->area=t.area();

  }

  cout << "Areas computed " << endl;

  return;
}




void integrals(Triangulation& T, std::ofstream& log_file ) {

  FT m=0,  m2 = 0;
  FT  k = 0;
  Vector_2 vm(0,0);
  
  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)    {
    FT vol = fv->vol();

    FT al = fv->alpha() ;

    //    cout << vol << " , " << al << endl;
    m += vol * al;
    m2 += vol * al * al;

    Vector_2 v= fv->U();
    FT v2 = v*v;

    k += v2;

    vm = vm + v;

//   cout << "M: "<<
//    m  << "  " <<
//    m2 << "  " <<
//    vm << "  " <<
//    k << endl ;
  }

  log_file <<
    m  << "  " <<
    m2 << "  " <<
    vm << "  " <<
    k ;
  return;

}



void volumes(Triangulation& T, bool FEM) {

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)    {
    fv->vol.reset();
    //    fv->vol.set(0);
  }

  for(F_f_it fc=T.finite_faces_begin();
      fc!=T.finite_faces_end();
      fc++) {

    std::vector<Vertex_handle> v(3);

    for(int i=0;i<3;i++) 
      v[i]=fc->vertex(i);

    FT a=fc->area;

    for(int i0=0; i0< 3 ;i0++) {
      v[i0]->vol += a / 3.0 ;

      if(FEM) continue;

      Fb::v_A coeffs=fc->coeffs[i0];
      //	Fb::v_vrtcs vvs=fc->connects[i0];

      for(int n=0;n< coeffs.size() ;n++) {
	FT A=coeffs[n];
	fc->connects[i0][n]->vol += A* a/12.0;
      }

    }
  }

  cout << "Volumes computed" << endl;
}





FT kinetic_E( Triangulation& T) {

  FT TT=0;

  for(F_v_it vit=T.finite_vertices_begin();
      vit != T.finite_vertices_end();
      vit++) {

    Vector_2 U  = vit->U.val();
    FT       v  = vit->vol.val();
    TT += v * U.squared_length();
  }

  return TT;
}

