#include"main.h"
#include"sim_pars.h"
#include"periodic.h"

extern sim_pars simu;

void nabla(Triangulation& Tm) {

  //  const bool FME_grad=false; // if true, use FEM approx to the gradient

  // for(F_v_it fv=Tm.finite_vertices_begin();
  //     fv!=Tm.finite_vertices_end();
  //     fv++)    fv->reset_gradp();

  for(F_f_it fc=Tm.finite_faces_begin();
      fc!=Tm.finite_faces_end();
      fc++) {

    std::vector<Vertex_handle> v(3);

    for(int i=0;i<3;i++) 
      v[i]=fc->vertex(i);

    //    FT a=fc->area;

    for(int i0=0; i0< 3 ;i0++) {

      int i1=(i0+1)%3;
      int i2=(i1+1)%3;

      Vector_2 v01=per_vect(v[i0]->point(),v[i1]->point());
      Vector_2 v02=per_vect(v[i0]->point(),v[i2]->point());

      // convention: vectors pointing outwards
      Vector_2 lp1( -v02.y() ,  v02.x() );
      Vector_2 lp2(  v01.y() , -v01.x() );

      Vector_2 lp0 = -lp1 - lp2;

      // IF p_i may be safely substracted from each i node.-

      //  v[i0]->add_gradp(-(p1-p0)/6* lp1);
      //  v[i1]->add_gradp(-(p0-p1)/6* lp0);

      // otherwise.-
      //        v[i0]->add_gradp(-p1/6*lp1);
      //        v[i1]->add_gradp( p0/6*(lp1+lp2));

      v[i0]->nabla(v[i0],  ( 1.0/6.0) * lp1);
      v[i0]->nabla(v[i1],  (-1.0/6.0) * lp1);
      v[i1]->nabla(v[i1],  ( 1.0/6.0) * lp0);
      v[i1]->nabla(v[i0],  (-1.0/6.0) * lp0);

      // cout << 
      // 	v[i0]->idx.val() << " , " <<
      // 	v[i1]->idx.val() << " : " <<
      // 	v[i0]->nabla()[ v[i1] ] << " , " <<
      // 	v[i1]->nabla()[ v[i0] ] << 
      // 	std::endl;

      // cout << 
      // 	v[i0]->idx.val() << " , " <<
      // 	v[i1]->idx.val() << " : " <<
      // 	(-1.0/6.0) * lp1 << " , " <<
      // 	(-1.0/6.0) * lp0  << 
      // 	std::endl;

      // cout << "--" << 	std::endl;


      // otherwise...
      // v[i0]->add_gradp( p0/6*(lp1-lp2));
      // v[i1]->add_gradp( p0/6*(lp1-lp2));

      // v[i0]->add_gradp(-p1/6*lp1);
      // v[i1]->add_gradp(-p1/6*lp1);


      // v[i0]->add_gradp( p0/6*lp1);
      // v[i1]->add_gradp( p0/6*lp1);
      // v[i0]->add_gradp(-p1/6*lp1);
      // v[i1]->add_gradp(-p1/6*lp1);

      if(simu.FEMm()) continue;

      // skip quad calculation.-
      //      if(FME_grad)  continue;

      Fb::v_A coeffs=fc->coeffs[i2];

      for(int n=0;n<coeffs.size();n++) {
	
	FT A=coeffs[n];

	Vertex_handle vn=fc->connects[i2][n];


	v[i0]->nabla( vn   , (-A/24) * ( 2*lp1 + lp0) );
	vn   ->nabla( v[i0], (-A/24) *           lp0  );


	// v[i1]->add_gradp( A/24 * pn * lp1 );
	// vn   ->add_gradp(-A/24 * p1 * lp1 );

	v[i1]->nabla( vn   , (-A/24) * ( lp1 + 2*lp0) );
	vn   ->nabla( v[i1], (-A/24) *   lp1  );

	v[i2]->nabla( vn   , ( A/24) *  lp2  );
	vn   ->nabla( v[i2], (-A/24) *  lp2  );

	//continue;

	//continue;

// self edge
// // // can be safely skipped ????
	for(int m=0;m<coeffs.size();m++) {
	  //	for(int m= n ;m<coeffs.size();m++) {

	  Vertex_handle vm=fc->connects[i2][m];

	  FT Am=coeffs[m];

	  vn->nabla( vm ,  A*Am/60 * lp2);


	}


// next edge

	Fb::v_A coeffs_m=fc->coeffs[i0];

	for(int m=0;m<coeffs_m.size();m++) {

	  FT Am=coeffs_m[m];

	  Vertex_handle vm=fc->connects[i0][m];

	  vn->nabla( vm , -A*Am/120 *  (lp1+2*lp2) );
	  vm->nabla( vn ,  A*Am/120 *  (lp1+2*lp2) );

	} //m

      } //n

    } // vertices of triangle

  } // faces

  // for(F_v_it fv=Tm.finite_vertices_begin();
  //     fv!=Tm.finite_vertices_end();
  //     fv++) {

  //   FT invV;
  //   //    FT invV= 1.0 /fv->vol();

  //   // if(simu.FEM())
  //   if(FME_grad)
  //     invV=1.0 /fv->fem_vol();
  //   else
  //     invV=1.0 /fv->vol();

  //   fv->scale_nabla( invV );

  // }

  return;

}

