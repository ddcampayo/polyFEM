#include"main.h"
//#include"sim_pars.h"
#include"linear.h"
//#include"periodic.h"

//extern sim_pars simu;

// computes the gradient (vector field) of a scalar field

//typedef Vertex::field field;

void linear::gradient(const kind::f scalarf, const kind::f vectorf, bool do_mass) {

  cout << "gradient of field " << scalarf << " --> " << vectorf << '\n';

  if(lambda_x.size()==0)  fill_lambda();

  VectorXd field = field_to_vctr( scalarf );

  VectorXd grad_x = lambda_x *field;
  VectorXd grad_y = lambda_y *field;

  vctr_to_vfield( grad_x  , vectorf , 0 );
  vctr_to_vfield( grad_y  , vectorf , 1 );

  if(do_mass)
    // full mass inversion.-
    mass_v(vectorf);

  return;
}

VectorXd linear::chempot2(const VectorXd& al ) {

  cout << "chemical potential plus phi of a field " << endl;

  VectorXd al3 = al.array() * al.array() * al.array() ; // To make sure signs are correct

  if(stiff.size()==0) fill_stiff();
  if(mass.size()==0)  fill_mass();

  VectorXd lapl = stiff * al;

  mass_s(lapl);
  
  return  al3 - 0.5 * lapl ;
  //vctr_to_field( al_al3  , chempot  );

}


void linear::chempot(const kind::f scalarf, const kind::f chempot) {

  cout << "chemical potential of field " << scalarf << " --> " << chempot << '\n';

  //  if(lambda_x.size()==0)  fill_lambda();

  VectorXd al  = field_to_vctr( scalarf );

  //  VectorXd al3 = al.array().pow(3);

  VectorXd al3 = al.array() * al.array().square(); // To make sure signs are correct

  VectorXd al_al3 = al3 - al ;  //al.array() * ( -1 + al.array() * al.array() ) ;

  if(stiff.size()==0) fill_stiff();
  //  if(mass.size()==0)  fill_mass();

  VectorXd lapl = stiff * al;

  mass_s(lapl);
  
  vctr_to_field( al_al3  - 0.5 * lapl , chempot  );

  // vctr_to_field( al3 - 0.5 * lapl , chempot  );
  //vctr_to_field( al - 0.5 * lapl , chempot  );
  //  vctr_to_field( al_al3  , chempot  );

  // // model A
  // vctr_to_field( - al3  , chempot  );

  // model B
  //vctr_to_field( al3  , chempot  );

  return;

}



// Laplacian of a scalar field
//template<typename TT>
//void laplacian_Delta(const kind::f ffield, const kind::f gradfield , bool incomplete);

void linear::laplacian_s(const kind::f ffield, const kind::f gradfield  ) {

  cout << "scalar Laplacian of field " << ffield << " --> " << gradfield << '\n';

//  laplacian_Delta( ffield,  gradfield  ) ;
  laplacian_stiff( ffield,  gradfield  ) ;
}

void linear::laplacian_v( const kind::f ffield, const kind::f gradfield  ) {

  cout << "vector Laplacian of field " << ffield << " --> " << gradfield << '\n';

//  laplacian_Delta_v( ffield,  gradfield  ) ;
  laplacian_stiff_v( ffield,  gradfield  ) ;
}

void linear::laplacian_stiff_v(const kind::f ffield, const kind::f gradfield ) {

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    typedef Vertex::scalar_link scalar_link;

    scalar_link stiff=fv->stiff();

    fv->vf(gradfield).reset();

    for(
	scalar_link::iterator nn= stiff.begin();
	nn!=stiff.end(); ++nn) {

      Vertex_handle v=nn->first;

      Vector_2 p=v->vf(ffield).val();

      FT sstiff = -nn->second; // sign fixed

      fv->vf(gradfield) += p * sstiff;

    }

    //    if(!incomplete) 
      // fv->vf(gradfield) /=  ( fv->vol.val() );
  }

  mass_v(gradfield);


}


// Delta-lumped procedure

//template<typename TT>
void linear::laplacian_Delta(const kind::f ffield, const kind::f gradfield   ) {
  //void gradient(int fsf, int fvf ) {

  //  const bool FME_grad=false; // if true, use FEM approx to the gradient

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    //    cout << fv->idx.val() << ":  \n";

    typedef Vertex::scalar_link scalar_link;
    //    typedef std::map<Vertex_handle,Vector_2> vector_link; 
    scalar_link Delta=fv->Delta();

    fv->sf(gradfield).set(0);


    //    fv->get_nabla(nabla) ; // copy

    //    Vb::vector_link& nabla = fv->nabla() ; // alias

    for(
	scalar_link::iterator nn= Delta.begin();
	nn!=Delta.end(); ++nn) {

      Vertex_handle v=nn->first;

      FT p=v->sf(ffield).val();

      FT ddelta = -nn->second; // sign fixed

      fv->sf(gradfield) += p * ddelta;

      // cout <<
      // 	"    " << v->idx.val() << 
      // 	"    " <<  ddelta <<
      // 	"    " << p
      // 	       << std::endl;

    }

    //    if(simu.FEM())  fv->sf(gradfield) /=  6*fv->vol.val();
    //else 

    fv->sf(gradfield) /=  ( fv->vol.val() );
  }

}

void linear::laplacian_Delta_v( const kind::f ffield, const kind::f gradfield   ) {

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    typedef Vertex::scalar_link scalar_link;

    scalar_link Delta=fv->Delta();

    fv->vf(gradfield).reset();

    for(
	scalar_link::iterator nn= Delta.begin();
	nn!=Delta.end(); ++nn) {

      Vertex_handle v=nn->first;

      Vector_2 p=v->vf(ffield).val();

      FT ddelta = -nn->second; // sign fixed

      fv->vf(gradfield) += p * ddelta;

    }

    //    if(!incomplete) 
      // fv->vf(gradfield) /=  ( fv->vol.val() );
  }

  mass_v(gradfield);


}



// another candidate for templating

void linear::laplacian_stiff( const kind::f ffield, const kind::f gradfield   ) {

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    typedef Vertex::scalar_link scalar_link;

    scalar_link stiff=fv->stiff();

    fv->sf(gradfield).reset();

    for(
	scalar_link::iterator nn= stiff.begin();
	nn!=stiff.end(); ++nn) {

      Vertex_handle v=nn->first;

      FT p=v->sf(ffield).val();

      FT ddelta = -nn->second; // sign fixed

      fv->sf(gradfield) += p * ddelta;

      // cout <<
      // 	"    " << v->idx.val() << 
      // 	"    " <<  ddelta <<
      // 	"    " << p
      // 	       << std::endl;


    }

    //    if(simu.FEM())  fv->sf(gradfield) /=  6*fv->vol.val();
    //else 

    //    stiff-lumped procedure
    //    fv->sf(gradfield) /=  ( fv->vol.val() );
  }
  
  mass_s(gradfield);
  
}







void linear::PPE(const kind::f velocity , const FT dt,  const kind::f pressure ) {

//  linear linear_problem;

//  laplace_div( kind::USTAR , kind::DIVU , kind::P);
  laplace_div( velocity , dt , kind::DIVU , pressure );

  // // divide by dt:

  // for(F_v_it fv=T.finite_vertices_begin();
  //     fv!=T.finite_vertices_end();
  //     fv++)  fv->sf(pressure) /= dt;

  return;

}
