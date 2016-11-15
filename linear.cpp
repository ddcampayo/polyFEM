#include"main.h"
#include"linear.h"
#include"sim_pars.h"
//#include"periodic.h"

#define TOL 1e-12

// quite dirty
//#define LOAD

extern sim_pars simu;

void linear::save_matrices(void){ 

  if(stiff.size()==0) fill_stiff();
  if(mass.size()==0)  fill_mass();
  if(lambda_x.size()==0)  fill_lambda();

  std::cout << "Saving matrices" << std::endl;

  saveMarket(stiff, "stiff.mtx");
  saveMarket(mass, "mass.mtx");
  saveMarket(stiffp1, "stiffp1.mtx");
  saveMarket(lambda_x, "lambda_x.mtx");
  saveMarket(lambda_y, "lambda_y.mtx");

  return;

}

void linear::load_matrices(void){ 

  loadMarket(stiff, "stiff.mtx");
  loadMarket(stiffp1, "stiffp1.mtx");
  loadMarket(mass, "mass.mtx");

  return;

}



void linear::fill_stiff(void){ 

#ifdef LOAD
  std::cout << " Reading stiffness matrix" << std::endl;
  loadMarket(stiff, "stiff.mtx");
  loadMarket(stiffp1, "stiffp1.mtx");
#else

  std::cout << " Filling stiffness matrix" << std::endl;

  //  int n=simu.no_of_points();

  std::vector<triplet> aa , bb ;            // list of non-zeros coefficients

  int N=0;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    ++N;

    int vi=fv->idx.val();

    typedef Vertex::scalar_link scalar_link;
    //    typedef std::map<Vertex_handle,Vector_2> vector_link; 
    scalar_link stiffs=fv->stiff();

    for(
	scalar_link::iterator nn= stiffs.begin();
	nn!=stiffs.end(); ++nn) {

      int vj=nn->first->idx.val();
	
      FT ddelta = -nn->second; // sign fixed here !

      //      cout << vi << " , " << vj << "  " << ddelta << endl;

      aa.push_back( triplet(vi,vj,  ddelta ));

      if( (vi!=0) && (vj!=0) ) bb.push_back( triplet(vi - 1 , vj -1 ,  -ddelta ));



    }
  }

  stiff.resize(N,N);
  stiff.setFromTriplets(aa.begin(), aa.end());
  std::cout << " Filled stiffness  matrix" << std::endl;
  cout << "matrix size " << stiff.rows() << " x " << stiff.cols() << endl;

  // F_v_it fv=T.finite_vertices_begin();
  // int i0=fv->idx.val();
  // aa.push_back( triplet( N , i0 , 1 ) );  // anchor particle 0 to have fixed pressure
  //    //aa.push_back( triplet( N , 0 , 1 ) );  // anchor particle 0 to have fixed pressure

  // stiffp1.resize(N+1,N+1);
  // stiffp1.setFromTriplets(aa.begin(), aa.end());

  stiffp1.resize(N-1,N-1);
  stiffp1.setFromTriplets(bb.begin(), bb.end());

  std::cout << " Filled stiffness plus 1 matrix" << std::endl;
  cout << "matrix size " << stiffp1.rows() << " x " << stiffp1.cols() << endl;
#endif

  // Only non-direct iterative solvers
#ifndef DIRECT_SOLVER
  solver_stiffp1.setTolerance( TOL );
  solver_stiffp1.setMaxIterations( 40 * N );
#endif

  solver_stiffp1.compute(stiffp1);

  if(solver_stiffp1.info()!=Eigen::Success) {
    std::cout << "Failure decomposing stiff plus 1 matrix\n";
  }

  return;

}

void linear::fill_lambda(){ 

#ifdef LOAD
  std::cout << " Reading labmda matrices" << std::endl;
  loadMarket(lambda_x, "lambda_x.mtx");
  loadMarket(lambda_y, "lambda_y.mtx");
#else

  //  div( vectorf,  inter_scalarf );

  std::cout << " Filling lambda matrices" << std::endl;

  std::vector<triplet> ax;            // list of non-zeros coefficients
  std::vector<triplet> ay;            // list of non-zeros coefficients

  int N=0;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    ++N;

    int vi=fv->idx.val();

    typedef Vertex::vector_link vector_link;

    vector_link nabla=fv->nabla();

    for(
	vector_link::iterator nn= nabla.begin();
	nn!=nabla.end(); ++nn) {

      int vj=nn->first->idx.val();

      Vector_2 nnabla=nn->second;

      ax.push_back( triplet(vi,vj, nnabla.x() ));
      ay.push_back( triplet(vi,vj, nnabla.y() ));

    }

  }

  lambda_x.resize(N,N);
  lambda_y.resize(N,N);

  lambda_x.setFromTriplets(ax.begin(), ax.end());
  lambda_y.setFromTriplets(ay.begin(), ay.end());

#endif

  std::cout << " Filled lambda matrices" << std::endl;
  cout << "matrix size " << lambda_x.rows() << " x " << lambda_x.cols() << endl;

  return;

}


void linear::fill_mass() { 

#ifdef LOAD
  std::cout << " Reading mass matrix" << std::endl;
  loadMarket(mass, "mass.mtx");
#else

  std::cout << " Filling mass matrix" << std::endl;

  std::vector<triplet> aa;            // list of non-zeros coefficients

  int N=0;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    ++N;

    int vi=fv->idx.val();

    typedef Vertex::scalar_link scalar_link;

    scalar_link masss=fv->mass();

    for(
	scalar_link::iterator nn= masss.begin();
	nn!=masss.end(); ++nn) {

      int vj=nn->first->idx.val();
	
      FT mm = nn->second;

      aa.push_back( triplet(vi,vj, mm ));
    }
  }

  mass.resize(N,N);

  mass.setFromTriplets(aa.begin(), aa.end());

  std::cout << " Filled mass matrix" << std::endl;
  cout << "matrix size " << mass.rows() << " x " << mass.cols() << endl;

#endif

  // Only non-direct iterative solvers
#ifndef DIRECT_SOLVER
    solver_mass.setTolerance( TOL );
#endif

    solver_mass.compute(mass);

  if(solver_mass.info()!=Eigen::Success) {
    std::cout << "Failure decomposing mass matrix\n";
  }


  return;

}

void linear::fill_mas( const FT& dt ){ 

#ifdef LOAD
  std::cout << " Reading mas matrix" << std::endl;
  loadMarket(mas, "mas.mtx");
  std::cout << " Read mas matrix" << std::endl;
#else

  if(stiff.size()==0) fill_stiff();
  if(mass.size()==0)  fill_mass();

  std::cout << " Building MAS matrix: " ;

  FT a = dt * simu.mu() ;

  cout << " (mass -  " << a << "  x stiff) " << endl;

  mas = mass - a * stiff;

#endif

  solver_mas.compute(mas);

  if(solver_mas.info()!=Eigen::Success) {
    std::cout << "Failure decomposing mass minus a stiff matrix\n";
  }

}



void linear::fill_mbs( const FT& b ){ 


  if(stiff.size()==0) fill_stiff();
  if(mass.size()==0)  fill_mass();

  std::cout << " Building MbS matrix: " ;

  cout << " (mass -  " << b << "  x stiff) " << endl;

  // Model B dynamics:

  mbs = mass + b * stiff;

  solver_mbs.compute(mbs);

  if(solver_mbs.info()!=Eigen::Success) {
    std::cout << "Failure decomposing mass minus b stiff matrix\n";
  }

}


VectorXd linear::field_to_vctr(const kind::f scalarf ) {

  std::vector<int> indices;

  std::vector<FT> ff;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    indices.push_back( fv->idx.val() );
    ff.push_back( fv->sf(scalarf).val() );

  }

  int N=indices.size();
  VectorXd vctr;

  vctr.resize(N);

  for(int nn=0; nn<indices.size(); nn++)  {
    vctr(indices[nn])=ff[nn];

    // cout << indices[nn] << " : "
    //  	 << ff[nn] << endl;
  }

  return vctr;

}

void linear::vctr_to_field(const VectorXd& vv, const kind::f scalarf ) {

  //  int nn=0;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    // check (only makes sense if vertices are re-labeled at every step)
    int index = fv->idx.val();
      
    // if(index!=nn) {
    //   std::cout << "Error transfering from vector onto field " << scalarf << std::endl;
    // }

    fv->sf(scalarf).set(vv[index]);
    //    ++nn;

    // cout << index << " : "
    // 	 << vv[index] << endl;

  }
  return;
}




void linear::vctr_to_vfield(const VectorXd& vv, const kind::f vectorf, const int comp ) {

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    int index = fv->idx.val();

// TODO: define changes on each component of a vector field
    if(comp==0) {
      FT vvy = fv->vf(vectorf).val().y() ;

      fv->vf(vectorf).set(  Vector_2 ( vv[index] , vvy ) );
    } else {
      FT vvx = fv->vf(vectorf).val().x() ;

      fv->vf(vectorf).set(  Vector_2 ( vvx , vv[index] ) );
    }
  }
  return;

}



VectorXd linear::vfield_to_vctr(const kind::f vectorf , int comp ) {

  std::vector<int> indices;

  std::vector<FT> ff;

  for(F_v_it fv=T.finite_vertices_begin();
      fv!=T.finite_vertices_end();
      fv++)  {

    indices.push_back( fv->idx.val() );

    if(comp==0)
      ff.push_back( fv->vf(vectorf).val().x() );
    else
      ff.push_back( fv->vf(vectorf).val().y() );
  }

  int N=indices.size();
  VectorXd vctr(N);

  for(int nn=0; nn<indices.size(); nn++) 
    vctr(indices[nn])=ff[nn];

  return vctr;

}




// solves for scalarf in
// stiffp1 x scalarf = lambda · vectorf
// inter_scalarf is lambda · vectorf

void linear::laplace_div(
			 const kind::f velocity ,
			 const FT dt,
			 const kind::f divvel,
			 const kind::f pressure ) {
  
  if(stiffp1.size()==0)  fill_stiff();

  cout << "Solving Lapl of field " << pressure 
       << " = div of field " << velocity
       << ", with intermediate field " << divvel
       << '\n';

  if(lambda_x.size()==0)  fill_lambda();

  VectorXd vx = vfield_to_vctr( velocity ,0 );
  
  VectorXd vy = vfield_to_vctr( velocity ,1 );

  VectorXd div1 = lambda_x * vx + lambda_y * vy ;

  vctr_to_field( div1 , divvel );

  //  divv.conservativeResize(N+1);

//  divv[N]=0;


  // {
  //   MatrixXd dstiff = MatrixXd(stiff);
  //   cout << "s =" << dstiff.format(OctaveFmt) ;
  // }

  cout << "matrix size " << stiffp1.rows() << " x " << stiffp1.cols() << endl;

  // cout << "d = " << divv.format(OctaveFmt) ;

  //  x = solver.solveWithGuess(b,x0);

  //  if(mass.size()==0)   fill_mass();

  //  VectorXd div2 = mass * div1;

  int N=div1.size();

  VectorXd divv = -(div1.tail( N-1 ));

  cout << "rhs vector of size " << divv.rows() << endl;

  VectorXd lapl= solver_stiffp1.solve( divv );

  // Only non-direct iterative solvers
  //  cout << "Solver iterations: "  << solver_stiffp1.iterations()
  //     << " , max: "  << solver_stiffp1.maxIterations()
  //     << endl ;

  if(solver_stiffp1.info()!=Eigen::Success) 
      cout << "Warning: unsucessful stiffp1 solve, error code " 
	   << solver_stiffp1.info() << endl ;

  // cout << "l = " << lapl.format(OctaveFmt) ;

//  vctr_to_field(lapl, scalarf);
//  mass_s(scalarf); //it's the other way around
  VectorXd lapl2(N);

  lapl2(0) = 0;
  lapl2.tail(N-1) = ( 1.0 / dt ) * lapl;

  vctr_to_field( lapl2 , pressure );

  return;
}


//#undef LUMPED
//#define LUMPED

#ifdef LUMPED

// lumped mass "inversion"

void linear::mass_v(const kind::f vectorf ) {

  for(F_v_it vit=T.vertices_begin();
      vit != T.vertices_end();
      vit++) {
    FT V=vit->vol();

    vit->vf(vectorf) /= V;

  }

  // VectorXd mm = field_to_vctr( kind::VOL  );

  // VectorXd lambda_v_x = vfield_to_vctr( vectorf , 0 );

  // lambda_v_x.cwiseQuotient( mm ) ;

  // VectorXd lambda_v_y = vfield_to_vctr( vectorf , 1 );

  // lambda_v_y.cwiseQuotient( mm ) ;

  // vctr_to_vfield( lambda_v_x  , vectorf , 0 );
  // vctr_to_vfield( lambda_v_y  , vectorf , 1 );


  return;
}

// solves mass x vectorf =  vectorf
// vectorf is I/O
// full mass matrix inversion

#else


// void linear::mass_v( VectorXd& fx , VectorXd& fy ) {

//   cout << "Mass - inverting vector x,y couple " << endl;

//   if(mass.size()==0)   fill_mass();

//   cout << "matrix  vector of size " << mass.rows() << " x " << mass.cols() << endl;

//   cout << "x component" << endl;

//   // x component.-

//   VectorXd grad_x= solver_mass.solve(fx);

//   if(solver_mass.info()!=Eigen::Success) 
//       cout << "Warning: unsucessful mass x solve, error code " 
// 	   << solver_mass.info() << endl ;

//   // y component.-
//   cout << "y component" << endl;

//   VectorXd grad_y= solver_mass.solve(fy);

//   if(solver_mass.info()!=Eigen::Success) 
//     cout << "Warning: unsucessful mass y solve, error code " 
// 	 << solver_mass.info() << endl ;

//   fx=grad_x;
//   fy=grad_y;

//   cout << "Mass invertion done " << endl;

//   return;
// }




void linear::mass_v(const kind::f vectorf ) {


  cout << "Mass - inverting vector field " << vectorf 
       << '\n';


  VectorXd lambda_v_x = vfield_to_vctr( vectorf , 0 );
  VectorXd lambda_v_y = vfield_to_vctr( vectorf , 1 );

  mass_s(lambda_v_x);
  mass_s(lambda_v_y);

  vctr_to_vfield( lambda_v_x , vectorf , 0 );
  vctr_to_vfield( lambda_v_y , vectorf , 1 );

  cout << "Mass invertion of vector field " << vectorf << "  done " << endl;

  return;
}

void linear::mass_s( VectorXd& lambda ) {

  cout << "Mass - inverting scalar field " << endl;

  if(mass.size()==0)   fill_mass();

  VectorXd grad= solver_mass.solve(lambda);

  if(solver_mass.info()!=Eigen::Success) 
      cout << "Warning: unsucessful mass x solve, error code " 
	   << solver_mass.info() << endl ;

  lambda = grad;

  return;
}



void linear::mass_s(const kind::f scalarf ) {

  cout << "Mass - inverting scalar field " << scalarf 
       << '\n';

  VectorXd lambda = field_to_vctr( scalarf );

  mass_s( lambda );
  
  vctr_to_field( lambda , scalarf );

  return;
}

#endif


// solves (mass - dt x mu x stiff)  ustar = mass ( dt x force + u0 )

void linear::ustar_inv_cp(
			  const kind::f Ustar,
			  const FT dt,
			  const kind::f U0 ,
			  const bool overdamped, const bool semi) {

  if(mas.size()==0)  fill_mas( dt );
  if(lambda_x.size()==0)  fill_lambda();

  //  FT eps=1e-1;
  
//// x
  VectorXd al = field_to_vctr( kind::ALPHA );
    
  VectorXd U0_x = vfield_to_vctr( kind::FORCE , 0 ); // = dt * 

  // gradient of chem pot .-

//  cout << " chem pot" << endl;

  VectorXd chemp = field_to_vctr( kind::CHEMPOT );

//  cout << "grad chem pot" << endl;

  VectorXd grad_cp_x = lambda_x * chemp;
  VectorXd grad_cp_y = lambda_y * chemp;

//  cout << "Inverting grad chem pot" << endl;

  mass_s( grad_cp_x );
  mass_s( grad_cp_y );

  cout << "Grad chem pot done, for u*" << endl;

  // optional .-
//  cout << "grad  chem pot back" << endl;

  vctr_to_vfield( grad_cp_x  , kind::GRADCHEMPOT , 0 );
  vctr_to_vfield( grad_cp_y  , kind::GRADCHEMPOT , 1 );

  //


  U0_x -=  al * grad_cp_x ;

  U0_x *= dt;
  
  if(!overdamped) {
    //    VectorXd force_x
    U0_x += vfield_to_vctr( U0 , 0 );
  }
  
  VectorXd massU0_x;

  if (semi) {

    VectorXd grad_x = vfield_to_vctr( kind::GRADP , 0 );

    massU0_x = mass * U0_x - dt * grad_x;

  } else
    massU0_x =  mass * U0_x ;

  VectorXd Ustarx= solver_mas.solve(massU0_x);
  if(solver_mas.info()!=Eigen::Success) 
      cout << "Warning: unsucessful mas x solve, error code " 
	   << solver_mas.info() << endl ;


  vctr_to_vfield( Ustarx , Ustar , 0 );

//// y

  VectorXd U0_y = vfield_to_vctr( kind::FORCE , 1 ); // = dt * 

  U0_y -= al *  grad_cp_y ;
  //  U0_y -= eps * al * vfield_to_vctr( kind::GRADALPHA , 1 );

  U0_y *= dt;

  if(!overdamped) {
    //  VectorXd force_y = vfield_to_vctr( kind::FORCE , 1 );

    U0_y += vfield_to_vctr( U0 , 1 );
  }

  
  VectorXd massU0_y;

  if (semi) {
    VectorXd grad_y = vfield_to_vctr( kind::GRADP , 1 );
    massU0_y =  mass * U0_y - dt * grad_y;

  } else
    massU0_y =  mass * U0_y ;

  VectorXd Ustary= solver_mas.solve(massU0_y);

  if(solver_mas.info()!=Eigen::Success) 
      cout << "Warning: unsucessful mas y solve, error code " 
	   << solver_mas.info() << endl ;


  vctr_to_vfield( Ustary , Ustar  , 1);

  return;
}



// solves (mass - dt x mu x stiff)  ustar = mass ( dt x force + u0 )

void linear::ustar_inv(
			  const kind::f Ustar,
			  const FT dt,
			  const kind::f U0 ,
			  const bool overdamped, const bool semi) {

  if(mas.size()==0)  fill_mas( dt );
    
  VectorXd U0_x = vfield_to_vctr( kind::FORCE , 0 ); // = dt * 

  U0_x *= dt;
  
  if(!overdamped) {
    //    VectorXd force_x
    U0_x += vfield_to_vctr( U0 , 0 );
  }
  
  VectorXd massU0_x;

  if (semi) {

    VectorXd grad_x = vfield_to_vctr( kind::GRADP , 0 );

    massU0_x = mass * U0_x - dt * grad_x;

  } else
    massU0_x =  mass * U0_x ;

  VectorXd Ustarx= solver_mas.solve(massU0_x);
  if(solver_mas.info()!=Eigen::Success) 
      cout << "Warning: unsucessful mas x solve, error code " 
	   << solver_mas.info() << endl ;


  vctr_to_vfield( Ustarx , Ustar , 0 );

//// y

  VectorXd U0_y = vfield_to_vctr( kind::FORCE , 1 ); // = dt * 

  U0_y *= dt;

  if(!overdamped) {
    //  VectorXd force_y = vfield_to_vctr( kind::FORCE , 1 );

    U0_y += vfield_to_vctr( U0 , 1 );
  }

  
  VectorXd massU0_y;

  if (semi) {
    VectorXd grad_y = vfield_to_vctr( kind::GRADP , 1 );
    massU0_y =  mass * U0_y - dt * grad_y;

  } else
    massU0_y =  mass * U0_y ;

  VectorXd Ustary= solver_mas.solve(massU0_y);

  if(solver_mas.info()!=Eigen::Success) 
      cout << "Warning: unsucessful mas y solve, error code " 
	   << solver_mas.info() << endl ;


  vctr_to_vfield( Ustary , Ustar  , 1);

  return;
}



// solves (mass - dt x mu x stiff)  alpha = alpha0

void linear::alpha_inv_cp(const kind::f alpha,
		       const FT dt,
		       const kind::f alpha0 ) {
  FT D=0.04;

  FT b = D*dt ;

  VectorXd al0 = field_to_vctr( alpha );
  VectorXd al  = field_to_vctr( alpha );

  //VectorXd al3 = al.array().pow(3);
  //VectorXd mass0 = mass * al0 + b * stiff * al3;
  // VectorXd al = solver_mbs.solve(mass0);


  // Model B dynamics:

  if(mbs.size()==0)  fill_mbs( b );
  
  VectorXd cp2 = field_to_vctr( kind::CHEMPOT ) + al ;

  VectorXd mass0 = mass * al0 + b * stiff * cp2;

  VectorXd new_al = solver_mbs.solve(mass0);

  if(solver_mbs.info()!=Eigen::Success) 
    cout << "Warning: unsucessful mbs x solve, error code " 
	 << solver_mbs.info() << endl ;
  ////
  
  // Model A dynamics:
  //  VectorXd cp2 = - al.array() * al.array() * al.array() ;
  //  VectorXd new_al = 1.0/(1.0-b)*( al0 + b*cp2 );
  ////
  
  vctr_to_field( new_al , alpha );

  return;
}


//#define compile_uhalf_inv
#ifdef compile_uhalf_inv

// solves   dt * mu * stiff vectorf =  mass ( ustar - U0 - dt * force )

void linear::uhalf_inv(const kind::f U, const FT dt, const kind::f U0 ) {


  if(stiff.size()==0) fill_stiff();
  if(mass.size()==0)  fill_mass();

//// x

  VectorXd Ustar_x = vfield_to_vctr( kind::USTAR , 0 );
  VectorXd U0_x    = vfield_to_vctr( U0    , 0 );
  VectorXd force_x = vfield_to_vctr( kind::FORCE , 0 );

  VectorXd mass_x1;

  FT inv_mu = 1.0/simu.mu();
  mass_x1 =  mass * ( inv_mu * (
				(1.0/dt) * ( Ustar_x - U0_x )
				- force_x
				)
		      );

  int N=mass_x1.size();

  VectorXd mass_x = mass_x1.tail( N-1 );

  VectorXd Ux= solver_stiffp1.solve(mass_x);

  if(solver_stiffp1.info()!=Eigen::Success) 
      cout << "Warning: unsucessful mas x solve, error code " 
	   << solver_stiffp1.info() << endl ;

  VectorXd Ux1(N);

  Ux1(0) = 0;
  Ux1.tail(N-1) = Ux;

  vctr_to_vfield( Ux1 , U , 0 );

//// y


  VectorXd Ustar_y = vfield_to_vctr( kind::USTAR , 1 );
  VectorXd U0_y    = vfield_to_vctr( U0    , 1 );
  VectorXd force_y = vfield_to_vctr( kind::FORCE , 1 );

  VectorXd mass_y1;

  mass_y1 =  mass * ( inv_mu * (
				(1.0/dt) * ( Ustar_y - U0_y )
				- force_y
				)
		      );

  VectorXd mass_y = mass_y1.tail( N-1 );

  VectorXd Uy= solver_stiffp1.solve(mass_y);

  if(solver_stiffp1.info()!=Eigen::Success) 
      cout << "Warning: unsucessful mas x solve, error code " 
	   << solver_stiffp1.info() << endl ;

  VectorXd Uy1(N);

  Uy1(0) = 0;
  Uy1.tail(N-1) = Uy;

  vctr_to_vfield( Uy1 , U , 1 );

  return;
}

#endif



