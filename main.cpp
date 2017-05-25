// polyFEM
// with mass matrix, and all
// plus, correction to quadratic consistency
// periodic boundary conditions

// Cahn-Hilliard solver, overdamped regime, reduced units
// check "D" in linear::chempot_inv for (I think) the only
// reduced parameter

// projected to and from mesh

#include <CGAL/Timer.h>

// write out matrices
//#define WRITE

//#define EXPLICIT

#include"main.h"


#include"sim_pars.h"

#include"linear.h"

#include"fields.h"

// Init global stuff.-

#include"periodic.h"

const FT LL=10; // length of original domain

Iso_rectangle domain(-LL/2, -LL/2, LL/2, LL/2);


// TODO: the two triangulations store different things.
//       specific bases and faces should be implemented for each

sim_pars simu;


//#define FULL
#define FULL_FULL
//#define FULL_LUMPED
//#define FLIP

#ifdef FULL_FULL
#define FULL
#endif

#ifdef FULL_LUMPED
#define FULL
#endif

#include"onto_from_mesh.h"



//const Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "];");


Triangulation Tp(domain); // particles
Triangulation Tm(domain); // mesh


int main() {

//  CGAL::Timer time;
//
//  time.start();
  
  cout << "Creating point cloud" << endl;

  simu.read();

  create();
  
  if(simu.create_points()) {

    //    set_alpha_circle( Tp , 2);
    //    set_alpha_under_cos(  Tp ) ;

    cout << "Creating alpha field " << endl;
    
    set_alpha_random(  Tm ) ;

    cout << "Numbering particles " << endl;

    number(Tp);
    number(Tm);
  }
  

  // every step
  areas(Tp);
  quad_coeffs(Tp , simu.FEMp() ); volumes(Tp, simu.FEMp() );

  // just once!
  linear algebra(Tm);

  areas(Tm);
  quad_coeffs(Tm , simu.FEMm() ); volumes(Tm, simu.FEMm() );

  cout << "Setting up diff ops " << endl;

  if(simu.create_points()) {
    nabla(Tm);
    Delta(Tm);
  }

  const std::string mesh_file("mesh.dat");
  const std::string particle_file("particles.dat");

  // step 0 draw.-
    draw(Tm, mesh_file     , true);
    draw(Tp, particle_file , true);
  
  cout << "Assigning alpha to particles " << endl;

  
// #if defined FULL_FULL
//   {
//     Delta(Tp);
//     linear algebra_p(Tp);
//     from_mesh_full( Tm , Tp ,  algebra_p,kind::ALPHA);
//   }
// #elif defined FULL_LUMPED
//   from_mesh_lumped( Tm , Tp , kind::ALPHA);
//  #elif defined FLIP
//   from_mesh(Tm , Tp , kind::ALPHA);
//  #else
//   from_mesh(Tm , Tp , kind::ALPHA);
// #endif

  cout << "Moving info" << endl;
  move_info( Tm );
  move_info( Tp );

  draw(Tm, mesh_file     , true);
  draw(Tp, particle_file , true);

  simu.advance_time();
  simu.next_step();

  //  bool first_iter=true;

  CGAL::Timer time;

  time.start();

  std::ofstream log_file;

  log_file.open("main.log");

  bool is_overdamped = ( simu.mu() > 1 ) ; // high or low Re

  for(;
      simu.current_step() <= simu.Nsteps();
      simu.next_step()) {

    cout
      << "Step " << simu.current_step() 
      << " . Time " <<  simu.time()
      << " ; t step " << simu.dt()
      << endl;

    FT dt=simu.dt();

    FT dt2 = dt / 2.0 ;

    int iter=0;
    FT displ=1e10;

    FT min_displ=1e10;
    int min_iter=0;

    const int max_iter=1; //10;
    const FT  max_displ=  1e-8; // < 0 : disable

//  leapfrog, special first step.-
//    if(simu.current_step() == 1) dt2 *= 0.5;

//    dt2 *= 0.5;

    move_info(Tm);
    move_info(Tp);

    // iter loop
    for( ; iter<max_iter ; iter++) {


//       cout << "Projecting U , alpha0 from mesh " << endl;

// #if defined FULL_FULL
//       {
// 	Delta(Tp);
// 	linear algebra_p(Tp);
// 	from_mesh_full_v(Tm, Tp, algebra_p , kind::U);
// 	from_mesh_full  (Tm, Tp, algebra_p , kind::ALPHA0);
// 	from_mesh_full  (Tm, Tp, algebra_p , kind::ALPHA);
//       }
// #elif defined FULL_LUMPED
//       from_mesh_lumped_v(Tm, Tp, kind::U);
//       from_mesh_lumped  (Tm, Tp, kind::ALPHA0);
//       from_mesh_lumped  (Tm, Tp, kind::ALPHA);
// #elif defined FLIP
//       from_mesh_v(Tm, Tp, kind::U);
//       from_mesh  (Tm, Tp, kind::ALPHA0);
//       from_mesh  (Tm, Tp, kind::ALPHA);
// #else
//       from_mesh_v(Tm, Tp, kind::U);
//       from_mesh  (Tm, Tp, kind::ALPHA0);
//       from_mesh  (Tm, Tp, kind::ALPHA);
// #endif

      // comment for no move.-
      //     displ = move( Tp , dt2 );

      cout << "Iter " << iter << " , moved avg " << displ << " to half point" << endl;

      if( displ < min_displ) {
	min_displ=displ;
	min_iter=iter;
      }

      if( (displ < max_displ) && (iter !=0) ) break;

//       areas(Tp);
//       quad_coeffs(Tp , simu.FEMp() ); volumes(Tp, simu.FEMp() );
      
//       cout << "Proj U0, alpha0 onto mesh " << endl;

// #if defined FULL
//       onto_mesh_full_v(Tp,Tm,algebra,kind::UOLD);
//       onto_mesh_full  (Tp,Tm,algebra,kind::ALPHA0);
//       onto_mesh_full  (Tp,Tm,algebra,kind::ALPHA);
// #elif defined FLIP
//       flip_volumes(Tp , Tm , simu.FEMm() );
//       onto_mesh_flip_v(Tp,Tm,simu.FEMm(),kind::UOLD);
//       onto_mesh_flip  (Tp,Tm,simu.FEMm(),kind::ALPHA0);
//       onto_mesh_flip  (Tp,Tm,simu.FEMm(),kind::ALPHA);
// #else
//       onto_mesh_delta_v(Tp,Tm,kind::UOLD);
//       onto_mesh_delta  (Tp,Tm,kind::ALPHA0);
//       onto_mesh_delta  (Tp,Tm,kind::ALPHA);
// #endif

      // partly explicit ( unstable ? ):
      cout << "Calculating chem pot explicitely" << endl;

	// inner iter loop

      for( int alpha_it=0 ; alpha_it < 1 ; alpha_it++) { // max_iter ; alpha_it++) {

	cout << "Alpha loop iter " << alpha_it << endl;

	if (iter==0)
	  algebra.chempot( kind::ALPHA0, kind::CHEMPOT );
	else
	  algebra.chempot( kind::ALPHA , kind::CHEMPOT );

	  //	  algebra.chempot( kind::ALPHA , kind::CHEMPOT );

	//	algebra.alpha_inv_cp(kind::ALPHA, dt2 , kind::ALPHA0 );

	algebra.alpha_explicit(kind::ALPHA, dt2 , kind::ALPHA0 );
      }

      // cout << "Getting chem pot force" << endl;

      // algebra.chem_pot_force();

      // // substract spurious overall movement.-
      
      // zero_mean_v( Tm , kind::FORCE);

      // cout << "Solving PPE" << endl;

      // comment for no move.-
      //      algebra.PPE( kind::FORCE , 1 , kind:: P ); // Dt set to 1

      //      algebra.u_inv_od(kind::U);

      //      zero_mean_v( Tm , kind::U);


      //      cout << "Evolving U " << endl;

      // comment for no move.-
      // // NOOO  u_new( Tm , dt2 );
      //      cout << "U evolved " << endl;

    } // iter loop

// #if defined FULL_FULL
//       {
// 	Delta(Tp);
// 	linear algebra_p(Tp);
// 	from_mesh_full_v(Tm, Tp, algebra_p , kind::U);
// 	from_mesh_full  (Tm, Tp, algebra_p , kind::ALPHA);
//       }
// #elif defined FULL_LUMPED
//       from_mesh_lumped_v(Tm, Tp, kind::U);
//       from_mesh_lumped  (Tm, Tp, kind::ALPHA);
// #elif defined FLIP
//       from_mesh_v(Tm, Tp, kind::U);
//       from_mesh  (Tm, Tp, kind::ALPHA);
// #else
//       from_mesh_v(Tm, Tp, kind::U);
//       from_mesh  (Tm, Tp, kind::ALPHA);
// #endif

      // comment for no move.-
      //displ=move( Tp , dt );

//    update_half_velocity( Tp , false ); 

    // comment for no move.-
      //    update_half_velocity( Tp , is_overdamped ); 

    update_half_alpha( Tm );

//     areas(Tp);

//     quad_coeffs(Tp , simu.FEMp() ); volumes(Tp, simu.FEMp() );

//     cout << "Proj U_t+1 , alpha_t+1 onto mesh " << endl;

// #if defined FULL
//     onto_mesh_full_v(Tp,Tm,algebra,kind::U);
//     onto_mesh_full  (Tp,Tm,algebra,kind::ALPHA0);
//     onto_mesh_full  (Tp,Tm,algebra,kind::ALPHA);
// #elif defined FLIP
//     flip_volumes(Tp , Tm , simu.FEMm() );
//     onto_mesh_flip_v(Tp,Tm,simu.FEMm(),kind::U);
//     onto_mesh_flip  (Tp,Tm,simu.FEMm(),kind::ALPHA0);
//     onto_mesh_flip  (Tp,Tm,simu.FEMm(),kind::ALPHA);
// #else
//     onto_mesh_delta_v(Tp,Tm,kind::U);
//     onto_mesh_delta  (Tp,Tm,kind::ALPHA);
//     onto_mesh_delta  (Tp,Tm,kind::ALPHA);
// #endif


    
    if(simu.current_step()%simu.every()==0) {
      draw(Tm, mesh_file     , true);
      draw(Tp, particle_file , true);
    }

    log_file
      << simu.current_step() << "  "
      <<  simu.time() << "  " ;

    // integrals( Tp , log_file);     log_file << "  ";
    // fidelity(  Tp , log_file );    log_file << endl;

    simu.advance_time();

  } // time loop

  time.stop();

  log_file.close();
  
  cout << "Total runtime: " << time.time() << endl;
  return 0;

}



void create(void) {

  int N=simu.no_of_particles();
  std::vector<Point> points;
  //  points.reserve(N);

  if(simu.create_points()) {
    if(simu.at_random()) {
      points.reserve(N);

      CGAL::Random_points_in_square_2<Point,Creator> g(LL/2.0-0.0001);
      CGAL::copy_n( g, N, std::back_inserter(points));

      cout << N << "  particles placed at random" << endl;
    } else {
      // if((plotting)&&(Nin%2==0)&&(spike)) {
      // 	cout << "Please enter an odd number of particles" << endl;
      // 	std::abort();
      // }

      int Nb=sqrt(N + 1e-12);

      N=Nb*Nb;

      simu.set_no_of_particles(N);

      points.reserve(N);
      cout << N << " particles placed on square lattice" << endl;

      FT spacing=LL/FT(Nb+0);
      FT side=LL-1*spacing;

      points_on_square_grid_2(side/2.0, N, std::back_inserter(points),Creator());;

      if(simu.perturb()) {
	CGAL::perturb_points_2(
			       points.begin(), points.end(),
			       simu.pert_rel()* spacing );//,Creator());
	cout << "each particle perturbed about " << simu.pert_rel()* spacing  << endl;
      }

    }

    cout << "Inserting" << endl;

    Tp.insert(points.begin(), points.end());

    
    points.clear();

    // int Nb = sqrt(N + 1e-12);
    // int nm = Nb* simu.mesh_factor() + 1 ;
    // int Nm = nm * nm;

    int Nm=simu.no_of_nodes();

    int nm=sqrt(Nm + 1e-12);

    Nm= nm * nm;

    simu.set_no_of_nodes(Nm);

    points.reserve(Nm);
    cout << Nm << " mesh on square lattice" << endl;

    FT spacing=LL/FT( nm +0);
    FT side=LL-1*spacing;

    points_on_square_grid_2(side/2.0, Nm , std::back_inserter(points),Creator());;

    // TODO: perfectly regular square grids are not too good, in fact
    CGAL::perturb_points_2(
    			   points.begin(), points.end(),
    			   0.001* spacing );

    Tm.insert(points.begin(), points.end());


    

  } else {

    int N=simu.no_of_particles();

    char part_file[]="particles.dat";

    cout << "reading from file : " << part_file << endl;

    std::ifstream main_data;
    main_data.open(part_file );

    if(main_data.fail()){
      cout << "part file not found "  << endl;
      abort();
    }
    
    for(int i=0;i<N;i++) {
      FT x,y;
      main_data >> x;
      main_data >> y;

      //      cout << x << "  " << y << endl;

      Vertex_handle vh=Tp.insert(Point(x,y));

#include"readin.h"


    }
  
    cout << "particles' data read" << endl;

    main_data.close();

  }

  // straight from the manual.-

  Triangulation::Covering_sheets cs = Tp.number_of_sheets();

  cout << "Original covering (particles): " << cs[0] << ' ' << cs[1] << endl;

//  return ;
  
  Tp.convert_to_1_sheeted_covering();

  cs = Tp.number_of_sheets();

  cout << "Current covering (particles): " << cs[0] << ' ' << cs[1] << endl;


  return ;

  if ( Tp.is_triangulation_in_1_sheet() ) // = true
    {
      bool is_extensible = Tp.is_extensible_triangulation_in_1_sheet_h1()
	|| Tp.is_extensible_triangulation_in_1_sheet_h2(); // = false
      Tp.convert_to_1_sheeted_covering();
      cs = Tp.number_of_sheets();
      cout << "Current covering: " << cs[0] << ' ' << cs[1] << endl;
      if ( is_extensible ) // = false
	cout << "It is safe to change the triangulation here." << endl;
      else {
	cout << "It is NOT safe to change the triangulation here!" << endl;
	abort();
      }
      //      T.convert_to_9_sheeted_covering();
      //      cs = T.number_of_sheets();
      //      cout << "Current covering: " << cs[0] << ' ' << cs[1] << endl;
    } else {
	cout << "Triangulation not on one sheet!" << endl;
	abort();    
  }


  //  cout << "It is (again) safe to modify the triangulation." << endl;


  return ;

}



void number(Triangulation& T) {

  int i=0;

  for(F_v_it vit=T.vertices_begin();
      vit != T.vertices_end();
      vit++) {
    //    vit->indx.set(i); //or
    vit->idx=i;
    ++i;
  }

  return;
}



  // for (
  //      Periodic_point_iterator pit=
  // 	 T.periodic_points_begin(stored_cover);
  //      pit != T.periodic_points_end(stored_cover);
  //      ++pit)
  //   {
  //     //  pt = *ptit;
  //     //      if (! (pt[0].second.is_null() && pt[1].second.is_null() && pt[2].second.is_null()) )
  //     //        {
  //         // Convert the current Periodic_triangle to a Triangle if it is
  //         // not strictly contained inside the original domain.
  //         // Note that this requires EXACT constructions to be exact!
  //     //          t_bd = T.triangle(pt);
  // 	  //}

  // 	  Point p=pit->first;
  // 	  Offset os=pit->second;

  // 	  interior
  // 	    << p.x()+os[0]*LL << "  " 
  // 	    << p.y()+os[1]*LL << "  " 
  // 	    //	    << vit->indx() 
  // 	    << endl;

  //   }


  // for(F_v_it vit=T.vertices_begin();
  //     vit != T.vertices_end();
  //     vit++)
  //   interior
  //     << vit->point().x() << "  " 
  //     << vit->point().y() << "  " 
  //     << vit->indx() << endl;

  // return;
