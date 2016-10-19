// with mass matrix, and all
// plus, correction to quadratic consistency
// periodic boundary conditions

#include <CGAL/Timer.h>

// write out matrices
//#define WRITE

//  #define EXPLICIT

#include"main.h"
#include"sim_pars.h"
#include"linear.h"


// Init global stuff.-

#include"periodic.h"

const FT LL=3; // length of original domain

Iso_rectangle domain(-LL/2, -LL/2, LL/2, LL/2);

#include"fields.h"

// TODO: the two triangulations store different things.
//       specific bases and faces should be implemented for each
Triangulation Tp(domain); // particles
Triangulation Tm(domain); // mesh

sim_pars simu;

//const Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "];");

// stuff only used here:

//#define FULL
//#define FULL_FULL
//#define FULL_LUMPED
//#define FLIP

#ifdef FULL_FULL
#define FULL
#endif

#ifdef FULL_LUMPED
#define FULL
#endif


#include"onto_from_mesh.h"

int main() {

//  CGAL::Timer time;
//
//  time.start();


  cout << "Creating point cloud" << endl;

  simu.read();

  create();

  if(simu.create_points()) {
    set_fields_Zalesak();
    number();
  }

  // // every step
  // areas(Tp);  quad_coeffs(Tp , simu.FEMp() ); volumes(Tp, simu.FEMp() );   Delta(Tp);

  // just once!

  areas(Tm);  quad_coeffs(Tm , simu.FEMm() );
  volumes(Tm, simu.FEMm() );
  Delta(Tm);

  linear algebra(Tm);

  // if(simu.create_points()) {
  //   nabla();
  //   Delta();
  // }

  
  // just for the looks of step 0:
   // onto_mesh_lumped();
// #ifdef FULL
//   onto_mesh_full(algebra);
// #else
//   onto_mesh_delta();
// #endif

  move_info();

  // /// Prev test begin
  //cout << "Calculating Lapl U" << endl;
  //algebra.laplacian_v(kind::UOLD,kind::LAPLU);

  //FT dt=simu.dt();

  //cout << "Calculating Ustar implicitely" << endl;
  //algebra.ustar_inv(kind::USTAR,  dt , kind::UOLD, false);

  //cout << "Solving PPE" << endl;
  //algebra.PPE( kind::USTAR, dt, kind:: P );

  //cout << "Calculating grad p" << endl;
  //algebra.gradient(kind::P, kind::GRADP);
  //algebra.mass_s(kind::DIVU);

  
//draw();
//  return 1;

  // /// Prev test end

#ifdef WRITE
  algebra.save_matrices();
#endif

  //  set_fields();

  //  set_vels();

  draw();  

  simu.advance_time();
  simu.next_step();

  //  bool first_iter=true;

  CGAL::Timer time;

  time.start();

  std::ofstream log_file;

  log_file.open("main.log");

  areas(Tp);  quad_coeffs(Tp , simu.FEMp() ); volumes(Tp, simu.FEMp() );
     
  for(;
      simu.current_step() <= simu.Nsteps();
      simu.next_step()) {

    cout
      << "Step " << simu.current_step() 
      << " . Time " <<  simu.time()
      << " ; t step " << simu.dt()
      << endl;

    FT displ;

    FT dt=simu.dt();

    FT dt2 = dt / 2.0 ;

//  leapfrog, special first step.-
//    if(simu.current_step() == 1) dt2 *= 0.5;

//    dt2 *= 0.5;

    cout << "Transfering info to particles" << endl; 

    areas(Tp);  quad_coeffs(Tp , simu.FEMp() ); volumes(Tp, simu.FEMp() );

#if defined FULL_FULL
    {
      Delta(Tp);
      linear algebra_p(Tp);
      from_mesh_full(algebra_p,kind::ALPHA);
    }
#elif defined FULL_LUMPED
    from_mesh_lumped(kind::ALPHA);
#elif defined FLIP
    from_mesh(kind::ALPHA);
#else
    from_mesh(kind::ALPHA);
#endif

    //    from_mesh(kind::ALPHA);

// leapfrog.-
//      displ=move( 0.5 * dt );

// half step.-
//    set_vels();

    displ=move( dt2 );

    cout << "Moved avg " << displ << " to half point" << endl;

    set_vels_rotating();

    displ=move( dt );

    cout << "Moved avg " << displ << " from half point" << endl;

    move_info();

    cout << "Transfering info onto mesh " << endl;

    areas(Tp);  quad_coeffs(Tp , simu.FEMp() ); volumes(Tp, simu.FEMp() );
    
      //onto_mesh_lumped();
#if defined FULL
    onto_mesh_full(algebra,kind::ALPHA);
#elif defined FLIP
    flip_volumes(Tp , Tm , simu.FEMm() );
    onto_mesh_flip(Tp,Tm,simu.FEMm(),kind::ALPHA);
#else
    onto_mesh_delta(kind::ALPHA);
#endif

    move_info();

    if(simu.current_step()%simu.every()==0) draw();

    log_file
      << simu.current_step() << "  "
      <<  simu.time() << "  " ;

    integrals( Tm , log_file);     log_file << "  ";
    fidelity(Tm,log_file );        log_file << "  ";
    integrals( Tp , log_file);     log_file << "  ";
    fidelity(  Tp , log_file );        log_file << endl;

    simu.advance_time();

  }

  time.stop();

  log_file.close();
  
  cout << "Total runtime: " << time.time() << endl;
  return 0;

}



void create(void) {

  int N=simu.no_of_particles();
  std::vector<Point> points;
  points.reserve(N);

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

    points.reserve(Nm);
    cout << Nm << " mesh on square lattice" << endl;

    int nm=sqrt(Nm + 1e-12);

    FT spacing=LL/FT( nm +0);
    FT side=LL-1*spacing;

    points_on_square_grid_2(side/2.0, Nm , std::back_inserter(points),Creator());;

    // TODO: perfectly regular square grids are not too good, in fact
    CGAL::perturb_points_2(
			   points.begin(), points.end(),
			   0.001* spacing );

    Tm.insert(points.begin(), points.end());

//    if(simu.initial_velocity())
//      setup_v();

//// Insert in circle only
//    N=0;
//    for( std::vector<Point>::iterator pp=points.begin();
//         pp<points.end(); pp++) {
//         FT x=pp->x();
//         FT y=pp->y();
//         FT rr=x*x+y*y;
//         if (rr < 0.5*0.5) {
//           T.insert(*pp);
//           N++;
//         }
//    }
//    simu.set_no_of_points(N);

  } else {

    int N=simu.no_of_particles();

    char part_file[]="particles.dat";

    cout << "reading from file : " << part_file << endl;

    std::ifstream main_data;
    main_data.open(part_file );

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

    char mesh_file[]="mesh.dat";

    cout << "reading from file : " << mesh_file << endl;

    main_data.open(mesh_file );

    int Nm=simu.no_of_nodes();

    for(int i=0;i<Nm;i++) {
      FT x,y;
      main_data >> x;
      main_data >> y;

      //      cout << x << "  " << y << endl;

      Vertex_handle vh=Tm.insert(Point(x,y));

#include"readin.h"


    }
  
    cout << "mesh data read" << endl;

    main_data.close();

  }

  // straight from the manual.-

  Triangulation::Covering_sheets cs = Tp.number_of_sheets();

  cout << "Original covering (particles): " << cs[0] << ' ' << cs[1] << endl;

//  return ;
  
  Tp.convert_to_1_sheeted_covering();

  cs = Tp.number_of_sheets();

  cout << "Current covering (particles): " << cs[0] << ' ' << cs[1] << endl;

  cs = Tm.number_of_sheets();

  cout << "Original covering (mesh): " << cs[0] << ' ' << cs[1] << endl;

//  return ;
  
  Tm.convert_to_1_sheeted_covering();

  cs = Tm.number_of_sheets();

  cout << "Current covering (mesh): " << cs[0] << ' ' << cs[1] << endl;

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



void number(void) {

  int i=0;

  for(F_v_it vit=Tp.vertices_begin();
      vit != Tp.vertices_end();
      vit++) {
    //    vit->indx.set(i); //or
    vit->idx=i;
    ++i;
  }

  i=0;

  for(F_v_it vit=Tm.vertices_begin();
      vit != Tm.vertices_end();
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
