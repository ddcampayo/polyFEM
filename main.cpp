// with mass matrix, and all
// plus, correction to quadratic consistency
// periodic boundary conditions

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

const FT LL=1; // length of original domain

Iso_rectangle domain(-LL/2, -LL/2, LL/2, LL/2);


// TODO: the two triangulations store different things.
//       specific bases and faces should be implemented for each

sim_pars simu;

//const Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "];");

Triangulation Tp(domain); // particles

int main() {

//  CGAL::Timer time;
//
//  time.start();
  
  cout << "Creating point cloud" << endl;

  simu.read();

  create();

  if(simu.create_points()) {

    set_alpha_circle( Tp );

    number(Tp);
  }
  
  // areas(Tp);
  // quad_coeffs(Tp , simu.FEMp() ); volumes(Tp, simu.FEMp() );

  // volumes(Tp, simu.FEMp() );
  // Delta(Tp);

  // linear algebra(Tp);

  // if(simu.create_points()) {
  //   nabla(Tp);
  //   Delta(Tp);
  // }
  
  move_info( Tp );

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

   /// Prev test end

#ifdef WRITE
  algebra.save_matrices();
#endif

      // areas(Tp);
      // quad_coeffs(Tp , simu.FEMp() ); volumes(Tp, simu.FEMp() );

      // nabla(Tp);
      // Delta(Tp);


      // linear algebra(Tp);

      // cout << "Calculating grad alpha" << endl;
      // algebra.gradient(kind::ALPHA, kind::GRADALPHA);
  
  const std::string particle_file("particles.dat");

  draw(Tp, particle_file , true);
  
  simu.advance_time();
  simu.next_step();

  //  bool first_iter=true;

  CGAL::Timer time;

  time.start();

  std::ofstream log_file;

  log_file.open("main.log");

  bool overdamped = ( simu.mu() > 1 ) ; // high or low Re

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

    const int max_iter=10;
    const FT  max_displ= 1e-8; // < 0 : disable

//  leapfrog, special first step.-
//    if(simu.current_step() == 1) dt2 *= 0.5;

//    dt2 *= 0.5;

    move_info(Tp);

    // iter loop
    for( ; iter<max_iter ; iter++) {

      cout << "Move iteration  " << iter << " of " << max_iter << " ";

      displ=move( Tp , dt2 );

      cout << "Moved avg " << displ << " to half point" << endl;

      if( (displ < max_displ) && (iter !=0) ) break;

      areas(Tp);
      quad_coeffs(Tp , simu.FEMp() ); volumes(Tp, simu.FEMp() );

      nabla(Tp);
      Delta(Tp);

      linear algebra(Tp);
      
      if( displ < min_displ) {
	min_displ=displ;
	min_iter=iter;
      }

      //      set_forces_Kolmo(Tp);

//  Reynolds number discrimination



#ifdef EXPLICIT

	cout << "Calculating Ustar explicitely" << endl;

	algebra.laplacian_v(kind::UOLD,kind::LAPLU);

	u_star(Tp, dt2 , false );

#else

	//	cout << "Calculating alpha implicitely" << endl;

	//	algebra.alpha_inv(kind::ALPHA,  dt2, kind::ALPHA0 );

	cout << "Calculating grad alpha" << endl;
	algebra.gradient(kind::ALPHA, kind::GRADALPHA);

	cout << "Calculating Ustar implicitely" << endl;

	//	algebra.ustar_inv(kind::USTAR,  dt2 , kind::UOLD, false , false);

	algebra.ustar_inv(kind::USTAR,  dt2 , kind::UOLD, overdamped , false);

#endif

	cout << "Solving PPE" << endl;
      
	algebra.PPE( kind::USTAR, dt2 , kind:: P );

	cout << "Calculating grad p" << endl;
	algebra.gradient(kind::P, kind::GRADP);

	cout << "Evolving U " << endl;

    //      u_new( dt );
	u_new( Tp , dt2 );

    } // iter loop

    displ=move( Tp , dt );
    
//    update_half_velocity( Tp , false ); 
    update_half_velocity( Tp , overdamped ); 

    areas(Tp);

    quad_coeffs(Tp , simu.FEMp() ); volumes(Tp, simu.FEMp() );

    if(simu.current_step()%simu.every()==0)
      draw(Tp, particle_file , true);

    log_file
      << simu.current_step() << "  "
      <<  simu.time() << "  " ;

    integrals( Tp , log_file);     log_file << "  ";
    fidelity(  Tp , log_file );    log_file << endl;

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
