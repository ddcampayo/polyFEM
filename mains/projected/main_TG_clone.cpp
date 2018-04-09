// plus, correction to quadratic consistency
// periodic boundary conditions

// particle re-creation!


#include <CGAL/Timer.h>

// write out matrices
//#define WRITE

//  #define EXPLICIT

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


// with mass matrix, and all
// stuff only used here:

//#define FULL
#define FULL_FULL
//#define FULL_LUMPED
//#define FLIP

///////////


#ifdef FULL_FULL
#define FULL
#endif

#ifdef FULL_LUMPED
#define FULL
#endif


#include"onto_from_mesh.h"

Triangulation Tm(domain); // mesh

void create(Triangulation&);
void clone(const Triangulation&,Triangulation&);


int main() {

//  CGAL::Timer time;
//
//  time.start();
  
  cout << "Creating point cloud" << endl;

  simu.read();

  create( Tm );

  if(simu.create_points()) {
    set_fields_TG(Tm);
    number(Tm);
  }



  // just once!
  linear algebra(Tm);

  areas(Tm);
  quad_coeffs(Tm , simu.FEMm() ); volumes(Tm, simu.FEMm() );

  if(simu.create_points()) {
    nabla(Tm);
    Delta(Tm);
  }
  
  move_info(Tm);

  const std::string mesh_file("mesh.dat");
  const std::string particle_file("particles.dat");

  draw(Tm, mesh_file     , true);
  
  simu.advance_time();
  simu.next_step();

  //  bool first_iter=true;

  CGAL::Timer time;

  time.start();

  std::ofstream log_file;

  log_file.open("main.log");

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

    move_info(Tm);

    Triangulation Tp(domain); // particles
    clone(Tm,Tp);
    
    // FLIP idea: set initial increment to null
    //    reset_v( Tm , kind::DELTAU);
    
    for( ; iter<max_iter ; iter++) {

      zero_mean_v( Tp, kind::U);


      // comment for no move.-
      cout << "Moving half step " << endl;
      FT d0;

      displ = move( Tp , dt2 , d0 );

      //      cout << "Moved half step " << endl;

      areas(Tp);
      quad_coeffs(Tp , simu.FEMp() ); volumes(Tp, simu.FEMp() );
       
      cout
	<< "Iter " << iter
	<< " , moved avg " << d0 << " to half point, "
	<< displ << " from previous"
	<< endl;

      
      if( displ < min_displ) {
	min_displ=displ;
	min_iter=iter;
      }

      cout << "iter " << iter << "  ; relative displacement:  " << displ << endl;

      if(displ < max_displ) break;

      cout << "Proj U0 onto mesh " << endl;

#if defined FULL
      onto_mesh_full_v(Tp,Tm,algebra,kind::UOLD);
#elif defined FLIP
      flip_volumes(Tp , Tm , simu.FEMm() );
      onto_mesh_flip_v(Tp,Tm,simu.FEMm(),kind::UOLD);
#else
      onto_mesh_delta_v(Tp,Tm,kind::UOLD);
#endif

      //      Tm.transfer_v(

      cout << "Calculating Ustar implicitely" << endl;

      algebra.ustar_inv(kind::USTAR,  dt2 , kind::UOLD,  false, false);
//      algebra.ustar_inv(kind::USTAR,  0 , kind::UOLD,  false, false);
      
      cout << "Solving PPE" << endl;
      
      algebra.PPE( kind::USTAR, dt2 , kind:: P );

      cout << "Calculating grad p" << endl;
      algebra.gradient(kind::P, kind::GRADP);

      cout << "Evolving U " << endl;

    //      u_new( dt );
      u_new( Tm , dt2 );


      cout << "Proj Delta U from mesh " << endl;

      // FLIP idea: tranfer velocity _increment_ to particles
      
#if defined FULL_FULL
      {
	Delta(Tp);
	linear algebra_p(Tp);
	from_mesh_full_v(Tm, Tp, algebra_p , kind::DELTAU);
	//	from_mesh_full_v(Tm, Tp, algebra_p , kind::UOLD);
      }
#elif defined FULL_LUMPED
      from_mesh_lumped_v(Tm, Tp, kind::DELTAU);
      //      from_mesh_lumped_v(Tm, Tp, kind::UOLD);
#elif defined FLIP
      from_mesh_v(Tm, Tp, kind::DELTAU);
      //      from_mesh_v(Tm, Tp, kind::UOLD);
#else
      from_mesh_v(Tm, Tp, kind::DELTAU);
      //      from_mesh_v(Tm, Tp, kind::UOLD);
#endif

      // FLIP idea: apply velocity _increment_ to particles

      incr_v( Tp ,  kind::UOLD , kind::DELTAU , kind::U );

    }

    zero_mean_v( Tp, kind::U);

    cout << "Moving whole step: relative ";

    FT d0;
    
    displ=move( Tp , dt , d0 );

    cout
      <<  displ << " from half point, "
      <<  d0    << " from previous point"
      << endl;
    
      // comment for no move.-

    update_half_velocity( Tp , false ); 

    //    update_half_velocity( Tm , false ); 

    zero_mean_v( Tp, kind::U);

    areas(Tp);

    quad_coeffs(Tp , simu.FEMp() ); volumes(Tp, simu.FEMp() );

    cout << "Proj U_t+1 onto mesh " << endl;

#if defined FULL
    onto_mesh_full_v(Tp,Tm,algebra,kind::U);
#elif defined FLIP
    flip_volumes(Tp , Tm , simu.FEMm() );
    onto_mesh_flip_v(Tp,Tm,simu.FEMm(),kind::U);
#else
    onto_mesh_delta_v(Tp,Tm,kind::U);
#endif

    
    
    move_info( Tm );
    move_info( Tp );

    if(simu.current_step()%simu.every()==0)
      {
	draw(Tm, mesh_file     , true);
	draw(Tp, particle_file , false);
      }

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



void create(Triangulation& Tp ) {

  int N=simu.no_of_particles();
  std::vector<Point> points;
  //  points.reserve(N);

  if(simu.create_points()) {
    if(simu.at_random()) {
      points.reserve(N);

      CGAL::Random_points_in_square_2<Point,Creator> g(LL/2.0-0.0001);
      CGAL::cpp11::copy_n( g, N, std::back_inserter(points));


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

      //      for(int i = 0 ; i < Nb ; ++i )
	
      
      if( false ) {
	CGAL::perturb_points_2(
			       points.begin(), points.end(),
			       simu.pert_rel()* spacing );//,Creator());
	cout << "each particle perturbed about " << simu.pert_rel()* spacing  << endl;
      }

    }

   }

    cout << "Inserting" << endl;

    Tp.insert(points.begin(), points.end());

    return;

}


void clone(const Triangulation& Tfrom,Triangulation& Tto) {

  for(F_v_it vit=Tfrom.vertices_begin();
      vit != Tfrom.vertices_end();
      vit++) {

    Periodic_point pp=Tfrom.periodic_point(vit);
    Point p=Tfrom.point(pp);

    Vertex_handle fv=Tto.insert( p );

    fv->rold.set( p );

    fv->U.set( vit->U.val() );
    fv->Uold.set( vit->U.val() );

    fv->idx.set( vit->idx.val() );

  }      

  return;
  
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
