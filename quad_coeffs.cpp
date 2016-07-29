#include"main.h"
#include"sim_pars.h"
#include"periodic.h"
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

extern sim_pars simu;

using Eigen::VectorXd;
using Eigen::MatrixXd;

void quad_coeffs( Triangulation& T , bool FEM ) {

  for(F_f_it fc=T.finite_faces_begin();
      fc!=T.finite_faces_end();
      fc++) {

    // all this should be in creator, but...
    fc->ll = std::vector<FT>(3,0);

    //    fc->connects= std::vector<My_face_base::v_idx>(3,My_face_base::v_idx);
    //  std::vector<int>
    //    Fb::v_idx vinull;

    typedef std::vector<Vertex_handle> v_vrtcs;
    typedef std::vector<v_vrtcs> v_connects;

    typedef std::vector<FT> v_A;
    typedef std::vector<v_A> v_coeffs;

    fc->connects=  v_connects( 3, v_vrtcs() );

    fc->coeffs= v_coeffs(3, v_A());

  }

  for(F_e_it ed=T.finite_edges_begin();
      ed!=T.finite_edges_end();
      ed++) {

    //    Periodic_segment ps=T.periodic_segment(*ed);

    Segment ss=T.segment(ed);

    FT l2=ss.squared_length();

    //    FT l2=per_dist2(v0->point(), v1->point());

    Face_handle f1=ed->first;
    int i2=ed->second;

    if(!T.is_infinite(f1)) {
      f1->ll[i2]=l2;
    }

    //    Duplicate info on neighboring face.-
    Face_handle f2=f1->neighbor(i2);
    if(!T.is_infinite(f2)) {
      int i2p=T.mirror_index(f1,i2);

      f2->ll[i2p]=l2;
    }
  }

  cout << "Perimeters computed " << endl;

  if ( FEM ) return; // TODO the former should actually be some other function

  for(F_e_it ed=T.finite_edges_begin();
      ed!=T.finite_edges_end();
      ed++) {

    Face_handle f1=ed->first;
    int i2=ed->second;
    int i0=(i2+1)%3;
    int i1=(i0+1)%3;

    std::vector<Vertex_handle> v;

    v.push_back( f1->vertex(i0) );           //0
    v.push_back( f1->vertex(i1) );           //1
    v.push_back( f1->vertex(i2) );           //2
    v.push_back( T.mirror_vertex(f1,i2) );   //3


    /////////////// TODO: loop this up

    Face_circulator fc = T.incident_faces( v[0], f1 ),
      done(fc);

    fc++;

    for( ; ; fc++ ) {
      int ii0=fc->index(v[0]);
      int ii2=(ii0+2)%3;

      if(fc->vertex(ii2)==v[3]) break;

      v.push_back( fc->vertex(ii2)  );

    }

    Face_handle f2=f1->neighbor(i2);
    int i2p=T.mirror_index(f1,i2);

    fc = T.incident_faces( v[1], f2 );
    done=fc;

    fc++;

    for( ;  ; fc++ ) {
      int ii0=fc->index(v[1]);
      int ii2=(ii0+2)%3;

      if(fc->vertex(ii2)==v[2]) break;
      v.push_back( fc->vertex(ii2)  );
    }




    // fc = T.incident_faces( v[2], f1 );
    // done=fc;

    // fc++;

    // for( ;  ; fc++ ) {
    //   int ii0=fc->index(v[2]);
    //   int ii2=(ii0+2)%3;

    //   if(fc->vertex(ii2)==v[0]) break;
    //   v.push_back( fc->vertex(ii2)  );
    // }


    // fc = T.incident_faces( v[3], f2 );
    // done=fc;

    // fc++;

    // for( ;  ; fc++ ) {
    //   int ii0=fc->index(v[3]);
    //   int ii2=(ii0+2)%3;

    //   if(fc->vertex(ii2)==v[1]) break;
    //   v.push_back( fc->vertex(ii2)  );
    // }


    // cout
    //   << "coeffs for edge "
    //   << "( " << v[0]->idx.f()  << " ) -  " 
    //   << "( " << v[1]->idx.f()  << " ) :" << endl;

//    prune(v);

    bool singular=false ;

    if(v.size()==5) {
      singular=true;
      cout << "Warning: only 5 equations" << endl;
    }

    std::vector<Periodic_point> vP;

    for(int i=0; i< v.size() ;i++)
      vP.push_back( T.periodic_point(v[i]) );

    std::vector<FT> A( v.size() );

    // control against singular values is probably superceded by the
    // correct treatment for condition number at solve_linear
    FT result=-1;

    if(!singular) {
      result=solve_linear(vP,A);
      singular=(result > 1e4); // probably singular matrix (???)
    }
    if(singular)
      cout << "Warning: ill conditioned problem , result =" << result << endl;

    // cout << "Large A check" << endl;

    for(int i=0;i<A.size();i++) {
      	// cout << A[i] << endl;
      if(fabs(A[i]) > 1e2) {
	singular=true;
        cout << "Warning: large coeffient found: " 
             << i << " : " << A[i] << endl;

// negligible improvement:
//	break;
      }
    }

    if(singular) {
      cout <<  "Skippable.  ";
      cout <<  "dgesv result: " << result;
      cout << "  .  Degs: ";
      for(int j=0;j<v.size();j++) 
	cout << T.degree(v[j]) << "  ";
      cout << endl;
      cout << ".  Vectors: " << endl;
      
      for(int j=0;j<vP.size();j++) {
	  
	  // Periodic_segment ps={{vP[0],vP[j]}};
	  // Segment ss=T.segment(ps);
	  // Vector_2 v=ss.to_vector();
	Vector_2 v=per_vect(vP[0].first,vP[j].first);

	cout << v.x() << "  " << v.y() << endl;
      }

      cout << endl;

	//	// skip actually:
      continue;
    }


    for(int i=0; i< v.size() ;i++) {
	  // cout
	  //   << v[i]->indx() << "  -->  "
	  //   << A[i]  << endl;

      f1->connects[i2].push_back(v[i]);
      f1->coeffs[i2].push_back(A[i]);

	  //Duplicate.-
      f2->connects[i2p].push_back(v[i]);
      f2->coeffs[i2p].push_back(A[i]);

    }

    // cout
    //   << f1->connects[i2].size() << " , "
    //   << f1->coeffs[i2].size() <<  " f2: "
    //   << f2->connects[i2p].size() << " , "
    //   << f2->coeffs[i2p].size() << endl;


  }

  cout << "Quad coeffs computed " << endl;

  return;
}





FT solve_linear(std::vector<Periodic_point>& vP,
		  std::vector<FT>& AA ) {

  const int m=6;
  int n=vP.size();

  VectorXd X=VectorXd::Zero( m );  // check

  Vector_2 v01=per_vect(vP[0].first,vP[1].first);

  FT dx2 = v01.x()*v01.x();
  FT dy2 = v01.y()*v01.y();
  FT dxy = v01.x()*v01.y();

  FT d2 = dx2 + dy2;
  FT dd = std::sqrt(d2);

  X(3)= -dx2/d2;
  X(4)= -dy2/d2;
  X(5)= -dxy/d2;


  MatrixXd M=MatrixXd::Zero( m, n );

  for(int j=0;j<n;j++) {
    //    M.set(0, j, 1);

    Vector_2  v0j=per_vect(vP[0].first,vP[j].first);
    {
      // Periodic_segment ps={{vP[0],vP[j]}};
      // Segment ss=T.segment(ps);
      // v0j=ss.to_vector();
    }

    Vector_2  v1j=per_vect(vP[1].first,vP[j].first); 
    {
      // Periodic_segment ps={{vP[1],vP[j]}};
      // Segment ss=T.segment(ps);
      // v1j=ss.to_vector();

    }
	
    FT x0=v0j.x();
    FT y0=v0j.y();

    FT x1=v1j.x();
    FT y1=v1j.y();

    M(0, j)= 1 ;
    M(1, j)= x0 / dd ;
    M(2, j)= y0 / dd;
    M(3, j)= x0*x1 / d2;
    M(4, j)= y0*y1 / d2;
    M(5, j)= x0*y1 / d2;

  }

  // M & X are modified! careful...
  // cout << "M = "     << M.format(OctaveFmt) << endl;
  // cout << "X = "     << X.format(OctaveFmt) << endl;

  Eigen::JacobiSVD<MatrixXd> svd(M, Eigen::ComputeThinU | Eigen::ComputeThinV);

  VectorXd svs= svd.singularValues();
  // cout << "Singular values are:" << endl << svs << endl;

  svd.setThreshold(1e-6);
  int rank=svd.rank();

  FT cn=  svs(0) / svs(5) ;

  // cout << "Rank: " << rank
  //      << " ; condition n: " << cn
  //      <<   " ; threshold :" << svd.threshold()
  //      << endl;

  // if(
  //    ( (rank==3) && (cn<5.0) ) ||   // "good" equations
  //    (nn<6)                         // not too picky: single solution (least-square for 3 and 4).
  //    ) {


  VectorXd A = svd.solve(X);

  // cout << "A = "     << A.format(OctaveFmt) << endl;

  for(int i=0; i< n ;i++) AA[i]=A(i);

  return cn;
}


