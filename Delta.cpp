#include"main.h"
#include"sim_pars.h"


//extern Triangulation T;
extern sim_pars simu;

// Laplacian couplings

// Parameters for the construction of Delta
struct dp {

 dp() {
    D  = 1; // leave this one set to 1...
    d  = 3; // relative to D
    db = 1;
  }

  FT D,d,db;

} Delta_par;


void Delta(Triangulation& T) {

  // MatrixXd mass  = MatrixXd::Zero(Nin,Nin);
  // MatrixXd stiff = MatrixXd::Zero(Nin,Nin);
  // VectorXd Vol   = VectorXd::Zero(Nin); // FEM volumes
  // MatrixXd Delta = MatrixXd::Zero(Nin,Nin); // approx to Laplacian

  // aliases

  // FT dd=Delta_par.d;
  // if(simu.FEM()) dd=1.0/2.0;

  const FT& dd=Delta_par.d/6 ;
  const FT& db=Delta_par.db/6;
  const FT& dD=Delta_par.D/6 ;

  for(F_f_it fc=T.finite_faces_begin();
      fc!=T.finite_faces_end();
      fc++) {

    std::vector<Vertex_handle> v(3);

    for(int i=0;i<3;i++) 
      v[i]=fc->vertex(i);

    FT a=fc->area;

    for(int i0=0; i0< 3 ;i0++) {
      // self.-
      //      int idx0=v[i0]->idx.f();

      //      cout << idx0 << " gets area = " << a << endl;

      FT l0=fc->ll[i0];

      int i1=(i0+1)%3;
      int i2=(i1+1)%3;

      FT l1=fc->ll[i1];
      FT l2=fc->ll[i2];

      //      int idx1=v[i1]->indx();
      //      int idx2=v[i2]->indx();

      FT mm=a/6.0;
      FT ss= l0/a/4.0;

      v[i0]->mass(v[i0],mm);

      v[i0]->stiff(v[i0],ss);

      FT del= dd * l0/a/2.0;

      //      Delta(idx0,idx0)+= del;

      v[i0]->Delta(v[i0],del);


      //      Delta(idx0,idx0)+= ss;

      mm=a/12.0;
      ss=(l2- l0 -l1 )/a/8.0;

      v[i0]->mass(v[i1],mm);
      v[i1]->mass(v[i0],mm);

      v[i0]->stiff(v[i1],ss);
      v[i1]->stiff(v[i0],ss);

      del = dd * (l2- l0 -l1 )/a/4.0;

      v[i0]->Delta(v[i1],del);
      v[i1]->Delta(v[i0],del);

      // Delta(idx0,idx1)+=  del;
      // Delta(idx1,idx0)+=  del;

      //      v[i0]->Vol( 2*a );

      if(simu.FEMm()) continue;

      typedef std::vector<Vertex_handle> v_vrtcs;
    //    typedef std::vector<v_vrtcs> v_connects;

      typedef std::vector<FT> v_A;
    //    typedef std::vector<v_A> v_coeffs;

      v_vrtcs  cvs=fc->connects[i0]; // cv: connected vertices
      v_A coeffs=fc->coeffs[i0];

      for(int n=0;n<cvs.size();n++) {

	FT mm,ss;
// coeffs, order 1


	FT A=coeffs[n];
	Vertex_handle cv=cvs[n];

	// this could be made as a loop over vertices, but

	mm=  A* a/60;

	v[i0]->mass(cv,mm);
	cv->mass(v[i0],mm);

	ss= -A*l0/(12*a);

	v[i0]->stiff(cv,ss);
	cv->stiff(v[i0],ss);

	del = A/(2*a)*(
		       (dD-db)*(l1+l2)
		       -l0*dD);

	v[i0]->Delta(cv,del);
	cv->Delta(v[i0],del);

	mm= A* a/30.0;

	v[i1]->mass(cv,mm);
	cv->mass(v[i1],mm);

	ss= A*(l1+l0-l2)/(24*a);

	v[i1]->stiff(cv,ss);
	cv->stiff(v[i1],ss);      

	del = A/(2*a)*(
		       (2*dD-dd)*(l1+l2-l0)/2
		       +l1*db);

	v[i1]->Delta(cv,del);
	cv->Delta(v[i1],del);

	mm= A* a/30.0;

	v[i2]->mass(cv,mm);
	cv->mass(v[i2],mm);

	ss= A*(l2+l0-l1)/(24*a);

	v[i2]->stiff(cv,ss);
	cv->stiff(v[i2],ss);      

	del = A/(2*a)*(
		       (2*dD-dd)*(l1+l2-l0)/2
		       +l2*db);

	v[i2]->Delta(cv,del);
	cv->Delta(v[i2],del);

	// if(fabs(mm)>1000) {
	//   cout << idx << " ( " ;
	//   cout << mm << " , ";
	//   cout << ss << " ) " << endl;
	// }

	  //	  cout  << endl << "    ";

// coeffs, order 2
// self edge
	for(int m=0;m<cvs.size();m++) {

	  Vertex_handle cv2=cvs[m];

	  FT Am=coeffs[m];

	  mm= A*Am* a/90;
	  ss= A*Am*(l0+l1+l2)/(48*a);

	  // if(fabs(mm)>1000) {
	  //   cout <<  "        ( " ;
	  //   cout << idxm << " , ";
	  //   cout << mm << " , ";
	  //   cout << ss << " ) " << endl;
	  // }

	  cv->mass(cv2,mm);
	  cv->stiff(cv2,ss);

	}


// next edge

	v_vrtcs  cvs2=fc->connects[i1]; // cv: connected vertices
	v_A coeffs2=fc->coeffs[i1];

	for(int m=0;m<cvs2.size();m++) {

	  Vertex_handle cv2=cvs2[m];

	  FT Am=coeffs2[m];

	  mm= A*Am* a/180.0;
	  ss= A*Am*(l2-l0-l1)/(48*a);

	  cv->mass(cv2,mm);
	  cv2->mass(cv,mm);

	  cv->stiff(cv2,ss);
	  cv2->stiff(cv,ss);

	}
	//	cout << " ]" << endl;
      }

      //      cout << endl;

    }
  }

}
