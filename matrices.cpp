void matrices(void) {

  // Quick check up.-
  for(F_f_it fc=T.finite_faces_begin();
      fc!=T.finite_faces_end();
      fc++) {
    std::vector<Vertex_handle> v(3);

    for(int i=0;i<3;i++) 
      v[i]=fc->vertex(i);
    
    for(int i2=0;i2<3;i2++) {
      int i0=(i2+1)%3;
      int i1=(i0+1)%3;

      cout
	<< "edge ( "
	<< v[i0]->indx() << " ,  "
	<< v[i1]->indx() << " )  involves  ";
      for(int i=0; i< fc->connects[i2].size() ;i++) 
	cout << fc->connects[i2][i] << " , ";

      cout << endl;
   }

  }

  std::ofstream matrices;
  matrices.open("matrices.m");

  matrices.precision(15);

  matrices << "x = [" << endl;
  for(F_v_it vit=T.finite_vertices_begin();
      vit != T.finite_vertices_end();
      vit++)
    matrices << vit->point().x() << endl;
  matrices << "];" << endl;

  matrices << "y = [" << endl;
  for(F_v_it vit=T.finite_vertices_begin();
      vit != T.finite_vertices_end();
      vit++)
      matrices << vit->point().y() << endl;
  matrices << "];" << endl;

  MatrixXd mass  = MatrixXd::Zero(Nin,Nin);
  MatrixXd stiff = MatrixXd::Zero(Nin,Nin);
  VectorXd Vol   = VectorXd::Zero(Nin); // FEM volumes
  MatrixXd Delta = MatrixXd::Zero(Nin,Nin); // approx to Laplacian

  // aliases
  const FT& dd=Delta_par.d;
  const FT& db=Delta_par.db;
  const FT& dD=Delta_par.D;

  for(F_f_it fc=T.finite_faces_begin();
      fc!=T.finite_faces_end();
      fc++) {

    std::vector<Vertex_handle> v(3);

    for(int i=0;i<3;i++) 
      v[i]=fc->vertex(i);

    FT a=fc->area;

    for(int i0=0; i0< 3 ;i0++) {
      // self.-
      int idx0=v[i0]->indx();

      //      cout << idx0 << " gets area = " << a << endl;

      FT l0=fc->ll[i0];

      int i1=(i0+1)%3;
      int i2=(i1+1)%3;

      FT l1=fc->ll[i1];
      FT l2=fc->ll[i2];

      int idx1=v[i1]->indx();
      int idx2=v[i2]->indx();

      FT mm=a/6.0;
      FT ss= l0/a/4.0;

      mass(idx0,idx0)+= mm;
      stiff(idx0,idx0)+= ss;

      FT del= dd * l0/a/2.0;

      Delta(idx0,idx0)+= del;
      //      Delta(idx0,idx0)+= ss;

      mm=a/12.0;
      ss=(l2- l0 -l1 )/a/8.0;

      mass(idx0,idx1)+=mm;
      mass(idx1,idx0)+=mm;

      stiff(idx0,idx1)+= ss;
      stiff(idx1,idx0)+= ss;

      del = dd * (l2- l0 -l1 )/a/4.0;

      Delta(idx0,idx1)+=  del;
      Delta(idx1,idx0)+=  del;

      Vol(idx0) += 2*a;

      if(simu.FEM()) continue;

      Fb::v_idx indices=fc->connects[i0];
      Fb::v_A coeffs=fc->coeffs[i0];
      
      cout
	<< "edge ( "
	<< idx0 << " ,  "
	<< idx1 << " )  involves by quad : " << endl;
   
      cout << "triangle: ( " 
	   << idx0 << " , "
	   << idx1 << " , "
	   << idx2 << " )  " << endl;

      for(int n=0;n<indices.size();n++) {

	FT mm,ss;
// coeffs, order 1

	int idx=indices[n];

	cout  << idx << " x [ ";

	FT A=coeffs[n];

	// this could be made as a loop over vertices, but

	mm=  A* a/60;

	mass(idx0,idx)+=mm;
	  //	if(idx0!=idx)
	mass(idx,idx0)+=mm;

	ss= -A*l0/(12*a);

	stiff(idx0,idx)+=ss;
      //	if(idx0!=idx)
	stiff(idx,idx0)+=ss;


	del = A/(2*a)*(
		       (dD-db)*(l1+l2)
		       -l0*dD);

	Delta(idx0,idx)+=  del;
      //	if(idx0!=idx)
	//	Delta(idx,idx0)+=ss;

	mm= A* a/30.0;

	mass(idx1,idx)+=mm;
	  //	if(idx1!=idx)
	mass(idx,idx1)+=mm;

	ss= A*(l1+l0-l2)/(24*a);
      
	stiff(idx1,idx)+=ss;
	  //	if(idx1!=idx)
	stiff(idx,idx1)+=ss;

	del = A/(2*a)*(
		       (2*dD-dd)*(l1+l2-l0)/2
		       +l1*db);

	Delta(idx1,idx)+=  del;
	  //	if(idx1!=idx)
	//	Delta(idx,idx1)+=ss;

	mm= A* a/30.0;

	mass(idx2,idx)+=mm;
	  //	if(idx2!=idx)
	mass(idx,idx2)+=mm;

	ss= A*(l2+l0-l1)/(24*a);
	  
	stiff(idx2,idx)+=ss;
	  //	if(idx2!=idx)
	stiff(idx,idx2)+=ss;

	del = A/(2*a)*(
		       (2*dD-dd)*(l1+l2-l0)/2
		       +l2*db);

	Delta(idx2,idx) += del;
	  //	if(idx2!=idx)
	//	Delta(idx,idx2)+=ss;


	if(fabs(mm)>1000) {
	  cout << idx << " ( " ;
	  cout << mm << " , ";
	  cout << ss << " ) " << endl;
	}

	  //	  cout  << endl << "    ";

// coeffs, order 2
// self edge
	for(int m=0;m<indices.size();m++) {

	  int idxm=indices[m];

	  cout  << idxm << " , ";
	    
	  FT Am=coeffs[m];

	  mm= A*Am* a/90;
	  ss= A*Am*(l0+l1+l2)/(48*a);

	  if(fabs(mm)>1000) {
	    cout <<  "        ( " ;
	    cout << idxm << " , ";
	    cout << mm << " , ";
	    cout << ss << " ) " << endl;
	  }

	  mass(idx,idxm)+=mm;
	  //if(idx!=idxm)
//	    mass(idxm,idx)+=mm;
	  stiff(idx,idxm)+=ss;
	  //if(idx!=idxm)
//	    stiff(idxm,idx)+=ss;
	}

	cout  << " ] & [ ";

// next edge
	Fb::v_idx indices_n=fc->connects[i1];
	Fb::v_A coeffs_n=fc->coeffs[i1];

	for(int m=0;m<indices_n.size();m++) {

	  int idxm=indices_n[m];

	  cout  << idxm << " , ";

	  FT Am=coeffs_n[m];

	  mm= A*Am* a/180.0;
	  ss= A*Am*(l2-l0-l1)/(48*a);

	  if(fabs(mm)>1000) {
	    cout <<  "        ( " ;
	    cout << idxm << " , ";
	    cout << mm << " , ";
	    cout << ss << " ) " << endl;
	  }

	  mass(idx,idxm)+=mm;
	   //	   if(idx!=idxm)
	  mass(idxm,idx)+=mm;

	  stiff(idx,idxm)+=ss;
	    //	   if(idx!=idxm)
	  stiff(idxm,idx)+=ss;
	}
	cout << " ]" << endl;
      }

      cout << endl;

    }
  }



  matrices << "m = " <<  mass.format(OctaveFmt) << endl;
  matrices << "s = " <<  stiff.format(OctaveFmt) << endl;


  matrices << "V = " <<  Vol.format(OctaveFmt) << endl;

  for(int i=0 ; i<Nin ; i++) 
    Delta.row(i) /= (Vol(i)/1.0);

  matrices << "D = " <<  Delta.format(OctaveFmt) << endl;

  matrices.close();

}

