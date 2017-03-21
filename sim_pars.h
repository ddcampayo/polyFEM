struct sim_pars
{
  bool create_points_;
  bool at_random_;
  bool perturb_;
  FT pert_rel_;
  bool FEMp_;
  bool FEMm_;

  int no_of_particles_;
  int no_of_nodes_;
  int mesh_factor_;
  FT h_;

  int Nsteps_;
  int every_;
  int current_step_;

  FT dt_;
  FT time_;

  FT mu_;
  FT v0_;
  FT alpha_;

  sim_pars() {

// do only straight FEM
    FEMm_ = FEMp_ = false;


    create_points_=true;

    no_of_particles_=100;
    mesh_factor_=4;

    no_of_nodes_=    no_of_particles_ * mesh_factor_ * mesh_factor_ ;
 
    at_random_=false;

    h_=0.2;

    perturb_=false;
    pert_rel_=0.1;

    Nsteps_=100;
    every_=10;
    dt_=0.1;
    time_=0;
    current_step_=0;

    mu_=1;
    v0_=0;
  }

public:
  void read() {
    std::ifstream params;
    params.open("simu.cfg");

    std::string word;

    params >> word;

    params >> FEMp_;

    if(FEMp_) cout << "Doing only FEM on particles " << endl;
    else cout << "Doing full quad on particles" << endl;

    params >> word;

    params >> FEMm_;

    if(FEMm_) cout << "Doing only FEM on mesh " << endl;
    else cout << "Doing full quad on mesh" << endl;


    int ii;
    FT ff;

    params >> word;

    params >> ii;
    set_no_of_particles(ii); // sets h_ too!

    cout << "Nominal no. of particles "<< no_of_particles_ << endl;

    params >> word;
    params >> mesh_factor_;

    set_no_of_nodes();

    cout << "Mesh_factor :" << mesh_factor_ << " ; Nominal no. of mesh nodes"<< no_of_nodes_ << endl;

    params >> word;
    params >> create_points_;
    if(create_points_) cout << "Points created anew" << endl;

    params >> word;
    params >>  at_random_;
    if(create_points_) {
      if(at_random_) cout << "Points at random" << endl;
      else  cout << "Points on regular lattice" << endl;
    }
    
    params >> word;
    params >> perturb_;
    params >> pert_rel_;

    if(perturb_)
      cout << "Perturbing points about " << pert_rel_ << endl;


    params >> word;
    params >> dt_;

    cout << "Timestep: " << dt_ << endl;

    params >> word;
    params >> ii;
    if (ii>=0) Nsteps_=ii;
    cout << "Will run " << Nsteps_ << " steps, ";

    params >> word;
    params >> ii;
    if (ii>=0) every_=ii;
    cout << "printing every " << every_ << " steps" << endl;


    params >> word;
    params >> ff;
    if (ii>=0) mu_=ff;
    cout << "Viscosity  " << mu_ << endl;


    params >> word;
    params >> ff;
    if (ff>=0) v0_=ff;
    cout << "Velocity amplitude " << v0_ << endl;

    params >> word;
    params >> ff;
    if (ff>=0) alpha_=ff;
    cout << "Relaxation parameter for iterative procedure " << alpha_ << endl;

    params.close();

    cout << "Reynold's number =" << Re() << endl;
    cout << "Courant's number =" << Co() << endl;
    cout << "Viscous Courant's number =" << Co_mu() << endl;


  }

  bool initial_velocity() const {return  v0_ > 1e-10;}

  FT Re() const {
    return 1.0*v0()/mu();
  }

  FT Co() const {
    return v0()*dt()/h();
  }

  FT Co_mu() const {
    FT hh=h();
    return mu()*dt()/(hh*hh);
  }


  FT v0() const {
    return v0_;
  }

  FT alpha() const {
    return alpha_;
  }


  FT dt() const {
    return dt_;
  }

  FT set_dt(FT ddt)  {
    return dt_=ddt;
  }


  FT mu() const {
    return mu_;
  }


  FT h() const {
    return h_;
  }


  bool create_points() const {return create_points_;}

  int no_of_particles() const {return no_of_particles_;}
  int no_of_nodes() const {return no_of_nodes_;}

  int set_no_of_particles(int nn)  {
    h_=std::sqrt(1.0/FT(nn));
    return no_of_particles_=nn;
  }

  // Dirty! If factor is too large, it is really just the no of nodes
  int set_no_of_nodes(void)  {
    if (mesh_factor_ < 50) {
      int nn = std::sqrt(no_of_particles_) * mesh_factor_ + 0 ;
      return no_of_nodes_= nn*nn;
    } else
      return no_of_nodes_= mesh_factor_;
  }

  int set_no_of_nodes(const int nm)  {
    return no_of_nodes_= nm;
  }

  int mesh_factor() const {return mesh_factor_;}

  bool at_random() const {return at_random_;}
  bool perturb() const {return perturb_;}
  FT pert_rel() const {
		    if(perturb_) return pert_rel_;
		    else return 0.0;
		    }

  bool FEMm() const {return FEMm_;}
  bool FEMp() const {return FEMp_;}

  int Nsteps() const {return Nsteps_;}
  int every() const {return every_;}

  int current_step() const {return current_step_;}

  int next_step()  {return ++current_step_;}

  FT advance_time()  {return  time_+=dt() ;}

  FT time() const  {return  time_ ;}

};

