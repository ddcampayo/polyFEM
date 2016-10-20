//#include <Eigen/Dense>

//#define CHOLMOD

#ifdef CHOLMOD
 #include <Eigen/CholmodSupport>
#else
// #include <Eigen/IterativeLinearSolvers>
 #include <Eigen/SparseCholesky>
//#include <Eigen/SparseQR>
#endif

#include <iostream>
#include <vector>

#include <unsupported/Eigen/SparseExtra>

using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
//using Eigen::ConjugateGradient;

const Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "];");

class linear {
 public:  // constructor
 linear(Triangulation& TT) : T(TT) {}

 private:

  Triangulation& T; // Its triangulation

  typedef   SparseMatrix<double>  SpMat;
  typedef Eigen::Triplet<double> triplet;

  // not inverted.-
  SpMat stiff;
  SpMat lambda_x;
  SpMat lambda_y;

  // inverted.-
  SpMat mass;
  SpMat stiffp1;
  SpMat mas;

  void fill_lambda();
  void fill_stiff();
  void fill_mass();
  void fill_mas( const FT& );

#define DIRECT_SOLVER
#ifdef DIRECT_SOLVER
  #ifdef CHOLMOD
    Eigen::CholmodSupernodalLLT<SpMat> solver_mass;
    Eigen::CholmodSupernodalLLT<SpMat> solver_stiffp1;
    Eigen::CholmodSupernodalLLT<SpMat> solver_mas;
    // Automatic choice of SN vs simplicial
    //Eigen::CholmodDecomposition<SpMat> solver_mass;
    //Eigen::CholmodDecomposition<SpMat> solver_stiffp1;
    //Eigen::CholmodDecomposition<SpMat> solver_mas;
  #else
    Eigen::SimplicialLDLT<SpMat> solver_mass;
    Eigen::SimplicialLDLT<SpMat> solver_stiffp1;
    Eigen::SimplicialLDLT<SpMat> solver_mas;
  #endif
#else
  Eigen::BiCGSTAB<SpMat> solver_mass;
  Eigen::BiCGSTAB<SpMat> solver_stiffp1;
  Eigen::BiCGSTAB<SpMat> solver_mas;
#endif

  //  ConjugateGradient<SpMat> solver_mass;
  //   ConjugateGradient<SpMat> solver_stiffp1;

  //Eigen::SparseQR<SpMat> solver_mass;
  //Eigen::SparseQR<SpMat> solver_stiffp1;

public:

  void ustar_inv(const kind::f Ustar, const FT a, const kind::f U0 , bool ) ;
  void laplacian_v(const kind::f ff1, const kind::f ff2);
  void gradient(const kind::f fsf, const kind::f fvf, bool mass=true );
  void PPE(const kind::f velocity , FT dt,  const kind::f pressure );
  void laplace_div(const kind::f vectorf, const kind::f inter_scalarf, const kind::f scalarf );
  void mass_v(const kind::f vectorf );
  void mass_s(const kind::f scalarf );
  void save_matrices(void);
  void load_matrices(void);


 private:

  void laplacian_stiff(const kind::f ffield, const kind::f gradfield  );
  void laplacian_stiff_v(const kind::f ffield, const kind::f gradfield  );
  void laplacian_Delta_v(const kind::f ffield, const kind::f gradfield  );
  void laplacian_Delta(const kind::f ffield, const kind::f gradfield   );
  void laplacian_s(const kind::f ff1, const kind::f ff2 ) ;
  void div(const kind::f fsf, const kind::f fvf ) ;

  VectorXd vfield_to_vctr(const kind::f vectorf , int comp ) ;
  VectorXd field_to_vctr(const kind::f scalarf    );
  void  vctr_to_field(const VectorXd& vv, const kind::f scalarf ) ;
  void vctr_to_vfield(const VectorXd& vv, const kind::f vectorf, const int comp );


};
