// Compilation:
// g++ -I ./fftwpp -fopenmp CH_exp.cc fftwpp/fftw++.cc -lfftw3 -lfftw3_omp

#include "fftwpp/Complex.h"
#include "fftwpp/Array.h"
#include "fftwpp/fftw++.h"


typedef unsigned int uint;

typedef double FT;

typedef Array::array1<Complex> c_vector;
typedef Array::array2<Complex> c_array;

typedef Array::array1<FT> r_vector;
typedef Array::array2<FT> r_array;

  

class CH_FFT {
  
  //#include "utils.h"

  //using namespace utils;
  //using namespace fftwpp;
  //using fftw::maxthreads;
  //  using fftwpp::fftw; //::get_max_threads;

  FT L; // size of system
  uint nx,ny; // size of system

  // TODO: read this
  static constexpr FT L_SD = 62.5;  // realistic value
  //static constexpr FT L_SD = 1.25;  // hydro regime
  const static  bool shift=false;
  const static bool quiet=false;
  const size_t align;
  const uint threads;



 public:

  // constructors
  // Huge list, because fftwpp::fft2d lacks a default constructor!

 CH_FFT(const FT L_, const uint nnx) :
  L(L_), nx(nnx), ny(nnx),
    freq_pref( 2 * M_PI / L ),
    align(sizeof(Complex)),
    fr( nx,ny,align ),
    f( nx,ny,align ),
    f3r( nx,ny,align ),
    f3( nx,ny,align ),
    mu_r( nx,ny,align ),
    mu( nx,ny,align ),
    press_r( nx,ny,align ),
    press( nx,ny,align ),
    grad_mu_x_r( nx,ny,align ),
    grad_mu_x( nx,ny,align ),
    grad_mu_y_r( nx,ny,align ),
    grad_mu_y( nx,ny,align ),
    force_x_r( nx,ny,align ),
    force_x( nx,ny,align ),
    force_y_r( nx,ny,align ),
    force_y( nx,ny,align ),
    v_x_r( nx,ny,align ),
    v_x( nx,ny,align ),
    v_y_r( nx,ny,align ),
    v_y( nx,ny,align ),
    q2(nx,ny,align),
    q_x(nx,  align),
    q_y(nx,  align),
    threads( get_max_threads() ),
    Forward2_f( -1 , fr , f , threads),
    Backward2_f( 1,   f  , fr , threads),
    Forward2_f3( -1 , fr , f , threads),
    Backward2_f3( 1,   f  , fr , threads),
    Forward2_mu( -1 , fr , f , threads),
    Backward2_mu( 1,   f  , fr , threads),
    Forward2_gmx( -1 , fr , f , threads),
    Backward2_gmx( 1,   f  , fr , threads),
    Forward2_gmy( -1 , fr , f , threads),
    Backward2_gmy( 1,   f  , fr , threads),
    Forward2_forcex( -1 , fr , f , threads),
    Backward2_forcex( 1,   f  , fr , threads),
    Forward2_forcey( -1 , fr , f , threads),
    Backward2_forcey( 1,   f  , fr , threads),
    Forward2_vx( -1 , fr , f , threads),
    Backward2_vx( 1,   f  , fr , threads),
    Forward2_vy( -1 , fr , f , threads),
    Backward2_vy( 1,   f  , fr , threads)    
    {
      init();
    }

  void printout(const std::string& name, const c_array& f );
  void draw(const std::string& name, const int time, const c_array& f ) ;
  void printout_x(const std::string& name, const unsigned int nx);
  void random( void );
  void all_fields(void);
  void all_fields_NS(const FT&);
  void vel_fields_NS(void);
  void p_fields( const FT& ) ;
  void mu_fields(void);
  void force_fields(void);
  void vel_fields(void);
  void evolve( const FT );
  void histogram(const std::string& name, const int time, const c_array& ff );

  void set_f(const uint i,const uint j, const FT val) {
    fr(i,j).re = val;
    fr(i,j).im = 0;
    return;
  }

  void set_f(const c_array& ); 
  void set_u(const c_array& , const c_array& );
  
  c_array field_f() const {return fr;}
  c_array field_fq() const {return f;}
  c_array field_mu() const {return mu_r;}
  c_array field_p() const {return press_r;}
  c_array field_grad_mu_x() const {return grad_mu_x_r;}
  c_array field_grad_mu_y() const {return grad_mu_y_r;}
  c_array field_force_x() const {return force_x_r;}
  c_array field_force_y() const {return force_y_r;}
  c_array field_vel_x() const {return v_x_r;}
  c_array field_vel_y() const {return v_y_r;}

  int Nx() const {return nx;}
  size_t alignment() const {return align;}

 private:

  FT freq_pref;

  struct hist_data{
    FT f;   // phi value
    FT q2;  // q^2 value
    long int count;  // frequency
  };
  
  std::map<long int,hist_data> histo;
  //  std::map<int,int> h_count;

  void cube( ) ;
  void blank( );
  void delta( );
  void gaussian( );
  void sine_wave( );
  void freq_scramble(c_array& f);
  void init( void );

  c_array fr; // real space
  c_array f;  // Fourier space
  
  c_array f3r;
  c_array f3;

  c_array press_r;
  c_array press;

  c_array mu_r;
  c_array mu;

  c_array grad_mu_x_r;
  c_array grad_mu_x;

  c_array grad_mu_y_r;
  c_array grad_mu_y;

  c_array force_x_r;
  c_array force_x;

  c_array force_y_r;
  c_array force_y;
  
  c_array v_x_r;
  c_array v_x;
  c_array v_y_r;
  c_array v_y;

  r_array q2;
  //  r_array q4(nx,ny,align);

  r_vector q_x;
  r_vector q_y;
 
      //  // Pointers, because these don't have void constructors !

  //  fftwpp::fft2d*  Forward2_f;
  //  fftwpp::fft2d* Backward2_f;

  /* fftwpp::fft2d  Forward2; */
  /* fftwpp::fft2d Backward2; */
  
  fftwpp::fft2d  Forward2_f;
  fftwpp::fft2d Backward2_f;

  fftwpp::fft2d  Forward2_f3;
  fftwpp::fft2d Backward2_f3;

  fftwpp::fft2d  Forward2_mu;
  fftwpp::fft2d Backward2_mu;

  fftwpp::fft2d  Forward2_gmx;
  fftwpp::fft2d Backward2_gmx;

  fftwpp::fft2d  Forward2_gmy;
  fftwpp::fft2d Backward2_gmy;

  fftwpp::fft2d  Forward2_forcex;
  fftwpp::fft2d Backward2_forcex;

  fftwpp::fft2d  Forward2_forcey;
  fftwpp::fft2d Backward2_forcey;
  
  fftwpp::fft2d  Forward2_vx;
  fftwpp::fft2d Backward2_vx;

  fftwpp::fft2d  Forward2_vy;
  fftwpp::fft2d Backward2_vy;



};
