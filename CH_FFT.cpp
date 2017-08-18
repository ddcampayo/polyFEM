// Compilation:
// g++ -I ./fftwpp -c CH_FFT.cpp
// // g++ -I ./fftwpp -fopenmp -c CH_FFT.cpp  fftwpp/fftw++.cc -lfftw3 -lfftw3_omp

#include "CH_FFT.h"

  //using namespace std;
using std::cout;
using std::endl;

void CH_FFT::blank( void )
{
  for(unsigned int i=0; i < f.Nx(); ++i)
    for(unsigned int j=0; j < f.Ny(); ++j)
      f(i,j).re = f(i,j).im = 0;

  return;
}

void CH_FFT::printout_x(const std::string& name, const unsigned int nx)
{
  cout << name << " = np.matrix( ' ";

  for(unsigned int i=0; i < nx; ++i)

    cout << L * ( FT(i)/nx - 0.5 ) << " ";

  cout << " ' ) " << endl;

  return;
}


void CH_FFT::delta( void )
{

  unsigned int nx = f.Nx() ;
  unsigned int ny = f.Ny() ;

  unsigned int centerx = nx / 2 ;
  unsigned int centery = ny / 2 ;

  FT strength = nx * ny ; // * (L*L); // deltas are special
  
  f(centerx,centery) = strength;
  return;
}


void CH_FFT::gaussian( void )
{

  unsigned int nx = f.Nx() ;
  unsigned int ny = f.Ny() ;

  const FT w= 0.1 * L;
  const FT w2= w*w;
  
  const FT norm= 1e-3;// 1.0/(2*M_PI*w2);

  for(unsigned int i=0; i < nx; ++i) {

    FT y = L * ( 0.5 - i /FT(nx) );

    for(unsigned int j=0; j < ny; ++j) {
      FT x = L * ( j/FT(ny) - 0.5 );

      // cout << i << " , " << j << " --- " << x << " , " << y << endl;
      f(i,j).re =  std::exp(  -  ( x * x + y * y) / ( 2 * w * w ) ) ;
      f(i,j).im = 0;
    }
  }
  return;
}




void CH_FFT::sine_wave( void )
{
  unsigned int nx = f.Nx() ;
  unsigned int ny = f.Ny() ;

  for(unsigned int i=0; i < nx; ++i) {
    FT y=i/FT(nx);

    for(unsigned int j=0; j < ny; ++j) {
      FT x= j /FT(ny);


      const FT k0= 2*M_PI ;
      const FT kx= 1 * k0;
      const FT ky= 0 * k0;

      const FT off=0;

      const FT ampl=1e-3;
      //int par = i + j;

      //int sign = ( par%2 == 0 ? 1 : -1 );
      
      //      f(i,j).re = off+ sign * std::sin(  kx * x + ky * y ) ;
      f(i,j).re = off + ampl * std::sin(  kx * x + ky * y ) ;
      f(i,j).im = 0;
    }
  }
}




void CH_FFT::freq_scramble( c_array& f)
{
  unsigned int nx = f.Nx() ;
  unsigned int ny = f.Ny() ;

  for(unsigned int i=0; i < nx; ++i) {

    for(unsigned int j=0; j < ny; ++j) {
      
      int par = i + j;

      if ( par%2 == 0 ) continue;
      
      f(i,j) *= -1;
    }
  }

  return;
}



void CH_FFT::draw(const std::string& name, const int time, const c_array& ff ) {

  std::stringstream  namefile;
  namefile << "field_" << name << "_" << time << ".dat";

  //  cout << "writing on file : " << namefile.str() << endl;

  std::ofstream main_data;

  main_data.open(namefile.str().c_str() );

  for(unsigned int i=0; i < nx ; ++i) {

    for(unsigned int j=0; j < ny; ++j) 
      main_data << real( ff( i,j ) ) << "  ";

    main_data << endl;
  }

  main_data.close();

  return;
}



// f3 = f^3. Just for reals currently

void CH_FFT::cube( void ) {
   
  for(unsigned int i=0; i < nx ; ++i)
    for(unsigned int j=0; j < ny; ++j)
      f3r(i,j) = fr(i,j) * ( fr(i,j) * fr(i,j) ) ;

  return;
	
}

void CH_FFT::printout( const std::string& name, const c_array& f ) {
   
  cout << name << " = np.matrix( ' ";//  << endl ;

  for(unsigned int i=0; i < nx ; ++i) {
    //      cout << " , " ;

    for(unsigned int j=0; j < ny; ++j) {
      cout << real(f(i,j)) << "  ";
    }

    if(i != nx - 1)  cout << " ; " ;//<< endl;
  }
  cout << " ') " << endl;

  cout << endl;

  cout << name << "_im = np.matrix( ' ";//  << endl ;

  for(unsigned int i=0; i < nx; ++i) {
    //      cout << " , " ;

    for(unsigned int j=0; j < ny; ++j) {
      cout << imag(f(i,j)) << "  ";
    }

    if(i != nx - 1)  cout << " ; " ;//<< endl;
  }
  cout << " ') " << endl;

  cout << endl;

  return;

}




#include <ctime>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>



void CH_FFT::random( void )
{

  boost::mt19937 randomNumbergenerator( time( 0 ) );

  const FT limit = 0.1 ;

  typedef boost::random::uniform_real_distribution< FT > uniform;

  uniform distribution( -limit , limit );

  boost::variate_generator< boost::mt19937&, uniform >  gen( randomNumbergenerator, distribution );

  FT mean=0;

  int NN=0;
  
  for(unsigned int i=0; i < nx; ++i)
    for(unsigned int j=0; j < ny; ++j) {

      FT rand = gen();

      
      fr(i,j).re = rand;
      fr(i,j).im = 0;

      mean += rand;
      ++NN;
    }

  mean /= NN;

  for(unsigned int i=0; i < nx; ++i)
    for(unsigned int j=0; j < ny; ++j)
       fr(i,j).re -= mean;


  freq_scramble( fr );

  Forward2_f.fft( fr , f );
  
  Forward2_f.Normalize( f );

  freq_scramble( fr );

  return;

}

void CH_FFT::init(void)
{
  for(uint i=0; i < nx; ++i)
    //    q_y( i ) = freq_pref * ( nx / 2.0  - FT(i) ) ;
    q_y( i ) = freq_pref * (  FT(i) - nx /2.0 ) ;

  for(uint j=0; j < ny; ++j)
    q_x( j ) = freq_pref *  ( FT(j)  - ny /2.0) ;


  for(uint i=0; i < nx; ++i)
    for(uint j=0; j < ny; ++j)
      //      FT q_x = freq_pref *  ( FT(j)  - ny /2.0) ;
      q2(i,j) = q_x(j) * q_x(j) + q_y(i) * q_y(i) ;


  
  
}  
  
  


void CH_FFT::mu_fields(void)
{

  // Backward2_f.fft( f , fr );
  // freq_scramble( fr );

  cube();

  freq_scramble( f3r );

  Forward2_f3.fft( f3r , f3 );
 
  Forward2_f3.Normalize( f3 );

  for(uint i=0; i < nx; i++)
    for(uint j=0; j < ny; j++) {
      const Complex imI(0,1);
	
      mu(i,j) = ( q2(i,j) /2.0 - 1) * f(i,j) + f3(i,j) ;
      grad_mu_x(i,j) = imI * q_x(j) * mu(i,j);
      grad_mu_y(i,j) = imI * q_y(i) * mu(i,j);

    }

  Backward2_mu.fft( mu , mu_r );
  freq_scramble( mu_r );

  Backward2_gmx.fft( grad_mu_x , grad_mu_x_r );
  freq_scramble( grad_mu_x_r );

  Backward2_gmy.fft( grad_mu_y , grad_mu_y_r );
  freq_scramble( grad_mu_y_r );


}

void CH_FFT::force_fields(void)
{

    // funny way, but a*b is not supported in Complex.h
  for(uint i=0; i < nx; i++)
    for(uint j=0; j < ny; j++) {
      force_x_r(i,j) =  -fr(i,j) * grad_mu_x_r(i,j);
      force_y_r(i,j) =  -fr(i,j) * grad_mu_y_r(i,j);

	//debug
	//force_x_r(i,j) = ( ( (i == nx / 2) && (j == ny / 2) ) ? nx * ny : 0) ;
	//force_y_r(i,j) = force_x_r(i,j) ;

    }

  freq_scramble( force_x_r );
  Forward2_forcex.fft( force_x_r , force_x );
  Forward2_forcex.Normalize( force_x );
  freq_scramble( force_x_r );
  
  freq_scramble( force_y_r );
  Forward2_forcey.fft( force_y_r , force_y );
  Forward2_forcey.Normalize( force_y );
  freq_scramble( force_y_r );
}

void CH_FFT::vel_fields(void)
{
  for(uint i=0; i < nx; i++)
    for(uint j=0; j < ny; j++) {
      FT qq2= q2(i,j);
      FT den = qq2 + sqrt( qq2 ) / L_SD ;

      if( den < 1e-16) {
	v_x(i,j) = 0;
	v_y(i,j) = 0;
	continue;
      }

      // readability.-
      FT qqx= q_x(j);
      FT qqy= q_y(i);

      Complex fx= force_x(i,j);
      Complex fy= force_y(i,j);

      Complex q_dot_F=  qqx * fx + qqy * fy ;

      //	Complex q_dot_F_over_q2=  q_dot_F / qq2;
	
      v_x(i,j) =  (fx - qqx * q_dot_F / qq2 ) / den ;
      v_y(i,j) =  (fy - qqy * q_dot_F / qq2 ) / den ;
    }

  Backward2_vx.fft( v_x , v_x_r );
  freq_scramble( v_x_r );

  
  Backward2_vy.fft( v_y , v_y_r );
  freq_scramble( v_y_r );

  return;
}


void CH_FFT::all_fields(void) {
  mu_fields();
  force_fields();
  vel_fields();

  return;
}


void CH_FFT::evolve(const FT b) {

  for(uint i=0; i < nx; i++)
    for(uint j=0; j < ny; j++) {
	
      Complex rhs = f(i,j) - b * ( q2(i,j) * f3(i,j) ) ;
      Complex lhs = 1  + b * q2(i,j) * ( q2(i,j) / 2.0 - 1.0 ) ;
	
      f(i,j) = rhs / lhs ;
      
    }

  Backward2_f.fft( f , fr );
  freq_scramble( fr );

}

void CH_FFT::histogram(const std::string& name,
		       const int time, const c_array& ff ) {

  histo.clear();
  //  h_count.clear();

  //  c_array gg = ff;
  //  freq_scramble( gg );
  
  for(uint i=0; i < nx; i++)
    for(uint j=0; j < ny; j++) {
      // FT qq=sqrt( q2(i,j) );

      FT re=real( ff(i,j) );
      FT im=imag( ff(i,j) );

      int ii = i - nx /2.0  ;

      int jj = j - ny /2.0 ;

      long int iijj= ii * ii + jj * jj;

      //      hist_data data;

      histo[iijj].f += re*re + im*im ;
      ++histo[iijj].count;

      if(histo[iijj].count == 1) { // 1st time

	FT qqx= q_x(j);
	FT qqy= q_y(i);

	histo[iijj].q2 =  qqx*qqx + qqy*qqy;
      }
      
      
      //      ++h_count[iijj];
//      cout << i << " " << j << "  " << iijj  // string (key)
//	   << "  "	   <<       histo[iijj].q2
//	   << "  "	   <<       histo[iijj].count
//	   << "  "	   <<       histo[iijj].f
//	   << "  "	   <<       re*re + im*im 
//	   << endl ;

    }

  std::map<long int,hist_data>::iterator it;

  //  std::map<int, double>::iterator it;
  //  std::map<int, int>::iterator c_it = h_count.begin();

  //   // average:
  // for ( it = histo.begin(); it != histo.end(); it++ ) {
  //   //    hist_data data = it->second ;
  //   it->second.f /= FT( it->second.count );    //FT(c_it-> second);
  //   //    c_it++;
  // }

  std::stringstream  namefile;
  namefile << time << "/histo_" << name << "_" << time << ".dat";

  //  cout << "writing on file : " << namefile.str() << endl;

  std::ofstream main_data;

  main_data.open(namefile.str().c_str() );

  for ( it = histo.begin(); it != histo.end(); it++ ) {
    //    FT qq2 = freq_pref*sqrt(FT (  it->first ) );

    FT qq = sqrt(it->second.q2 );
    FT qf= qq * it->second.f / FT(  it->second.count );
    
    main_data << qq
              << "  "
              << qf
	      // << "  "
              // << it->second.count
              << std::endl ;
  }
  
  main_data.close();

  return;


}