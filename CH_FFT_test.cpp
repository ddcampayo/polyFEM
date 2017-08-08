// Compilation:
// g++ -I ./fftwpp -fopenmp CH_FFT_test.cpp CH_FFT.o fftwpp/fftw++.cc -lfftw3 -lfftw3_omp

#include "CH_FFT.h"

using std::cout;
using std::endl;


int main(int argc, char* argv[])
{

  FT L=128;
  unsigned int nx=512;
  unsigned int ny=nx;


  CH_FFT CH( L, nx );
  
  CH.random( );
  //  CH.all_fields();

  //      q4[i][j] = q2[i][j] * q2[i][j];

  FT Co_D=1;

  FT Dx = L/FT(nx) ;
  FT Dx2 = Dx * Dx;

  FT D=1;
  
  FT Dt = Dx2 * Dx2 /D * Co_D ;

  cout << "# Time step for Co " << Co_D << " :  "  << Dt << endl;

  // Dt=0.01;

  // cout << "# Time step set to " << Dt << endl;
  
  FT b = D*Dt;

  FT t;

  int n_steps = 100000;
  int draw_every=100;
  
  // c_array one_plus_b_q2_q4(nx,ny,align);

  // for(uint i=0; i < nx; i++)
  //   for(uint j=0; j < ny; j++)
  //     one_plus_b_q2_q4[i][j] = 1  + b * q2[i][j] * ( q2[i][j] / 2.0 - 1.0 )   ;

  for(int step=0; step< n_steps ; ++step) {

    t= step*Dt;

    CH.all_fields();

    CH.evolve( b );

    if( step % draw_every == 0){
      cout << "Drawing on step " << step << endl;

      CH.draw( "phi", step, CH.field_f() );

      CH.draw( "mu", step, CH.field_mu() );

      CH.draw( "grad_mu_x", step, CH.field_grad_mu_x() );

      CH.draw( "grad_mu_y", step, CH.field_grad_mu_y() );

      CH.draw( "force_x", step, CH.field_force_x() );

      CH.draw( "force_y", step, CH.field_force_y() );

      CH.draw( "vel_x", step, CH.field_vel_x() );

      CH.draw( "vel_y", step, CH.field_vel_y() );

      CH.histogram( "phi",  step, CH.field_fq() );
    }
  }

  cout << "#Run " << n_steps << " , final time  " << t << endl;
  
  //  Backward2_f.fft( f );
  
  // freq_scramble( f );

  //  printout( "f" , f );

  return 0;
}



