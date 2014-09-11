#include <iostream>
#include <fstream>
#include <complex>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/resource.h>
#include <omp.h>

using namespace std;

#define PI 3.141592653589793
#define QE 1.602177e-19
#define E0 8.85419e-12
#define U0 1.2566e-6
#define ME 9.10939e-31
#define KB 1.38066e-23
#define C 3.0000e8 
//2.9979246e8

///////////////////////////////////////////////////////////////////////////////////
// 
// Meteor FDTD model
//
// version 1.0: January 2014: first implementation in 3D
// January 20, 2014: added -x face PML; added J equations
// version 1.1: replacing PML in incident field with one-way wave equation (1st order Mur)
// version 1.2: adding TF/SF in both polarizations; reading meteor plasma from file (2D)
// and extrapolating in z dimension with gaussian size
// Feb 13, 2014: adding different source types, polarizations
// version 1.3: Feb 15, 2014: changing incident source to 1D
// version 2.0: July 10, 2014: includes near-to-far field transformation, 
// and ability to use a scattering "plate" rather than a meteor to calibrate RCS. 
// Results show that it matches theory within 5-10%.
// version 2.1: September 3, 2014: added more inputs from Matlab, and generate better log file.
// "doplate" can be 1 or 2; 1 for plate, 2 for sphere. Updated sphere to use a real 
// conductor (s = 1e7) rather than PEC walls.
// fork: working on adding magnetic field back in! flag to turn it off and use Young DI method.
//
///////////////////////////////////////////////////////////////////////////////////

int main()
{

  // time the whole ting
  double tStart = omp_get_wtime();

  //------------------------------------------------------------------------//
  //-------READ INPUTS FROM FILE -------------------------------------------//

  FILE * inputFile;
  inputFile = fopen("inputs.dat","rb");

  // inputs:
  // source wavelength or frequency
  // source modulation (gaussian params?)
  // meteor size (min and max of elongated shape)
  // meteor maximum density
  // meteor location (x,y,z coordinates of peak density)
  // meteor rotation angle
  // meteor type (gaussian or other from sigrid)
  // nntf: number of NTF frequencies; fntf: frequencies
  // doplate: flag to turn on 1 m^2 plate scatterer, to test NTF
  // doB0: flag to use magnetic field and L&K method, versus simple DI method

  int dopml;
  int savefields;
  double f0, fsig;
  int sourcetype, sourcepol;
  double gridfactor;
  double dt;
  int tsteps;
  int Nx, Ny, Nz;
  double metsize [2];
  double wpmax, metangle;
  double numax;
  int mettype;
  double metloc [3];
  int nprobes;
  int numfiles;
  int nntf;
  int doplate;
  double platesize;
  int doB0;

  fread(&dopml,sizeof(int),1,inputFile);
  fread(&savefields,sizeof(int),1,inputFile);
  fread(&f0,sizeof(double),1,inputFile);
  fread(&fsig,sizeof(double),1,inputFile);
  fread(&sourcetype,sizeof(int),1,inputFile);
  fread(&sourcepol,sizeof(int),1,inputFile);
  fread(&gridfactor,sizeof(double),1,inputFile);
  fread(&dt,sizeof(double),1,inputFile);
  fread(&tsteps,sizeof(int),1,inputFile);
  fread(&Nx,sizeof(int),1,inputFile);
  fread(&Ny,sizeof(int),1,inputFile);
  fread(&Nz,sizeof(int),1,inputFile);
  fread(&metsize,sizeof(double),2,inputFile);
  fread(&wpmax,sizeof(double),1,inputFile);
  fread(&numax,sizeof(double),1,inputFile);
  fread(&metloc,sizeof(double),3,inputFile);
  fread(&metangle,sizeof(double),1,inputFile);
  fread(&mettype,sizeof(int),1,inputFile);
  fread(&numfiles,sizeof(int),1,inputFile);
  fread(&nprobes,sizeof(int),1,inputFile);
  int probex[nprobes], probey[nprobes], probez[nprobes];
  fread(&probex,sizeof(int),nprobes,inputFile);
  fread(&probey,sizeof(int),nprobes,inputFile);
  fread(&probez,sizeof(int),nprobes,inputFile);
  fread(&nntf,sizeof(int),1,inputFile);
  double fntf [nntf];
  fread(&fntf,sizeof(double),nntf,inputFile);
  fread(&doplate,sizeof(int),1,inputFile);
  fread(&platesize,sizeof(double),1,inputFile);
  fread(&doB0,sizeof(int),1,inputFile);

  fclose(inputFile);

  double lambda = C/f0;
     
  // turn on/off TF/SF (for testing)
  int dotfsf = 1;

  //------------------------------------------------------------------------//
  //-------GRID PARAMETERS -------------------------------------------------//

  double dx, dy, dz;
    
  dx = lambda/gridfactor;
  dy = lambda/gridfactor;
  dz = lambda/gridfactor;
    
  // define based on plasma frequency, source frequency, and meteor size

  double x[Nx], y[Ny], z[Nz];

  for (int i = 0; i < Nx; i++) {
    x[i] = -(Nx-1)/2*dx + dx*i;
    y[i] = -(Ny-1)/2*dy + dy*i;
    z[i] = -(Nz-1)/2*dz + dz*i;
  }

  // fractional time to write outputs to files
  double partialtime;

  // source location = slice locations

  int ssx = (Nx+1)/2;
  int ssy = (Ny+1)/2;
  int ssz = (Nz+1)/2;

  //------------------------------------------------------------------------//
  //-------DEFINE SOURCE----------------------------------------------------//

  double EsourceA [tsteps];
  double EsourceB [tsteps];
  double rdelayA = 1.6*sqrt(6)/PI/f0;
  double rdelayB = (1.6+1/sqrt(3)) * sqrt(6)/PI/f0;

  for (int t = 0; t < tsteps; t++) {
    if (sourcetype == 1) {
      // sinusoid
      EsourceA[t] = sin(2*PI*f0*t*dt);
      EsourceB[t] = cos(2*PI*f0*t*dt);
    } else if (sourcetype == 2) {
      // modulated gaussian
      EsourceA[t] = sin(2*PI*f0*(t-3*fsig)*dt) * exp(-(t-3*fsig)*(t-3*fsig)/(2*fsig*fsig));
      EsourceB[t] = cos(2*PI*f0*(t-3*fsig)*dt) * exp(-(t-3*fsig)*(t-3*fsig)/(2*fsig*fsig));
    } else if (sourcetype == 3) {
      // ricker wavelet
      EsourceA[t] = (1 - 2*PI*PI*f0*f0*(t*dt-rdelayA)*(t*dt-rdelayA)) * exp(-PI*PI*f0*f0*(t*dt-rdelayA)*(t*dt-rdelayA));
      EsourceB[t] = (1 - 2*PI*PI*f0*f0*(t*dt-rdelayB)*(t*dt-rdelayB)) * exp(-PI*PI*f0*f0*(t*dt-rdelayB)*(t*dt-rdelayB));
    } else {
      // default back to sinusoid
      EsourceA[t] = sin(2*PI*f0*t*dt);
      EsourceB[t] = cos(2*PI*f0*t*dt);
     }
  }

  //------------------------------------------------------------------------//
  //-------MULTIPLY FACTORS FOR FDTD ---------------------------------------//

  double sigm = 0;
  double c1h = (2*U0 - sigm*dt)/(2*U0 + sigm*dt);
  double c2h = 2*dt/(2*U0 + sigm*dt);

  double sig = 0;
  double c1e0, c2e0;
  double c1e [Nx][Ny][Nz];
  double c2e [Nx][Ny][Nz];

  double dist;

  c1e0 = (2*E0 - sig*dt)/(2*E0 + sig*dt);
  c2e0 = 2*dt/(2*E0 + sig*dt);

  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nz; k++) {
	c1e[i][j][k] = c1e0;
	c2e[i][j][k] = c2e0;
	
	if (doplate == 2) { // spherical scatterer with high eps, sigma
	  dist = sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k]);
	  if (dist < platesize) {
	    c1e[i][j][k] = (2*E0 - 1e7*dt)/(2*E0 + 1e7*dt);
	    c2e[i][j][k] = 2*dt/(2*E0 + 1e7*dt);
	  }
	}

      }
    }
  }

  //------------------------------------------------------------------------//
  //-------MAGNETIC FIELD --------------------------------------------------//

  // read 0D magnetic field (varying B doesn't buy much in 3D)
  double Bx, By, Bz;
  FILE * BFile;
  BFile = fopen("B0.dat","rb");
  fread(&Bx,sizeof(double),1,BFile);
  fread(&By,sizeof(double),1,BFile);
  fread(&Bz,sizeof(double),1,BFile);
  fclose(BFile);

  // compute gyrofrequencies, store in 2D arrays

  double wcex, wcey, wcez, wce0;

  wcex = -QE * Bx / ME;
  wcey = -QE * By / ME;
  wcez = -QE * Bz / ME;
  wce0 = sqrt(wcex*wcex + wcey*wcey + wcez*wcez);

 // these are A and K matrices for Lee&Kalluri current solution
  double Ae [3][3];
  double Ke [3][3];
  double Emid [3];

  // intermediate calculations in Lee&Kalluri
  double Ee1, Ee2, Se1, Ce1, Ce2, Ce3, Ce4, Jex0, Jey0;

  //------------------------------------------------------------------------//
  //------- FIELD INITIALIZATION --------------------------------------------//

  // electric field vector
  double Ex [Nx][Ny][Nz];
  double Ey [Nx][Ny][Nz];
  double Ez [Nx][Ny][Nz];
  // magnetic field vector
  double Hx [Nx][Ny][Nz];
  double Hy [Nx][Ny][Nz];
  double Hz [Nx][Ny][Nz];
  // spatial-averaged fields
  double Exm, Eym, Ezm, Hxm, Hym, Hzm, Sx, Sy, Sz;

  // incident fields for TF/SF: only need y/z components
  double Eyinc [Nx];
  double Ezinc [Nx];
  double Hyinc [Nx];
  double Hzinc [Nx];

  // probe fields
  double Exprobe [tsteps][nprobes];
  double Eyprobe [tsteps][nprobes];
  double Ezprobe [tsteps][nprobes];
  double Hxprobe [tsteps][nprobes];
  double Hyprobe [tsteps][nprobes];
  double Hzprobe [tsteps][nprobes];
  int ix, iy, iz;  // probe point indices

  // electron current
  double Jex [Nx][Ny][Nz];
  double Jey [Nx][Ny][Nz];
  double Jez [Nx][Ny][Nz];

  double S [Nx][Ny][Nz];

  double Emag [Nx][Ny][Nz];

  // slices for output to file
  double Eslicez [2][Nx][Ny];
  double Eslicey [2][Nx][Nz];
  double Eslicex [2][Ny][Nz];
  double Jslicez [2][Nx][Ny];
  double Jslicey [2][Nx][Nz];
  double Jslicex [2][Ny][Nz];
  double Hslicez [2][Nx][Ny];
  double Hslicey [2][Nx][Nz];
  double Hslicex [2][Ny][Nz];

  double Emagslice [Nx][Ny];

  // initialize 

  memset(Ex,0,sizeof(Ex));
  memset(Ey,0,sizeof(Ey));
  memset(Ez,0,sizeof(Ez));
  memset(Hx,0,sizeof(Hx));
  memset(Hy,0,sizeof(Hy));
  memset(Hz,0,sizeof(Hz));
  memset(Jex,0,sizeof(Jex));
  memset(Jey,0,sizeof(Jey));
  memset(Jez,0,sizeof(Jez));

  memset(Eyinc,0,sizeof(Eyinc));
  memset(Ezinc,0,sizeof(Ezinc));
  memset(Hyinc,0,sizeof(Hyinc));
  memset(Hzinc,0,sizeof(Hzinc));

  memset(Emag,0,sizeof(Emag));

  memset(Eslicez,0,sizeof(Eslicez));
  memset(Eslicey,0,sizeof(Eslicey));
  memset(Eslicex,0,sizeof(Eslicex));
  memset(Jslicez,0,sizeof(Jslicez));
  memset(Jslicey,0,sizeof(Jslicey));
  memset(Jslicex,0,sizeof(Jslicex));
  memset(Hslicez,0,sizeof(Hslicez));
  memset(Hslicey,0,sizeof(Hslicey));
  memset(Hslicex,0,sizeof(Hslicex));

  memset(Emagslice,0,sizeof(Emagslice));


  //------------------------------------------------------------------------//
  //-------PLASMA DENSITY and COLL FREQUENCY -------------------------------//

  // set up plasma density and collision frequency

  double wpe [Nx][Ny][Nz];     // plasma frequency
  double wp2 [Nx][Ny];         // 2D slice of wp to be read from file
  double nue [Nx][Ny][Nz];     // collision frequency
  double nu2 [Nx][Ny];         // 2D slice of nu to be read from file

  // read wp2 from file

  FILE * wpFile;
  wpFile = fopen("wp.dat","rb");
  fread(&wp2,sizeof(double),Nx*Ny,wpFile);
  fclose(wpFile);

  FILE * nuFile;
  nuFile = fopen("nu.dat","rb");
  fread(&nu2,sizeof(double),Nx*Ny,nuFile);
  fclose(nuFile);

  double msig0 = metsize[0];

  // extrapolate 3D arrays
  if (!doplate) {
    for (int i = 0; i < Nx; i++) {
      for (int j = 0; j < Ny; j++) {
	for (int k = 0; k < Nz; k++) {
	  wpe[i][j][k] = wp2[i][j] * exp(-pow((z[k]-metloc[2]),2)/(msig0*msig0));
	  nue[i][j][k] = nu2[i][j] * exp(-pow((z[k]-metloc[2]),2)/(msig0*msig0));
	}    
      }
    }
  } else {
    // if plate scatterer, set wpe to zero everywhere
    memset(wpe,0,sizeof(wpe));
    memset(nue,0,sizeof(nue));
  }


  //------------------------------------------------------------------------//
  //------- SET UP PML -----------------------------------------------------//

  // initialize A, B, and P arrays even if we don't use them

  double pmlm = 4;
  int pmllen = 10;

  double A_exy [2*pmllen];
  double A_exz [2*pmllen];
  double A_eyx [2*pmllen];
  double A_eyz [2*pmllen];
  double A_ezx [2*pmllen];
  double A_ezy [2*pmllen];
  double A_hxy [2*pmllen];
  double A_hxz [2*pmllen];
  double A_hyx [2*pmllen];
  double A_hyz [2*pmllen];
  double A_hzx [2*pmllen];
  double A_hzy [2*pmllen];

  double B_exy [2*pmllen];
  double B_exz [2*pmllen];
  double B_eyx [2*pmllen];
  double B_eyz [2*pmllen];
  double B_ezx [2*pmllen];
  double B_ezy [2*pmllen];
  double B_hxy [2*pmllen];
  double B_hxz [2*pmllen];
  double B_hyx [2*pmllen];
  double B_hyz [2*pmllen];
  double B_hzx [2*pmllen];
  double B_hzy [2*pmllen];

  double P_exy [Nx][2*pmllen][Nz];
  double P_exz [Nx][Ny][2*pmllen];
  double P_eyx [2*pmllen][Ny][Nz];
  double P_eyz [Nx][Ny][2*pmllen];
  double P_ezx [2*pmllen][Ny][Nz];
  double P_ezy [Nx][2*pmllen][Nz];

  double P_hxy [Nx][2*pmllen][Nz];
  double P_hxz [Nx][Ny][2*pmllen];
  double P_hyx [2*pmllen][Ny][Nz];
  double P_hyz [Nx][Ny][2*pmllen];
  double P_hzx [2*pmllen][Ny][Nz];
  double P_hzy [Nx][2*pmllen][Nz];

  // initialize
  memset(P_exz,0,sizeof(P_exz));
  memset(P_eyz,0,sizeof(P_eyz));
  memset(P_hxz,0,sizeof(P_hxz));
  memset(P_hyz,0,sizeof(P_hyz));
  memset(P_exy,0,sizeof(P_exy));
  memset(P_ezy,0,sizeof(P_ezy));
  memset(P_hxy,0,sizeof(P_hxy));
  memset(P_hzy,0,sizeof(P_hzy));
  memset(P_eyx,0,sizeof(P_eyx));
  memset(P_ezx,0,sizeof(P_ezx));
  memset(P_hyx,0,sizeof(P_hyx));
  memset(P_hzx,0,sizeof(P_hzx));

  // Psi fields for incident wave at +x boundary. Can use same A, B since they don't update
  double P_eyxI [2*pmllen];
  double P_ezxI [2*pmllen];
  double P_hyxI [2*pmllen];
  double P_hzxI [2*pmllen];

  memset(P_eyxI,0,sizeof(P_eyxI));
  memset(P_ezxI,0,sizeof(P_ezxI));
  memset(P_hyxI,0,sizeof(P_hyxI));
  memset(P_hzxI,0,sizeof(P_hzxI));

  // end initialization

  if (dopml) {
    double simax = 1.5 * (double)(pmlm + 1) / (150 * PI * dx);
    double kamax = 1;
    double almax = 0;

    double sx [pmllen];
    double sxm [pmllen];
    double sy [pmllen];
    double sym [pmllen];
    double sz [pmllen];
    double szm [pmllen];

    double kx [pmllen];
    double kxm [pmllen];
    double ky [pmllen];
    double kym [pmllen];
    double kz [pmllen];
    double kzm [pmllen];

    double ax [pmllen];
    double axm [pmllen];
    double ay [pmllen];
    double aym [pmllen];
    double az [pmllen];
    double azm [pmllen];
    double mf;

    for (int m = 0; m < pmllen; m++) {
      // note shifts here compared to matlab, due to zero index
      mf = (double)m;
      sx[m] = simax * pow(((mf+0.5)/pmllen),pmlm);
      sxm[m] = simax * pow(((mf+1)/pmllen),pmlm);
      sy[m] = simax * pow(((mf+0.5)/pmllen),pmlm);
      sym[m] = simax * pow(((mf+1)/pmllen),pmlm);
      sz[m] = simax * pow(((mf+0.5)/pmllen),pmlm);
      szm[m] = simax * pow(((mf+1)/pmllen),pmlm);

      kx[m] = 1 + (kamax-1) * pow(((mf+0.5)/pmllen),pmlm);
      kxm[m] = 1 + (kamax-1) * pow(((mf+1)/pmllen),pmlm);
      ky[m] = 1 + (kamax-1) * pow(((mf+0.5)/pmllen),pmlm);
      kym[m] = 1 + (kamax-1) * pow(((mf+1)/pmllen),pmlm);
      kz[m] = 1 + (kamax-1) * pow(((mf+0.5)/pmllen),pmlm);
      kzm[m] = 1 + (kamax-1) * pow(((mf+1)/pmllen),pmlm);

      ax[m] = almax * pow(((pmllen-mf-0.5)/pmllen),pmlm);
      axm[m] = almax * pow(((pmllen-mf-1)/pmllen),pmlm);
      ay[m] = almax * pow(((pmllen-mf-0.5)/pmllen),pmlm);
      aym[m] = almax * pow(((pmllen-mf-1)/pmllen),pmlm);
      az[m] = almax * pow(((pmllen-mf-0.5)/pmllen),pmlm);
      azm[m] = almax * pow(((pmllen-mf-1)/pmllen),pmlm);

      // do I need z components too?
    }

    // update / initialize A and B vectors
    for (int m = pmllen; m < 2*pmllen; m++) {
      int n = m - pmllen;
      B_exy[m] = exp(-((sy[n]/ky[n]) + ay[n])*dt/E0);
      A_exy[m] = sy[n] / (sy[n]*ky[n] + pow(ky[n],2)*ay[n]) * (B_exy[m] - 1);
      B_eyx[m] = exp(-((sx[n]/kx[n]) + ax[n])*dt/E0);
      A_eyx[m] = sx[n] / (sx[n]*kx[n] + pow(kx[n],2)*ax[n]) * (B_eyx[m] - 1);
      B_ezy[m] = exp(-((sy[n]/ky[n]) + ay[n])*dt/E0);
      A_ezy[m] = sy[n] / (sy[n]*ky[n] + pow(ky[n],2)*ay[n]) * (B_ezy[m] - 1);
      B_ezx[m] = exp(-((sx[n]/kx[n]) + ax[n])*dt/E0);
      A_ezx[m] = sx[n] / (sx[n]*kx[n] + pow(kx[n],2)*ax[n]) * (B_ezx[m] - 1);

      B_hxy[m] = exp(-((sym[n]/kym[n]) + aym[n])*dt/E0);
      A_hxy[m] = sym[n] / (sym[n]*kym[n] + pow(kym[n],2)*aym[n]) * (B_hxy[m] - 1);
      B_hyx[m] = exp(-((sxm[n]/kxm[n]) + axm[n])*dt/E0);
      A_hyx[m] = sxm[n] / (sxm[n]*kxm[n] + pow(kxm[n],2)*axm[n]) * (B_hyx[m] - 1);
      B_hzy[m] = exp(-((sym[n]/kym[n]) + aym[n])*dt/E0);
      A_hzy[m] = sym[n] / (sym[n]*kym[n] + pow(kym[n],2)*aym[n]) * (B_hzy[m] - 1);
      B_hzx[m] = exp(-((sxm[n]/kxm[n]) + axm[n])*dt/E0);
      A_hzx[m] = sxm[n] / (sxm[n]*kxm[n] + pow(kxm[n],2)*axm[n]) * (B_hzx[m] - 1);  

      // using theta conductivities for phi direction. 

      B_exz[m] = exp(-((sz[n]/kz[n]) + az[n])*dt/E0);
      A_exz[m] = sz[n] / (sz[n]*kz[n] + pow(kz[n],2)*az[n]) * (B_exz[m] - 1);
      B_eyz[m] = exp(-((sz[n]/kz[n]) + az[n])*dt/E0);
      A_eyz[m] = sz[n] / (sz[n]*kz[n] + pow(kz[n],2)*az[n]) * (B_eyz[m] - 1);
      B_hxz[m] = exp(-((szm[n]/kzm[n]) + azm[n])*dt/E0);
      A_hxz[m] = szm[n] / (szm[n]*kzm[n] + pow(kzm[n],2)*azm[n]) * (B_hxz[m] - 1);
      B_hyz[m] = exp(-((szm[n]/kzm[n]) + azm[n])*dt/E0);
      A_hyz[m] = szm[n] / (szm[n]*kzm[n] + pow(kzm[n],2)*azm[n]) * (B_hyz[m] - 1);
    }
    for (int m = 0; m < pmllen; m++) {
      B_exy[m] = B_exy[2*pmllen-1-m];
      A_exy[m] = A_exy[2*pmllen-1-m];
      B_eyx[m] = B_eyx[2*pmllen-1-m];
      A_eyx[m] = A_eyx[2*pmllen-1-m];
      B_ezy[m] = B_ezy[2*pmllen-1-m];
      A_ezy[m] = A_ezy[2*pmllen-1-m];
      B_ezx[m] = B_ezx[2*pmllen-1-m];
      A_ezx[m] = A_ezx[2*pmllen-1-m];

      B_hxy[m] = B_hxy[2*pmllen-1-m];
      A_hxy[m] = A_hxy[2*pmllen-1-m];
      B_hyx[m] = B_hyx[2*pmllen-1-m];
      A_hyx[m] = A_hyx[2*pmllen-1-m];
      B_hzy[m] = B_hzy[2*pmllen-1-m];
      A_hzy[m] = A_hzy[2*pmllen-1-m];
      B_hzx[m] = B_hzx[2*pmllen-1-m];
      A_hzx[m] = A_hzx[2*pmllen-1-m];

      B_exz[m] = B_exz[2*pmllen-1-m];
      A_exz[m] = A_exz[2*pmllen-1-m];
      B_eyz[m] = B_eyz[2*pmllen-1-m];
      A_eyz[m] = A_eyz[2*pmllen-1-m];
      B_hxz[m] = B_hxz[2*pmllen-1-m];
      A_hxz[m] = A_hxz[2*pmllen-1-m];
      B_hyz[m] = B_hyz[2*pmllen-1-m];
      A_hyz[m] = A_hyz[2*pmllen-1-m];
    }

    // PML: everything up to here seems gtg


  } // if dopml

  // offsets for PML
  int xshift = Nx-2*pmllen-1;
  int yshift = Ny-2*pmllen-1;
  int zshift = Nz-2*pmllen-1;

  // ----------------------------------
  // TF/SF indices
  int bx1 = pmllen + 10;
  int bx2 = Nx - pmllen - 10;
  int by1 = pmllen + 10;
  int by2 = Ny - pmllen - 10;
  int bz1 = pmllen + 10;
  int bz2 = Nz - pmllen - 10;

  // ----------------------------------
  // NTF indices and equivalent currents
  // "front" and "back" are in x-direction, direction of source propagation;
  // "top" and "bottom" are y-direction; "left" and right" are z-direction
  // (imagine looking down the direction of the radar beam)

  int ntfb = pmllen + 5;

  // DFTs of source functions
  double Ezdft [2*nntf];
  double Eydft [2*nntf];

  // I'm going to make these arrays cover the whole wall, then carve out the useful stuff in Matlab later
  // front and back: front is -xhat
  double My_front [Ny][Nz][2*nntf]; // M = -n x E, so xhat x Ez -> -My
  double Mz_front [Ny][Nz][2*nntf]; // xhat x Ey -> +Mz
  double Jy_front [Ny][Nz][2*nntf]; // J = n x H, so -xhat x Hz -> +Jy
  double Jz_front [Ny][Nz][2*nntf]; // -xhat x Hy -> -Jz
  double My_back [Ny][Nz][2*nntf];  // -xhat x Ez -> +My
  double Mz_back [Ny][Nz][2*nntf];  // -xhat x Ey -> -Mz
  double Jy_back [Ny][Nz][2*nntf];  // xhat x Hz -> -Jy
  double Jz_back [Ny][Nz][2*nntf];  // xhat x Hy -> +Jz

  // top and bottom: top is +yhat
  double Mx_top [Nx][Nz][2*nntf]; // M = -n x E, so -yhat x Ez -> -Mx
  double Mz_top [Nx][Nz][2*nntf]; // -yhat x Ex -> +Mz
  double Jx_top [Nx][Nz][2*nntf]; // J = n x H, so yhat x Hz -> +Jx
  double Jz_top [Nx][Nz][2*nntf]; // yhat x Hx -> -Jz
  double Mx_bot [Nx][Nz][2*nntf]; // yhat x Ez -> +Mx
  double Mz_bot [Nx][Nz][2*nntf]; // yhat x Ex -> -Mz
  double Jx_bot [Nx][Nz][2*nntf]; // -yhat x Hz -> -Jx
  double Jz_bot [Nx][Nz][2*nntf]; // -yhat x Hx -> +Jz

  // left and right: left is -zhat
  double Mx_left [Nx][Ny][2*nntf];  // M = -n x E, so zhat x Ey -> -Mx
  double My_left [Nx][Ny][2*nntf];  // zhat x Ex -> +My
  double Jx_left [Nx][Ny][2*nntf];  // J = n x H, so -zhat x Hy -> +Jx
  double Jy_left [Nx][Ny][2*nntf];  // -zhat x Hx -> -Jy
  double Mx_right [Nx][Ny][2*nntf]; // -zhat x Ey -> +Mx
  double My_right [Nx][Ny][2*nntf]; // -zhat x Ex -> -My
  double Jx_right [Nx][Ny][2*nntf]; // zhat x Hy -> -Jx
  double Jy_right [Nx][Ny][2*nntf]; // zhat x Hx -> +Jy

  memset(My_front,0,sizeof(My_front));
  memset(Mz_front,0,sizeof(Mz_front));
  memset(Jy_front,0,sizeof(Jy_front));
  memset(Jz_front,0,sizeof(Jz_front));
  memset(My_back,0,sizeof(My_back));
  memset(Mz_back,0,sizeof(Mz_back));
  memset(Jy_back,0,sizeof(Jy_back));
  memset(Jz_back,0,sizeof(Jz_back));

  memset(Mx_top,0,sizeof(Mx_top));
  memset(Mz_top,0,sizeof(Mz_top));
  memset(Jx_top,0,sizeof(Jx_top));
  memset(Jz_top,0,sizeof(Jz_top));
  memset(Mx_bot,0,sizeof(Mx_bot));
  memset(Mz_bot,0,sizeof(Mz_bot));
  memset(Jx_bot,0,sizeof(Mx_bot));
  memset(Jz_bot,0,sizeof(Mz_bot));

  memset(Mx_left,0,sizeof(Mx_left));
  memset(My_left,0,sizeof(My_left));
  memset(Jx_left,0,sizeof(Jx_left));
  memset(Jy_left,0,sizeof(Jy_left));
  memset(Mx_right,0,sizeof(Mx_right));
  memset(My_right,0,sizeof(My_right));
  memset(Jx_right,0,sizeof(Jx_right));
  memset(Jy_right,0,sizeof(Jy_right));


  //////////////////////////////////////////////////

  // measure memory usage, just for kicks                                                                    
  struct rusage ru;
  double memusage;
  getrusage(RUSAGE_SELF, &ru);
  memusage = (double)ru.ru_maxrss / 1024 / 1024;


  //------------------------------------------------------------------------//
  //------- OUTPUT TO LOG FILE ---------------------------------------------//

  ofstream logfile;
  logfile.open("log.txt");

  logfile << "------------------------------------------\n";
  logfile << "3D Meteor Scattering Simulation\n";
  logfile << "------------------------------------------\n";
  logfile << "dx = " << dx << "\n";
  logfile << "source type = " << sourcetype << "\n";
  logfile << "source polarization = " << sourcepol << "\n";
  logfile << "Time step = " << dt << "\n";
  logfile << "Total time steps = " << tsteps << "\n";
  logfile << "Total NTF frequencies = " << nntf << "\n";
  logfile << "last NTF frequency = " << fntf[nntf-1] << "\n";
  logfile << "Total memory usage = " << memusage << " GB\n";
  logfile.close();


  ////////////////////////////////////////////////////////////////////////////
  //                                                                      ////
  //------- MAIN TIME LOOP -----------------------------------------------////

  double tloopStart = omp_get_wtime();

  for (int t = 0; t < tsteps; t++) {

    // ----------------------------------------
    // Psi updates for H field
    if (dopml) {

      // Psi's applied to Hx
      for (int i = 0; i < Nx; i++) {
	for (int j = 0; j < pmllen; j++) {
	  for (int k = 0; k < Nz-1; k++) {
	    P_hxy[i][j][k] = B_hxy[j] * P_hxy[i][j][k] \
	      + A_hxy[j] * ( Ez[i][j+1][k] - Ez[i][j][k] ) / dy;
	  }
	}
	for (int j = Ny-pmllen-1; j < Ny-1; j++) {
	  for (int k = 0; k < Nz-1; k++) {
	    P_hxy[i][j-yshift][k] = B_hxy[j-yshift] * P_hxy[i][j-yshift][k] \
	      + A_hxy[j-yshift] * ( Ez[i][j+1][k] - Ez[i][j][k] ) / dy;
	  }
	}
      }

      for (int i = 0; i < Nx; i++) {
	for (int j = 0; j < Ny-1; j++) {
	  for (int k = 0; k < pmllen; k++) {
	    P_hxz[i][j][k] = B_hxz[k] * P_hxz[i][j][k] \
	      + A_hxz[k] * ( Ey[i][j][k+1] - Ey[i][j][k] ) / dz;
	  }
	  for (int k = Nz-pmllen-1; k < Nz-1; k++) {
	    P_hxz[i][j][k-zshift] = B_hxz[k-zshift] * P_hxz[i][j][k-zshift] \
	      + A_hxz[k-zshift] * ( Ey[i][j][k+1] - Ey[i][j][k] ) / dz;
	  }
	}
      }

      // Psi's applied to Hy
      for (int i = 0; i < pmllen; i++) {
	for (int j = 0; j < Ny; j++) {
	  for (int k = 0; k < Nz-1; k++) {
	    P_hyx[i][j][k] = B_hyx[i] * P_hyx[i][j][k] \
	      + A_hyx[i] * ( Ez[i+1][j][k] - Ez[i][j][k] ) / dx; 
	  }
	}
	/* P_hyxI[i] = B_hyx[i] * P_hyxI[i]		\
	  + A_hyx[i] * ( Ezinc[i+1] - Ezinc[i] ) / dx;  */
      }
      for (int i = Nx-pmllen-1; i < Nx-1; i++) {
	for (int j = 0; j < Ny; j++) {
	  for (int k = 0; k < Nz-1; k++) {
	    P_hyx[i-xshift][j][k] = B_hyx[i-xshift] * P_hyx[i-xshift][j][k] \
	      + A_hyx[i-xshift] * ( Ez[i+1][j][k] - Ez[i][j][k] ) / dx; 
	  }
	}
	P_hyxI[i-xshift] = B_hyx[i-xshift] * P_hyxI[i-xshift] \
	      + A_hyx[i-xshift] * ( Ezinc[i+1] - Ezinc[i] ) / dx; 
      }
 
      for (int i = 0; i < Nx-1; i++) {
	for (int j = 0; j < Ny; j++) {
	  for (int k = 0; k < pmllen; k++) {
	    P_hyz[i][j][k] = B_hyz[k] * P_hyz[i][j][k]			\
	      + A_hyz[k] * ( Ex[i][j][k+1] - Ex[i][j][k] ) / dz;
	  }
	  for (int k = Nz-pmllen-1; k < Nz-1; k++) {
	    P_hyz[i][j][k-zshift] = B_hyz[k-zshift] * P_hyz[i][j][k-zshift] \
	      + A_hyz[k-zshift] * ( Ex[i][j][k+1] - Ex[i][j][k] ) / dz;
	  }
	}
      }

      // Psi's applied to Hz
      for (int i = 0; i < pmllen; i++) {
	for (int j = 0; j < Ny-1; j++) {
	  for (int k = 0; k < Nz; k++) {
	    P_hzx[i][j][k] = B_hzx[i] * P_hzx[i][j][k] \
	      + A_hzx[i] * ( Ey[i+1][j][k] - Ey[i][j][k] ) / dx;
	  }
	}
	/*P_hzxI[i] = B_hzx[i] * P_hzxI[i]		\
	  + A_hzx[i] * ( Eyinc[i+1] - Eyinc[i] ) / dx; */
      }
      
      for (int i = Nx-pmllen-1; i < Nx-1; i++) {
	for (int j = 0; j < Ny-1; j++) {
	  for (int k = 0; k < Nz; k++) {
	    P_hzx[i-xshift][j][k] = B_hzx[i-xshift] * P_hzx[i-xshift][j][k] \
	      + A_hzx[i-xshift] * ( Ey[i+1][j][k] - Ey[i][j][k] ) / dx;
	  }
	}
	P_hzxI[i-xshift] = B_hzx[i-xshift] * P_hzxI[i-xshift] \
	  + A_hzx[i-xshift] * ( Eyinc[i+1] - Eyinc[i] ) / dx;
      }

      for (int i = 0; i < Nx-1; i++) {
	for (int j = 0; j < pmllen; j++) {
	  for (int k = 0; k < Nz; k++) {
	    P_hzy[i][j][k] = B_hzy[j] * P_hzy[i][j][k] \
	      + A_hzy[j] * ( Ex[i][j+1][k] - Ex[i][j][k] ) / dy;
	  }
	}
	for (int j = Ny-pmllen-1; j < Ny-1; j++) {
	  for (int k = 0; k < Nz; k++) {
	    P_hzy[i][j-yshift][k] = B_hzy[j-yshift] * P_hzy[i][j-yshift][k] \
	      + A_hzy[j-yshift] * (Ex[i][j+1][k] - Ex[i][j][k] ) / dy;
	  }
	}
      }

    } // if dopml


    // ----------------------------------------
    // Hx update

#pragma omp parallel for
    for (int i = 0; i < Nx; i++ ) {
      for (int j = 0; j < Ny-1; j++) {
	for (int k = 0; k < Nz-1; k++) {
	  Hx[i][j][k] = c1h * Hx[i][j][k] - c2h * \
	    ( ( Ez[i][j+1][k] - Ez[i][j][k] ) / dy - ( Ey[i][j][k+1] - Ey[i][j][k] ) / dz );
	}
      }
      // No Hxinc, since it is 1D!
    }

    // TFSF correction: Hx corrected on y, z walls of TFSF box
    if (dotfsf) {
    for (int i = bx1; i < bx2; i++) {
      for (int k = bz1; k < bz2; k++) {
	Hx[i][by1-1][k] = Hx[i][by1-1][k] + c2h * Ezinc[i] / dy;
	Hx[i][by2-1][k] = Hx[i][by2-1][k] - c2h * Ezinc[i] / dy;
      }
    }
    for (int i = bx1; i < bx2; i++) {
      for (int j = by1; j < by2; j++) {
	Hx[i][j][bz1-1] = Hx[i][j][bz1-1] - c2h * Eyinc[i] / dz;
	Hx[i][j][bz2-1] = Hx[i][j][bz2-1] + c2h * Eyinc[i] / dz;
      }
    }
    }

    // pml corrections
    if (dopml) {
      for (int i = 0; i < Nx; i++) {
	for (int j = 0; j < pmllen; j++) {
	  for (int k = 0; k < Nz-1; k++) {
	    Hx[i][j][k] = Hx[i][j][k] - c2h * P_hxy[i][j][k];
	  }
	}
	for (int j = Ny-pmllen-1; j < Ny-1; j++) {
	  for (int k = 0; k < Nz-1; k++) {
	    Hx[i][j][k] = Hx[i][j][k] - c2h * P_hxy[i][j-yshift][k];
	  }
	}
      }

      for (int i = 0; i < Nx; i++) {
	for (int j = 0; j < Ny-1; j++) {
	  for (int k = 0; k < pmllen; k++) {
	    Hx[i][j][k] = Hx[i][j][k] + c2h * P_hxz[i][j][k];
	  }
	  for (int k = Nz-pmllen-1; k < Nz-1; k++) {
	    Hx[i][j][k] = Hx[i][j][k] + c2h * P_hxz[i][j][k-zshift];
	  }
	}
      }

    } // if dopml

    // ----------------------------------------
    // Hy update
#pragma omp parallel for
    for (int i = 0; i < Nx-1; i++) {
      for (int j = 0; j < Ny; j++) {
	for (int k = 0; k < Nz-1; k++) {
	  // TF/SF
	  Hy[i][j][k] = c1h * Hy[i][j][k] - c2h * \
	    ( ( Ex[i][j][k+1] - Ex[i][j][k] ) / dz - ( Ez[i+1][j][k] - Ez[i][j][k] ) / dx );
	}
      }
      // incident field
      Hyinc[i] = c1h * Hyinc[i] + c2h * ( Ezinc[i+1] - Ezinc[i] ) / dx;
    }

    if (dotfsf) {
    // TFSF correction: Hy corrected on x, z walls of TFSF box, but there is no Exinc, so no z correction
    for (int j = by1; j < by2; j++) {
      for (int k = bz1; k < bz2; k++) {
	Hy[bx1-1][j][k] = Hy[bx1-1][j][k] - c2h * Ezinc[bx1] / dx;
	Hy[bx2-1][j][k] = Hy[bx2-1][j][k] + c2h * Ezinc[bx2-1] / dx; 
      }
    }
    }

     // pml corrections
    if (dopml) {
      for (int i = 0; i < pmllen; i++) {
	for (int j = 0; j < Ny; j++) {
	  for (int k = 0; k < Nz-1; k++) {
	    Hy[i][j][k] = Hy[i][j][k] + c2h * P_hyx[i][j][k];
	  }
	}
	//Hyinc[i] = Hyinc[i] + c2h * P_hyxI[i];
      }
      
      for (int i = Nx-pmllen-1; i < Nx-1; i++) {
	for (int j = 0; j < Ny; j++) {
	  for (int k = 0; k < Nz-1; k++) {
	    Hy[i][j][k] = Hy[i][j][k] + c2h * P_hyx[i-xshift][j][k];
	  }
	}
        Hyinc[i] = Hyinc[i] + c2h * P_hyxI[i-xshift];
      }

      for (int i = 0; i < Nx-1; i++) {
	for (int j = 0; j < Ny; j++) {
	  for (int k = 0; k < pmllen; k++) {
	    Hy[i][j][k] = Hy[i][j][k] - c2h * P_hyz[i][j][k];
	  }
	  for (int k = Nz-pmllen-1; k < Nz-1; k++) {
	    Hy[i][j][k] = Hy[i][j][k] - c2h * P_hyz[i][j][k-zshift];
	  }
	}
      }

    } // if dopml


    // ----------------------------------------
    // Hz update
#pragma omp parallel for
    for (int i = 0; i < Nx-1; i++) {
      for (int j = 0; j < Ny-1; j++) {
	for (int k = 0; k < Nz; k++) {
	  // TF/SF
	  Hz[i][j][k] = c1h * Hz[i][j][k] - c2h * \
	    ( ( Ey[i+1][j][k] - Ey[i][j][k] ) / dx - ( Ex[i][j+1][k] - Ex[i][j][k] ) / dy );
	}
      }
      // incident field
      Hzinc[i] = c1h * Hzinc[i] - c2h * ( Eyinc[i+1] - Eyinc[i] ) / dx;
    }
    
    if (dotfsf) {
    // TFSF correction: Hz corrected on x, y walls of TFSF box, but there is no Exinc, so no y correction
    for (int j = by1; j < by2; j++) {
      for (int k = bz1; k < bz2; k++) {
	Hz[bx1-1][j][k] = Hz[bx1-1][j][k] + c2h * Eyinc[bx1] / dx;
	Hz[bx2-1][j][k] = Hz[bx2-1][j][k] - c2h * Eyinc[bx2-1] / dx; 
      }
    }
    }

    // pml corrections
    if (dopml) {
      for (int i = 0; i < pmllen; i++) {
	for (int j = 0; j < Ny-1; j++) {
	  for (int k = 0; k < Nz; k++) {
	    Hz[i][j][k] = Hz[i][j][k] - c2h * P_hzx[i][j][k];
	  }
	}
	//Hzinc[i] = Hzinc[i] - c2h * P_hzxI[i];
      }
      for (int i = Nx-pmllen-1; i < Nx-1; i++) {
	for (int j = 0; j < Ny-1; j++) {
	  for (int k = 0; k < Nz; k++) {
	    Hz[i][j][k] = Hz[i][j][k] - c2h * P_hzx[i-xshift][j][k];
	  }
	}
	Hzinc[i] = Hzinc[i] - c2h * P_hzxI[i-xshift];
      }

      for (int i = 0; i < Nx-1; i++) {
	for (int j = 0; j < pmllen; j++) {
	  for (int k = 0; k < Nz; k++) {
	    Hz[i][j][k] = Hz[i][j][k] + c2h * P_hzy[i][j][k];
	  }
	}
	for (int j = Ny-pmllen-1; j < Ny-1; j++) {
	  for (int k = 0; k < Nz; k++) {
	    Hz[i][j][k] = Hz[i][j][k] + c2h * P_hzy[i][j-yshift][k];
	  }
	}
      }

    } // if dopml

    ////////////////////////////////////////////
    // NTF transformation

    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nz; k++) {
	for (int m = 0; m < nntf; m++) {
	  // front
	  Jy_front[j][k][2*m] += 0.5*(Hz[ntfb][j][k] + Hz[ntfb-1][j][k]) * sin(2*PI*fntf[m]*t*dt);
	  Jy_front[j][k][2*m+1] += 0.5*(Hz[ntfb][j][k] + Hz[ntfb-1][j][k]) * cos(2*PI*fntf[m]*t*dt);
	  Jz_front[j][k][2*m] += -0.5*(Hy[ntfb][j][k] + Hy[ntfb-1][j][k]) * sin(2*PI*fntf[m]*t*dt);
	  Jz_front[j][k][2*m+1] += -0.5*(Hy[ntfb][j][k] + Hy[ntfb-1][j][k]) * cos(2*PI*fntf[m]*t*dt);
	  // back
	  Jy_back[j][k][2*m] += -0.5*(Hz[Nx-ntfb][j][k] + Hz[Nx-ntfb-1][j][k]) * sin(2*PI*fntf[m]*t*dt);
	  Jy_back[j][k][2*m+1] += -0.5*(Hz[Nx-ntfb][j][k] + Hz[Nx-ntfb-1][j][k]) * cos(2*PI*fntf[m]*t*dt);
	  Jz_back[j][k][2*m] += 0.5*(Hy[Nx-ntfb][j][k] + Hy[Nx-ntfb-1][j][k]) * sin(2*PI*fntf[m]*t*dt);
	  Jz_back[j][k][2*m+1] += 0.5*(Hy[Nx-ntfb][j][k] + Hy[Nx-ntfb-1][j][k]) * cos(2*PI*fntf[m]*t*dt);
	}
      }
    }    

    for (int i = 0; i < Nx; i++) {
      for (int k = 0; k < Nz; k++) {
	for (int m = 0; m < nntf; m++) {
	  // top 
	  Jx_top[i][k][2*m] += 0.5*(Hz[i][Ny-ntfb][k] + Hz[i][Ny-ntfb-1][k]) * sin(2*PI*fntf[m]*t*dt);
	  Jx_top[i][k][2*m+1] += 0.5*(Hz[i][Ny-ntfb][k] + Hz[i][Ny-ntfb-1][k]) * cos(2*PI*fntf[m]*t*dt);
	  Jz_top[i][k][2*m] += -0.5*(Hx[i][Ny-ntfb][k] + Hx[i][Ny-ntfb-1][k]) * sin(2*PI*fntf[m]*t*dt);
	  Jz_top[i][k][2*m+1] += -0.5*(Hx[i][Ny-ntfb][k] + Hx[i][Ny-ntfb-1][k]) * cos(2*PI*fntf[m]*t*dt);
	  // bottom
	  Jx_bot[i][k][2*m] += -0.5*(Hz[i][ntfb][k] + Hz[i][ntfb-1][k]) * sin(2*PI*fntf[m]*t*dt);
	  Jx_bot[i][k][2*m+1] += -0.5*(Hz[i][ntfb][k] + Hz[i][ntfb-1][k]) * cos(2*PI*fntf[m]*t*dt);
	  Jz_bot[i][k][2*m] += 0.5*(Hx[i][ntfb][k] + Hx[i][ntfb-1][k]) * sin(2*PI*fntf[m]*t*dt);
	  Jz_bot[i][k][2*m+1] += 0.5*(Hx[i][ntfb][k] + Hx[i][ntfb-1][k]) * cos(2*PI*fntf[m]*t*dt);

	}
      }
    }

    for (int i = 0; i < Nx; i++) {
      for (int j = 0; j < Ny; j++) {
	for (int m = 0; m < nntf; m++) {
	  // left
	  Jx_left[i][j][2*m] += 0.5*(Hy[i][j][ntfb] + Hy[i][j][ntfb-1]) * sin(2*PI*fntf[m]*t*dt);
	  Jx_left[i][j][2*m+1] += 0.5*(Hy[i][j][ntfb] + Hy[i][j][ntfb-1]) * cos(2*PI*fntf[m]*t*dt);
	  Jy_left[i][j][2*m] += -0.5*(Hx[i][j][ntfb] + Hx[i][j][ntfb-1]) * sin(2*PI*fntf[m]*t*dt);
	  Jy_left[i][j][2*m+1] += -0.5*(Hx[i][j][ntfb] + Hx[i][j][ntfb-1]) * cos(2*PI*fntf[m]*t*dt);
	  // right
	  Jx_right[i][j][2*m] += -0.5*(Hy[i][j][Nz-ntfb] + Hy[i][j][Nz-ntfb-1]) * sin(2*PI*fntf[m]*t*dt);
	  Jx_right[i][j][2*m+1] += -0.5*(Hy[i][j][Nz-ntfb] + Hy[i][j][Nz-ntfb-1]) * cos(2*PI*fntf[m]*t*dt);
	  Jy_right[i][j][2*m] += 0.5*(Hx[i][j][Nz-ntfb] + Hx[i][j][Nz-ntfb-1]) * sin(2*PI*fntf[m]*t*dt);
	  Jy_right[i][j][2*m+1] += 0.5*(Hx[i][j][Nz-ntfb] + Hx[i][j][Nz-ntfb-1]) * cos(2*PI*fntf[m]*t*dt);
	}
      }
    }

    // ---------------------------------------

    // spit it out to file
    for (int i = 0; i < Nx; i++) {
      for (int j = 0; j < Ny-1; j++) {
	Hslicez[0][i][j] = Hy[i][j][ssz];
	Hslicez[1][i][j] = Hz[i][j][ssz];
      }
    }
    for (int i = 0; i < Nx; i++) {
      for (int k = 0; k < Nz-1; k++) {
	Hslicey[0][i][k] = Hy[i][ssy][k];
	Hslicey[1][i][k] = Hz[i][ssy][k];
      }
    }
    for (int j = 0; j < Ny-1; j++) {
      for (int k = 0; k < Nz-1; k++) {
	Hslicex[0][j][k] = Hy[ssx][j][k];
	Hslicex[1][j][k] = Hz[ssx][j][k];
      }
    }

    ////////////////////////////////////////////////////////////////////////////////
    // ---------------------------------------
    // Psi updates for E field

    if (dopml) {

      // Psi's applied to Ex
      for (int i = 0; i < Nx-1; i++) {
	for (int j = 1; j < pmllen+1; j++) {
	  for (int k = 1; k < Nz-1; k++) {
	    P_exy[i][j-1][k] = B_exy[j-1] * P_exy[i][j-1][k] \
	      + A_exy[j-1] * ( Hz[i][j][k] - Hz[i][j-1][k] ) / dy;
	  }
	}
	for (int j = Ny-pmllen-1; j < Ny-1; j++) {
	  for (int k = 1; k < Nz-1; k++) {
	    P_exy[i][j-yshift][k] = B_exy[j-yshift] * P_exy[i][j-yshift][k] \
	      + A_exy[j-yshift] * ( Hz[i][j][k] - Hz[i][j-1][k] ) / dy;
	  }
	}
      }

      for (int i = 0; i < Nx; i++) {
	for (int j = 1; j < Ny-1; j++) {
	  for (int k = 1; k < pmllen+1; k++) {
	    P_exz[i][j][k-1] = B_exz[k-1] * P_exz[i][j][k-1] \
	      + A_exz[k-1] * ( Hy[i][j][k] - Hy[i][j][k-1] ) / dz;
	  }
	  for (int k = Nz-pmllen-1; k < Nz-1; k++) {
	    P_exz[i][j][k-zshift] = B_exz[k-zshift] * P_exz[i][j][k-zshift] \
	      + A_exz[k-zshift] * ( Hy[i][j][k] - Hy[i][j][k-1] ) / dz;
	  }
	}
      }

      // Psi's applied to Ey
      for (int i = 1; i < pmllen+1; i++) {
	for (int j = 0; j < Ny-1; j++) {
	  for (int k = 1; k < Nz-1; k++) {
	    P_eyx[i-1][j][k] = B_eyx[i-1] * P_eyx[i-1][j][k] \
	      + A_eyx[i-1] * ( Hz[i][j][k] - Hz[i-1][j][k] ) / dx; 
	  }
	}
	/*P_eyxI[i-1] = B_eyx[i-1] * P_eyxI[i-1]	\
	  + A_eyx[i-1] * ( Hzinc[i] - Hzinc[i-1] ) / dx; */
      }      
      
      for (int i = Nx-pmllen-1; i < Nx-1; i++) {
	for (int j = 0; j < Ny-1; j++) {
	  for (int k = 1; k < Nz-1; k++) {
	    P_eyx[i-xshift][j][k] = B_eyx[i-xshift] * P_eyx[i-xshift][j][k] \
	      + A_eyx[i-xshift] * ( Hz[i][j][k] - Hz[i-1][j][k] ) / dx; 
	  }
	}
	P_eyxI[i-xshift] = B_eyx[i-xshift] * P_eyxI[i-xshift] \
	  + A_eyx[i-xshift] * ( Hzinc[i] - Hzinc[i-1] ) / dx; 
      }

      for (int i = 1; i < Nx-1; i++) {
	for (int j = 0; j < Ny-1; j++) {
	  for (int k = 1; k < pmllen+1; k++) {
	    P_eyz[i][j][k-1] = B_eyz[k-1] * P_eyz[i][j][k-1] \
	      + A_eyz[k-1] * ( Hx[i][j][k] - Hx[i][j][k-1] ) / dz;
	  }
	  for (int k = Nz-pmllen-1; k < Nz-1; k++) {
	    P_eyz[i][j][k-zshift] = B_eyz[k-zshift] * P_eyz[i][j][k-zshift] \
	      + A_eyz[k-zshift] * ( Hx[i][j][k] - Hx[i][j][k-1] ) / dz;
	  }
	}
      }

      // Psi's applied to Ez
      for (int i = 1; i < pmllen+1; i++) {
	for (int j = 1; j < Ny-1; j++) {
	  for (int k = 0; k < Nz-1; k++) {
	    P_ezx[i-1][j][k] = B_ezx[i-1] * P_ezx[i-1][j][k] \
	      + A_ezx[i-1] * ( Hy[i][j][k] - Hy[i-1][j][k] ) / dx;
	  }
	}
	/*P_ezxI[i-1] = B_ezx[i-1] * P_ezxI[i-1]		\
	  + A_ezx[i-1] * ( Hyinc[i] - Hyinc[i-1] ) / dx; */
      }
      
      for (int i = Nx-pmllen-1; i < Nx-1; i++) {
	for (int j = 1; j < Ny-1; j++) {
	  for (int k = 0; k < Nz-1; k++) {
	    P_ezx[i-xshift][j][k] = B_ezx[i-xshift] * P_ezx[i-xshift][j][k] \
	      + A_ezx[i-xshift] * ( Hy[i][j][k] - Hy[i-1][j][k] ) / dx;
	  }
	}
	P_ezxI[i-xshift] = B_ezx[i-xshift] * P_ezxI[i-xshift] \
	  + A_ezx[i-xshift] * ( Hyinc[i] - Hyinc[i-1] ) / dx;
      }
 
      for (int i = 1; i < Nx-1; i++) {
	for (int j = 1; j < pmllen+1; j++) {
	  for (int k = 0; k < Nz-1; k++) {
	    P_ezy[i][j-1][k] = B_ezy[j-1] * P_ezy[i][j-1][k] \
	      + A_ezy[j-1] * ( Hx[i][j][k] - Hx[i][j-1][k] ) / dy;
	  }
	}
	for (int j = Ny-pmllen-1; j < Ny-1; j++) {
	  for (int k = 0; k < Nz-1; k++) {
	    P_ezy[i][j-yshift][k] = B_ezy[j-yshift] * P_ezy[i][j-yshift][k] \
	      + A_ezy[j-yshift] * (Hx[i][j][k] - Hx[i][j-1][k] ) / dy;
	  }
	}
      }

    } // if dopml


    // ----------------------------------------
    // Ex update
#pragma omp parallel for
    for (int i = 0; i < Nx-1; i++) {
      for (int j = 1; j < Ny-1; j++) {
	for (int k = 1; k < Nz-1; k++) {
	  // TF/SF
	  Ex[i][j][k] = c1e[i][j][k] * Ex[i][j][k] + c2e[i][j][k] * \
	    ( ( Hz[i][j][k] - Hz[i][j-1][k] ) / dy - ( Hy[i][j][k] - Hy[i][j][k-1] ) / dz ) \
	    - c2e[i][j][k] * Jex[i][j][k];
	}
      }
      // no incident Ex in 1D
    }

    if (dotfsf) {
    // TF/SF correction: x fields on y, z walls
    for (int i = bx1; i < bx2-1; i++) {
      for (int j = by1; j < by2; j++) {
	Ex[i][j][bz1] = Ex[i][j][bz1] + c2e0 * Hyinc[i] / dz;
	Ex[i][j][bz2] = Ex[i][j][bz2] - c2e0 * Hyinc[i] / dz; 
      }
    }
    for (int i = bx1; i < bx2-1; i++) {
      for (int k = bz1; k < bz2; k++) {
	Ex[i][by1][k] = Ex[i][by1][k] - c2e0 * Hzinc[i] / dy;
	Ex[i][by2][k] = Ex[i][by2][k] + c2e0 * Hzinc[i] / dy;
      }
    }
    }

    // pml corrections for Ex
    if (dopml) {
      for (int i = 0; i < Nx-1; i++) {
	for (int j = 1; j < pmllen+1; j++) {
	  for (int k = 1; k < Nz-1; k++) {
	    Ex[i][j][k] = Ex[i][j][k] + c2e0 * P_exy[i][j-1][k];
	  }
	}
	for (int j = Ny-pmllen-1; j < Ny-1; j++) {
	  for (int k = 1; k < Nz-1; k++) {
	    Ex[i][j][k] = Ex[i][j][k] + c2e0 * P_exy[i][j-yshift][k];
	  }
	}
      }

      for (int i = 0; i < Nx-1; i++) {
	for (int j = 1; j < Ny-1; j++) {
	  for (int k = 1; k < pmllen+1; k++) {
	    Ex[i][j][k] = Ex[i][j][k] - c2e0 * P_exz[i][j][k-1];
	  }
	  for (int k = Nz-pmllen-1; k < Nz-1; k++) {
	    Ex[i][j][k] = Ex[i][j][k] - c2e0 * P_exz[i][j][k-zshift];
	  }
	}
      }

    } // if dopml


    // ----------------------------------------
    // Ey update
#pragma omp parallel for
    for (int i = 1; i < Nx-1; i++) {
      for (int j = 0; j < Ny; j++) {
	for (int k = 1; k < Nz; k++) {
	  // TF/SF
	  Ey[i][j][k] = c1e[i][j][k] * Ey[i][j][k] + c2e[i][j][k] *	\
	    ( ( Hx[i][j][k] - Hx[i][j][k-1] ) / dz - ( Hz[i][j][k] - Hz[i-1][j][k] ) / dx ) \
	    - c2e[i][j][k] * Jey[i][j][k];
	}
      }
      // incident field
      Eyinc[i] = c1e0 * Eyinc[i] - c2e0 * ( Hzinc[i] - Hzinc[i-1] ) / dx;
    }

    if (dotfsf) {
      // TF/SF correction: y fields on x, z walls: where are z walls??
      for (int j = by1; j < by2; j++) {
	for (int k = bz1; k < bz2; k++) {
	  Ey[bx1][j][k] = Ey[bx1][j][k] + c2e0 * Hzinc[bx1-1] / dx;
	  Ey[bx2-1][j][k] = Ey[bx2-1][j][k] - c2e0 * Hzinc[bx2-1] / dx; 
	}
      }
    }

    // pml corrections
    if (dopml) {
      for (int i = 1; i < pmllen+1; i++) {
	for (int j = 0; j < Ny; j++) {
	  for (int k = 1; k < Nz; k++) {
	    Ey[i][j][k] = Ey[i][j][k] - c2e0 * P_eyx[i-1][j][k];
	  }
	}
	//Eyinc[i] = Eyinc[i] - c2e * P_eyxI[i-1];
      }
      for (int i = Nx-pmllen-1; i < Nx-1; i++) {
	for (int j = 0; j < Ny; j++) {
	  for (int k = 1; k < Nz; k++) {
	    Ey[i][j][k] = Ey[i][j][k] - c2e0 * P_eyx[i-xshift][j][k];
	  }
	}
	Eyinc[i] = Eyinc[i] - c2e0 * P_eyxI[i-xshift];
      }

      for (int i = 1; i < Nx-1; i++) {
	for (int j = 0; j < Ny; j++) {
	  for (int k = 1; k < pmllen+1; k++) {
	    Ey[i][j][k] = Ey[i][j][k] + c2e0 * P_eyz[i][j][k-1];
	  }
	  for (int k = Nz-pmllen-1; k < Nz-1; k++) {
	    Ey[i][j][k] = Ey[i][j][k] + c2e0 * P_eyz[i][j][k-zshift];
	  }
	}
      }

    } // if dopml

    // ----------------------------------------
    // Ez update
#pragma omp parallel for
    for (int i = 1; i < Nx-1; i++) {
      for (int j = 1; j < Ny; j++) {
	for (int k = 0; k < Nz; k++) {
	  // TF/SF
	  Ez[i][j][k] = c1e[i][j][k] * Ez[i][j][k] + c2e[i][j][k] * \
	    ( ( Hy[i][j][k] - Hy[i-1][j][k] ) / dx - \
	      ( Hx[i][j][k] - Hx[i][j-1][k] ) / dy ) - c2e[i][j][k] * Jez[i][j][k];
	}
      }
      // incident field
      Ezinc[i] = c1e0 * Ezinc[i] + c2e0 * ( Hyinc[i] - Hyinc[i-1] ) / dx;
    }


    if (dotfsf) {
      // TF/SF correction: z fields on x, y walls?
      for (int j = by1; j < by2; j++) {
	for (int k = bz1; k < bz2; k++) {
	  Ez[bx1][j][k] = Ez[bx1][j][k] - c2e0 * Hyinc[bx1-1] / dx;
	  Ez[bx2-1][j][k] = Ez[bx2-1][j][k] + c2e0 * Hyinc[bx2-1] / dx; 
	}
      }
    }

    // Ez source?
    //Ez[ssx][ssy][ssz] = cos(2*PI*f0*(t-3*fsig)*dt) * exp(-(t-3*fsig)*(t-3*fsig)/(2*fsig*fsig));
    //Ey[ssx][ssy][ssz] = cos(2*PI*f0*(t-3*fsig)*dt) * exp(-(t-3*fsig)*(t-3*fsig)/(2*fsig*fsig));

    // choose incident field source based on polarization
    if (1) {
	  if (sourcepol == 1) {
	    // linear, Ez polarization
	    Ezinc[0] = EsourceA[t];
	    Eyinc[0] = 0;
	  } else if (sourcepol == 2) {
	    // linear, Ey polarization
	    Ezinc[0] = 0;
	    Eyinc[0] = EsourceA[t];
	  } else if (sourcepol == 3) { 
	    // 45 degrees linear
	    Ezinc[0] = EsourceA[t];
	    Eyinc[0] = EsourceA[t];
	  } else if (sourcepol == 4) {
	    // RHCP
	    Ezinc[0] = EsourceA[t];
	    Eyinc[0] = EsourceB[t];
	  } else if (sourcepol == 5) {
	    // LHCP
	    Ezinc[0] = EsourceB[t];
	    Eyinc[0] = EsourceA[t];
	  } else {
	    // default to linear Ez
	    Ezinc[0] = EsourceA[t];
	    Eyinc[0] = 0;
	  }
    }

    // plate scatterer
    if (doplate == 1) {
      int lowerlim = (int)(ssy-platesize/2/dy);
      int upperlim = (int)(ssy+platesize/2/dy);
      for (int j = lowerlim; j < upperlim-1; j++) {
	for (int k = lowerlim; k < upperlim; k++) {
	  Ey[ssx][j][k] = 0;
	}
      }
      for (int j = lowerlim; j < upperlim; j++) {
	for (int k = lowerlim; k < upperlim-1; k++) {
	  Ez[ssx][j][k] = 0;
	}
      }
    }
    
        // spit it out to file
    for (int i = 0; i < Nx; i++) {
      for (int j = 0; j < Ny; j++) {
	Eslicez[0][i][j] = Ey[i][j][ssz];
	Eslicez[1][i][j] = Ez[i][j][ssz];
      }
    }
    for (int i = 0; i < Nx; i++) {
      for (int k = 0; k < Nz; k++) {
	Eslicey[0][i][k] = Ey[i][ssy][k];
	Eslicey[1][i][k] = Ez[i][ssy][k];
      }
    }
    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nz; k++) {
	Eslicex[0][j][k] = Ey[ssx][j][k];
	Eslicex[1][j][k] = Ez[ssx][j][k];
      }
    }

    // pml corrections
    if (dopml) {
      for (int i = 1; i < pmllen+1; i++) {
	for (int j = 1; j < Ny; j++) {
	  for (int k = 0; k < Nz; k++) {
	    Ez[i][j][k] = Ez[i][j][k] + c2e0 * P_ezx[i-1][j][k];
	  }
	}
	//Ezinc[i] = Ezinc[i] + c2e * P_ezxI[i-1];
      }
      for (int i = Nx-pmllen-1; i < Nx-1; i++) {
	for (int j = 1; j < Ny; j++) {
	  for (int k = 0; k < Nz; k++) {
	    Ez[i][j][k] = Ez[i][j][k] + c2e0 * P_ezx[i-xshift][j][k];
	  }
	}
	Ezinc[i] = Ezinc[i] + c2e0 * P_ezxI[i-xshift];
      }

      for (int i = 1; i < Nx-1; i++) {
	for (int j = 1; j < pmllen+1; j++) {
	  for (int k = 0; k < Nz-1; k++) {
	    Ez[i][j][k] = Ez[i][j][k] - c2e0 * P_ezy[i][j-1][k];
	  }
	}
	for (int j = Ny-pmllen-1; j < Ny-1; j++) {
	  for (int k = 0; k < Nz-1; k++) {
	    Ez[i][j][k] = Ez[i][j][k] - c2e0 * P_ezy[i][j-yshift][k];
	  }
	}
      }

    } // if dopml

    // need to calculate and save DFT of source at NTF frequencies, to be used in NTF calculation.

    for (int m = 0; m < nntf; m++) {
      Ezdft[2*m] += Ezinc[0] * sin(2*PI*fntf[m]*(t+1/2)*dt);
      Ezdft[2*m+1] += Ezinc[0] * cos(2*PI*fntf[m]*(t+1/2)*dt);
      Eydft[2*m] += Eyinc[0] * sin(2*PI*fntf[m]*(t+1/2)*dt);
      Eydft[2*m+1] += Eyinc[0] * cos(2*PI*fntf[m]*(t+1/2)*dt);
    }

    ////////////////////////////////////////////                                                     
    // NTF transformation                                                                            

    for (int j = 0; j < Ny; j++) {
      for (int k = 0; k < Nz; k++) {
        for (int m = 0; m < nntf; m++) {
          // front                                                                                   
          My_front[j][k][2*m] += -Ez[ntfb][j][k] * sin(2*PI*fntf[m]*(t+1/2)*dt);
          My_front[j][k][2*m+1] += -Ez[ntfb][j][k] * cos(2*PI*fntf[m]*(t+1/2)*dt);
          Mz_front[j][k][2*m] += Ey[ntfb][j][k] * sin(2*PI*fntf[m]*(t+1/2)*dt);
          Mz_front[j][k][2*m+1] += Ey[ntfb][j][k] * cos(2*PI*fntf[m]*(t+1/2)*dt);
          // back                                                                                    
          My_back[j][k][2*m] += Ez[Nx-ntfb][j][k] * sin(2*PI*fntf[m]*(t+1/2)*dt);
          My_back[j][k][2*m+1] += Ez[Nx-ntfb][j][k] * cos(2*PI*fntf[m]*(t+1/2)*dt);
          Mz_back[j][k][2*m] += -Ey[Nx-ntfb][j][k] * sin(2*PI*fntf[m]*(t+1/2)*dt);
          Mz_back[j][k][2*m+1] += -Ey[Nx-ntfb][j][k] * cos(2*PI*fntf[m]*(t+1/2)*dt);
        }
      }
    }

    for (int i = 0; i < Nx; i++) {
      for (int k = 0; k < Nz; k++) {
        for (int m = 0; m < nntf; m++) {
          // top                                                                                     
          Mx_top[i][k][2*m] += -Ez[i][Ny-ntfb][k] * sin(2*PI*fntf[m]*(t+1/2)*dt);
          Mx_top[i][k][2*m+1] += -Ez[i][Ny-ntfb][k] * cos(2*PI*fntf[m]*(t+1/2)*dt);
          Mz_top[i][k][2*m] += Ex[i][Ny-ntfb][k] * sin(2*PI*fntf[m]*(t+1/2)*dt);
          Mz_top[i][k][2*m+1] += Ex[i][Ny-ntfb][k] * cos(2*PI*fntf[m]*(t+1/2)*dt);
          // bottom                                                                                  
          Mx_bot[i][k][2*m] += Ez[i][ntfb][k] * sin(2*PI*fntf[m]*(t+1/2)*dt);
          Mx_bot[i][k][2*m+1] += Ez[i][ntfb][k] * cos(2*PI*fntf[m]*(t+1/2)*dt);
          Mz_bot[i][k][2*m] += -Ex[i][ntfb][k] * sin(2*PI*fntf[m]*(t+1/2)*dt);
          Mz_bot[i][k][2*m+1] += -Ex[i][ntfb][k] * cos(2*PI*fntf[m]*(t+1/2)*dt);

        }
      }
    }

    for (int i = 0; i < Nx; i++) {
      for (int j = 0; j < Ny; j++) {
        for (int m = 0; m < nntf; m++) {
          // left                                                                                    
          Mx_left[i][j][2*m] += -Ey[i][j][ntfb] * sin(2*PI*fntf[m]*(t+1/2)*dt);
          Mx_left[i][j][2*m+1] += -Ey[i][j][ntfb] * cos(2*PI*fntf[m]*(t+1/2)*dt);
          My_left[i][j][2*m] += Ex[i][j][ntfb] * sin(2*PI*fntf[m]*(t+1/2)*dt);
          My_left[i][j][2*m+1] += Ex[i][j][ntfb] * cos(2*PI*fntf[m]*(t+1/2)*dt);
          // right                                                                                   
          Mx_right[i][j][2*m] += Ey[i][j][Nz-ntfb] * sin(2*PI*fntf[m]*(t+1/2)*dt);
	  Mx_right[i][j][2*m+1] += Ey[i][j][Nz-ntfb] * cos(2*PI*fntf[m]*(t+1/2)*dt);
          My_right[i][j][2*m] += -Ex[i][j][Nz-ntfb] * sin(2*PI*fntf[m]*(t+1/2)*dt);
          My_right[i][j][2*m+1] += -Ex[i][j][Nz-ntfb] * cos(2*PI*fntf[m]*(t+1/2)*dt);
        }
      }
    }

    // ----------------------------------------------------
    // output E and H probe points at desired points

    for (int i = 0; i < nprobes; i++) {
      ix = probex[i];
      iy = probey[i];
      iz = probez[i];
      Exprobe[t][i] = 0.5*(Ex[ix][iy][iz]+Ex[ix-1][iy][iz]);
      Eyprobe[t][i] = 0.5*(Ey[ix][iy][iz]+Ey[ix][iy-1][iz]);
      Ezprobe[t][i] = 0.5*(Ez[ix][iy][iz]+Ez[ix][iy][iz-1]);

      Hxprobe[t][i] = 0.25*(Hx[ix][iy][iz]+Hx[ix][iy-1][iz]+Hx[ix][iy][iz-1]+Hx[ix][iy-1][iz-1]);
      Hyprobe[t][i] = 0.25*(Hy[ix][iy][iz]+Hy[ix-1][iy][iz]+Hy[ix][iy][iz-1]+Hy[ix-1][iy][iz-1]);
      Hzprobe[t][i] = 0.25*(Hz[ix][iy][iz]+Hz[ix][iy-1][iz]+Hz[ix-1][iy][iz]+Hz[ix-1][iy-1][iz]);

    }

    // ---------------------------------------
    // update E_eff, for use in ionization updates

#pragma omp parallel for private(Sx,Sy,Sz,Exm,Eym,Ezm,Hxm,Hym,Hzm)
    for (int i = 1; i < Nx-1; i++) {
      for (int j = 1; j < Ny-1; j++) {
	for (int k = 1; k < Nz-1; k++) {

	  Exm = (Ex[i][j][k] + Ex[i-1][j][k]) / 2;
	  Eym = (Ey[i][j][k] + Ey[i][j-1][k]) / 2;
	  Ezm = (Ez[i][j][k] + Ez[i][j][k-1]) / 2;
	  Hxm = (Hx[i][j][k] + Hx[i][j-1][k] + Hx[i][j][k-1] + Hx[i][j-1][k-1]) / 4;
	  Hym = (Hy[i][j][k] + Hy[i-1][j][k] + Hy[i][j][k-1] + Hy[i-1][j][k-1]) / 4;
	  Hzm = (Hz[i][j][k] + Hz[i-1][j][k] + Hz[i][j-1][k] + Hz[i-1][j-1][k]) / 4;
	  Emag[i][j][k] = sqrt( Exm*Exm + Eym*Eym + Ezm*Ezm );

	  Sx = Eym*Hzm - Ezm*Hym;
	  Sy = Ezm*Hxm - Exm*Hzm;
	  Sz = Exm*Hym - Eym*Hxm;
	  S[i][j][k] = sqrt( Sx*Sx + Sy*Sy + Sz*Sz );

	}
      }
    }


    for (int i = 0; i < Nx; i++) {
      for (int j = 0; j < Ny; j++) {
	Emagslice[i][j] = Emag[i][j][ssz];
      }	
    }

    // -----------------------------------------
    // J updates -

    if (!doB0) { // Young's DI method

#pragma omp parallel for
    for (int i = 1; i < Nx-1; i++) {
      for (int j = 1; j < Ny-1; j++) {
	for (int k = 1; k < Nz-1; k++) {
	  Jex[i][j][k] = ((2-nue[i][j][k]*dt)/(2+nue[i][j][k]*dt)) * Jex[i][j][k] + \
	    2 * ((E0*wpe[i][j][k]*wpe[i][j][k]*dt)/(2+nue[i][j][k]*dt)) * Ex[i][j][k];
	  Jey[i][j][k] = ((2-nue[i][j][k]*dt)/(2+nue[i][j][k]*dt)) * Jey[i][j][k] + \
	    2 * ((E0*wpe[i][j][k]*wpe[i][j][k]*dt)/(2+nue[i][j][k]*dt)) * Ey[i][j][k];
	  Jez[i][j][k] = ((2-nue[i][j][k]*dt)/(2+nue[i][j][k]*dt)) * Jez[i][j][k] + \
	    2 * ((E0*wpe[i][j][k]*wpe[i][j][k]*dt)/(2+nue[i][j][k]*dt)) * Ez[i][j][k];
	}
      }
    }

    } else { // Lee and Kalluri method

#pragma omp parallel for private(Ee1,Ee2,Se1,Ce1,Ce2,Ce3,Ce4,Ae,Ke,Emid,Jex0,Jey0)
      for (int i = 1; i < Nx-1; i++) {
	for (int j = 1; j < Ny-1; j++) {
	  for (int k = 1; k < Nz-1; k++) {
	    
	    Ee1 = exp(-nue[i][j][k]*dt);
	    Ee2 = 1/(wce0*wce0 + nue[i][j][k]*nue[i][j][k]);
	    Se1 = sin(wce0*dt)/wce0;
	    Ce1 = (1 - cos(wce0*dt))/(wce0*wce0);
	    Ce2 = (1 - Ee1)/nue[i][j][k] - Ee1*nue[i][j][k]*Ce1 - Ee1*Se1;
	    Ce3 = nue[i][j][k] * (1 - Ee1*cos(wce0*dt)) + Ee1*wce0*sin(wce0*dt);
	    Ce4 = 1 - Ee1*cos(wce0*dt) - Ee1*nue[i][j][k]*Se1;
	    
	    // A and K matrices
	    Ae[0][0] = Ee1 * ( Ce1*wcex*wcex + cos(wce0*dt) );
	    Ae[0][1] = Ee1 * ( Ce1*wcex*wcey - Se1*wcez );
	    Ae[0][2] = Ee1 * ( Ce1*wcex*wcez + Se1*wcey );
	    Ae[1][0] = Ee1 * ( Ce1*wcey*wcex + Se1*wcez );     
	    Ae[1][1] = Ee1 * ( Ce1*wcey*wcey + cos(wce0*dt) );
	    Ae[1][2] = Ee1 * ( Ce1*wcey*wcez - Se1*wcex );
	    Ae[2][0] = Ee1 * ( Ce1*wcez*wcex - Se1*wcey );
	    Ae[2][1] = Ee1 * ( Ce1*wcez*wcey + Se1*wcex );
	    Ae[2][2] = Ee1 * ( Ce1*wcez*wcez + cos(wce0*dt) );
            
	    Ke[0][0] = Ee2 * ( Ce2*wcex*wcex + Ce3 );
	    Ke[0][1] = Ee2 * ( Ce2*wcex*wcey - Ce4*wcez );
	    Ke[0][2] = Ee2 * ( Ce2*wcex*wcez + Ce4*wcey );
	    Ke[1][0] = Ee2 * ( Ce2*wcey*wcex + Ce4*wcez );
	    Ke[1][1] = Ee2 * ( Ce2*wcey*wcey + Ce3 );
	    Ke[1][2] = Ee2 * ( Ce2*wcey*wcez - Ce4*wcex );
	    Ke[2][0] = Ee2 * ( Ce2*wcez*wcex - Ce4*wcey );
	    Ke[2][1] = Ee2 * ( Ce2*wcez*wcey + Ce4*wcex );
	    Ke[2][2] = Ee2 * ( Ce2*wcez*wcez + Ce3 );

	    // okay, updates for J finally.

	    Emid[0] = (Ex[i][j][k] + Ex[i-1][j][k])/2;
	    Emid[1] = (Ey[i][j][k] + Ey[i][j-1][k])/2;
	    Emid[2] = (Ez[i][j][k] + Ez[i][j][k-1])/2;
	    
	    Jex0 = Ae[0][0] * Jex[i][j][k] + Ae[0][1] * Jey[i][j][k] + Ae[0][2] * Jez[i][j][k] \
	      + E0 * wpe[i][j][k]*wpe[i][j][k] * ( Ke[0][0] * Emid[0] +  Ke[0][1] * Emid[1] + Ke[0][2] * Emid[2] );
	    
	    Jey0 = Ae[1][0] * Jex[i][j][k] + Ae[1][1] * Jey[i][j][k] + Ae[1][2] * Jez[i][j][k] \
	      + E0 * wpe[i][j][k]*wpe[i][j][k] * ( Ke[1][0] * Emid[0] +  Ke[1][1] * Emid[1] + Ke[1][2] * Emid[2] );
	    
	    Jez[i][j][k] = Ae[2][0] * Jex[i][j][k] + Ae[2][1] * Jey[i][j][k] + Ae[2][2] * Jez[i][j][k] \
	      + E0 * wpe[i][j][k]*wpe[i][j][k] * ( Ke[2][0] * Emid[0] +  Ke[2][1] * Emid[1] + Ke[2][2] * Emid[2] );
	    
	    Jex[i][j][k] = Jex0;
	    Jey[i][j][k] = Jey0;
	    
	  }
	}
      }  

    }


    ////////////////////////////////////////////////////////////////////////
    // ----------------------------------------

    // Figure out run time

    if (t == 20) {
      double runtimemin = (omp_get_wtime() - tloopStart) / 60.0;
      double totaltime = runtimemin * tsteps / 20.0;

      logfile.open("log.txt",std::fstream::app);
      logfile << "t = 20, time taken = " << runtimemin << " minutes.\n";
      logfile << "You can expect the total simulation to take " << totaltime << " minutes.\n";
      logfile.close();

    }

    if (t % (int)(tsteps/numfiles) == 0) {
      partialtime = 100 * t/tsteps;    
      logfile.open("log.txt",std::fstream::app);
      logfile << "t = " << t << ": " << partialtime << " % Done...  ";

      double Emax = 0;
      for (int i = 0; i < Nx; i++) {
	for (int j = 0; j < Ny; j++) {
	  for(int k = 0; k < Nz; k++) {
	    if (fabs(Ex[i][j][k]) > Emax) {
	      Emax = fabs(Ex[i][j][k]);
	    }
	  }
	}
      }
      logfile << "Maximum abs(Ex) is " << Emax << "\n";
      logfile.close();

      // write to files

      if (savefields) {
	FILE * filePtr;

	// E, J, and H are sliced in three directions, using the component selected.
	filePtr = fopen("output_E.dat","ab");
	fwrite(Eslicez,sizeof(double),Nx*Ny*2,filePtr);
	fwrite(Eslicey,sizeof(double),Nx*Nz*2,filePtr);
	fwrite(Eslicex,sizeof(double),Ny*Nz*2,filePtr);
	fwrite(Eyinc,sizeof(double),Nx,filePtr);
	fwrite(Ezinc,sizeof(double),Nx,filePtr);
	fclose(filePtr);

	filePtr = fopen("output_J.dat","ab");
	fwrite(Jslicez,sizeof(double),Nx*Ny*2,filePtr);
	fwrite(Jslicey,sizeof(double),Nx*Nz*2,filePtr);
	fwrite(Jslicex,sizeof(double),Ny*Nz*2,filePtr);
	fclose(filePtr);

	filePtr = fopen("output_H.dat","ab");
	fwrite(Hslicez,sizeof(double),Nx*Ny*2,filePtr);
	fwrite(Hslicey,sizeof(double),Nx*Nz*2,filePtr);
	fwrite(Hslicex,sizeof(double),Ny*Nz*2,filePtr);
	fwrite(Hyinc,sizeof(double),Nx,filePtr);
	fwrite(Hzinc,sizeof(double),Nx,filePtr);
	fclose(filePtr);

      }

    }

    // end of big time loop
  }

  /////////////////////////////////////////////////////////////////////////////

  // write out NTF fields
  int fbsize = Ny*Nz*2*nntf;
  int tbsize = Nx*Nz*2*nntf;
  int lrsize = Nx*Ny*2*nntf;

  FILE * ntfFile;
  ntfFile = fopen("output_NTF.dat","wb");
  fwrite(&nntf,sizeof(int),1,ntfFile);
  fwrite(&ntfb,sizeof(int),1,ntfFile);
  fwrite(&fntf,sizeof(double),nntf,ntfFile);
  fwrite(&Eydft,sizeof(double),2*nntf,ntfFile);
  fwrite(&Ezdft,sizeof(double),2*nntf,ntfFile);

  fwrite(&My_front,sizeof(double),fbsize,ntfFile);
  fwrite(&Mz_front,sizeof(double),fbsize,ntfFile);
  fwrite(&Jy_front,sizeof(double),fbsize,ntfFile);
  fwrite(&Jz_front,sizeof(double),fbsize,ntfFile);
  fwrite(&My_back,sizeof(double),fbsize,ntfFile);
  fwrite(&Mz_back,sizeof(double),fbsize,ntfFile);
  fwrite(&Jy_back,sizeof(double),fbsize,ntfFile);
  fwrite(&Jz_back,sizeof(double),fbsize,ntfFile);

  // test to verify dimensions
  cout << "My-front at y = 40, z = 80 is " << My_front[39][79][0] << "\n";
  cout << "Mz-front at y = 40, z = 80 is " << Mz_front[39][79][0] << "\n";

  fwrite(&Mx_top,sizeof(double),tbsize,ntfFile);
  fwrite(&Mz_top,sizeof(double),tbsize,ntfFile);
  fwrite(&Jx_top,sizeof(double),tbsize,ntfFile);
  fwrite(&Jz_top,sizeof(double),tbsize,ntfFile);
  fwrite(&Mx_bot,sizeof(double),tbsize,ntfFile);
  fwrite(&Mz_bot,sizeof(double),tbsize,ntfFile);
  fwrite(&Jx_bot,sizeof(double),tbsize,ntfFile);
  fwrite(&Jz_bot,sizeof(double),tbsize,ntfFile);

  fwrite(&Mx_left,sizeof(double),lrsize,ntfFile);
  fwrite(&My_left,sizeof(double),lrsize,ntfFile);
  fwrite(&Jx_left,sizeof(double),lrsize,ntfFile);
  fwrite(&Jy_left,sizeof(double),lrsize,ntfFile);
  fwrite(&Mx_right,sizeof(double),lrsize,ntfFile);
  fwrite(&My_right,sizeof(double),lrsize,ntfFile);
  fwrite(&Jx_right,sizeof(double),lrsize,ntfFile);
  fwrite(&Jy_right,sizeof(double),lrsize,ntfFile);
  fclose(ntfFile);

  // write to probe file

  FILE * ProbeFile;
  ProbeFile = fopen("Probe.dat","wb");
  fwrite(&nprobes,sizeof(int),1,ProbeFile);
  fwrite(&probex,sizeof(int),nprobes,ProbeFile);
  fwrite(&probey,sizeof(int),nprobes,ProbeFile);
  fwrite(&probez,sizeof(int),nprobes,ProbeFile);
  fwrite(&Exprobe,sizeof(double),tsteps*nprobes,ProbeFile);
  fwrite(&Eyprobe,sizeof(double),tsteps*nprobes,ProbeFile);
  fwrite(&Ezprobe,sizeof(double),tsteps*nprobes,ProbeFile);
  fwrite(&Hxprobe,sizeof(double),tsteps*nprobes,ProbeFile);
  fwrite(&Hyprobe,sizeof(double),tsteps*nprobes,ProbeFile);
  fwrite(&Hzprobe,sizeof(double),tsteps*nprobes,ProbeFile);
  fclose(ProbeFile);

  logfile.open("log.txt",std::fstream::app);
  logfile << "All done!\n";

  double totalruntime = (omp_get_wtime() - tStart) / 60.0;

  logfile << "Total run time = " << totalruntime << "minutes.\n";
  logfile.close();

  // end program - do not type below this!

}

//////////////////////////////////////////////////////////
