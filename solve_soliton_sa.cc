//Simulated Annealing kink 
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include "mersenne.h"

#define PI 3.1443926584 //  
#define MESH 70 // number of division of mesh
#define MIN -10.0 // minimum of x, y
#define MAX 10.0 // maximum of x, y  
#define START 200001 // total number of the calculation
#define TSTART 5 // Permit the transition to the direction of raising cost in every TSTART
#define KAKUNIN 1000 // Output Energy and Charge in every kakunin

using namespace std;

int M1_2=(MESH+1)*(MESH+1);
double h; // 
double rh; // inverse of h 
double tF[MESH],dx_F[MESH]; // average and derivative of the Phi field
double v=PI;
double mass=1.0;
double lam=mass/v/v; // lambda
double Ed[MESH]; // energy density
double cd[MESH]; // charge density
double memo_n[START/KAKUNIN+1]; // store charge history
double memo_E[START/KAKUNIN+1]; // store energy history

struct{
  double F[MESH+1]; // configuration of Phi field
} od, nd;

// definition of the initial configuration
double func_F(double x); 

// tarapezoidal formula for the multiple integration
double trapezoidal2(void); 

// integration for calculations of the tpological charge
double trapezoidal_N(void); //

// call random number routine
void ransu(void);

int main(void)
{
  double oE, nE; // old energy, new energy     
  double on, nn; // old charge, new charge
  int count = START; // number of calculation
  int tcount = TSTART; // 
  int dcount = 0; // the number of replacement of configurations
  int i; // for loop
  FILE *data, *fp; 
  char name[80]; // file name of field configurations
  char name2[80]; // information of the energy and charge density
  int number=0; // the number of save

  h = (MAX -MIN)/MESH;

  rh=1.0/h;

  // set the initial configuration
  for (i = MESH; i >= 0; i--){
    od.F[i] = func_F(MIN +i*h);
  }
  // boundary condition
  od.F[0] = -v; 
  od.F[MESH/2] = 0; 
  od.F[MESH] = v;

  nd = od;

  // define the paritial derivative of the phi field 
  for(i = MESH -1; i >= 0; i--){
    dx_F[i] = (nd.F[i+1] -nd.F[i])*rh;
  }

  // define the average of the phi field
  for(i = MESH -1; i >= 0; i--){
    tF[i] = (nd.F[i+1] +nd.F[i])*0.5;
  }
  
  // calculate the energy
  nn = trapezoidal_N();
  nE = trapezoidal2();
  oE = nE;
  on = nn;

  // naming for storing temporary results of calculation
  sprintf(name,"kink_%d.txt" ,number);
  sprintf(name2,"kink_%d_Bz.txt" ,number);
  printf(name);

  // for storing of the calculation results
  data = fopen(name,"w");
  for (i = MESH; i >= 0; i--){
    fprintf(data,"%lf\t%lf\n", MIN +i*h, nd.F[i]);
  }
  fclose(data);

  data = fopen(name2,"w");
  for (i = MESH -1; i >= 0; i--){
    fprintf(data,"%lf\t%lf\n", MIN+i*h, Ed[i]);
  }
  fclose(data);

  // change scaler field according to random numbers
  init_genrand((unsigned)time(NULL)); 
  do{

    // for check
    if(count%KAKUNIN == 0){
      memo_E[(START-count)/KAKUNIN] = oE;
      memo_n[(START-count)/KAKUNIN] = on;
      printf("%lf %lf %d\n",oE,on,dcount);
      dcount=0;
    }
    if(count%10000 == 0) printf("\n  %d\n",count);

    // perturb configurations according to random numbers       
    if(count > (START*0.4) ){ double tmp = M1_2*1.2;
      for(i = tmp; i >= 0; i--){
	ransu();
      }
    }
    else if ( (count > (START*0.1) ) && ( (count < (START*0.4) ) ) ) {
      double tmp = M1_2*1.0;
      for (i = tmp; i >= 0; i--){
	ransu();
      }
    }
    else {
      double tmp = M1_2*0.8;
      for (i=tmp; i>=0; i--){
	ransu();
      }
    }

    // boundary condition
    nd.F[0] = -v; 
    nd.F[MESH/2] = 0; 
    nd.F[MESH] = v;

    // define partial derivative of phi field
    for(i=MESH-1; i>=0; i--){
      dx_F[i] = (nd.F[i+1] -nd.F[i])*rh;
    }

    // define average of phi field
    for(i=MESH-1; i>=0; i--){
      tF[i] = (nd.F[i+1] +nd.F[i])*0.5;
    }

    // change of the topological charge
    nn = trapezoidal_N();//b/(4*pi)  
    // calcualte energy with the changed configuration
    nE = trapezoidal2();

    // soften the comparing condition of energy
    if((tcount == 0) && (count > (START*0.1))){
      tcount = TSTART;
      oE*=1.1;
    }
    if(nE < oE){
      oE = nE;
      on = nn;
      tcount = TSTART;
      dcount++;
      od = nd; // update all the configurations here         
    }
    if(oE != nE)  tcount--;
    if(count%5000 == 0){
      number++;
      sprintf(name,"kink_%d.txt" ,number);
      sprintf(name2,"kink_%d_Bz.txt" ,number);

      // for saving calculation results
      data = fopen(name,"w"); 
      for(i=MESH; i>=0; i--){
	fprintf(data,"%lf\t%lf\n", MIN+i*h, nd.F[i]);
      }
      fclose(data);
      data = fopen(name2,"w");
      for(i=MESH-1;i>=0;i--){
	fprintf(data,"%lf\t%lf\t%lf\n", MIN+i*h, Ed[i], cd[i]);
      }
      fclose(data);
    }
    count--;
  }while(count>0);

  // draw oE,on vs. count graph
  data = fopen("kink_c-E.txt","w"); // for saving calculation results 
  for (i=START/KAKUNIN;i>0;i--){
    fprintf(data,"%d\t%lf\t%lf\n", i*KAKUNIN, memo_E[i], memo_n[i]);
  }
  fclose(data);
  return 0;
}


// trapezoidal formula of multiple integral
double trapezoidal2(void) {
  double S = 0.0;
  double s; // coefficient for calculation of boundary points of trapezoidal integral
  
  int i; // loop variable
  for (i=MESH-1;i>=0;i--){
    if ((i==0)||(i==MESH-1)) s = 0.5;
    else s=1.0;
    Ed[i] = (0.5*dx_F[i]*dx_F[i] +lam/4.0*(tF[i]*tF[i] -v*v)*(tF[i]*tF[i] -v*v))*s;
    S += Ed[i]*h;
  }
  return S; 
}

// integral for calculating the topological charge
double trapezoidal_N(void)
{
  double S = 0.0;
  double s; // coefficient for calculation of boundary points of trapezoidal integral
  int i;
  for (i=MESH-1;i>=0;i--){
    if ((i==0)||(i==MESH-1)) s = 0.5;
    else s = 1.0;
    cd[i] = (dx_F[i]/(2*v))*s*h;
    S += cd[i];
  }
  return S;
}

void ransu(void)
{
  int rk,ri,rj; // new random variables
  int random; // for storing random variables
  double leng = 0.001; // initial amplitude parameter

  // change random variables discretely        
  random = genrand_int31()%((MESH-1)*(MESH-1)*2);
  rk = random%2;
  random = (random -rk)*0.5;
  ri = random%(MESH -1);
  ri += 1;
  switch(rk){
  case 0: nd.F[ri]=od.F[ri]+leng; break;
  case 1: nd.F[ri]=od.F[ri]-leng;
  }
}

// define initial configuration
double func_F(double x) {
  return x;
}
