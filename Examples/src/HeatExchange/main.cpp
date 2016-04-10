#include <iostream> // input output
#include <cmath> // (for sqrt)
#include <vector>
#include <tuple>
#include <string>
#include "readParameters.hpp"
#include "GetPot.hpp"
#include "solve.hpp"
#include "thomas.hpp"
//#include "gnuplot-iostream.hpp"// interface with gnuplot
/*!
  @file main.cpp
  @brief Temperature distribution in a 1D bar.

  @detail
    We solve  \f$ -T^{\prime\prime}(x)+act*(T(x)-T_e)=0, 0<x<L \f$ with 
    boundary conditions \f$ T(0)=To; T^\prime(L)=0\f$
    
    **************************************************
    Linear finite elements
    Iterative resolution by Gauss Siedel.
    **************************************************
    
    Example adapted by Luca Formaggia from  a code found in 
    "Simulation numerique an C++" di I. Danaila, F. Hecht e
    O. Pironneau.
*/
//! helper function
void printHelp()
{
  std::cout<<"USAGE: main [-h] [-v] -p parameterFile (default: parameters.pot)"<<std::endl;
  std::cout<<"-h this help"<<std::endl;
  std::cout<<"-v verbose output"<<std::endl;
}

//! main program
int main(int argc, char** argv)
{
  using namespace std; // avoid std::
  int status(0); // final program status
  GetPot   cl(argc, argv);
  if( cl.search(2, "-h", "--help") )
    {
      printHelp();
      return 0;
    }
  // check if we want verbosity
  bool verbose=cl.search(1,"-v");
  // Get file with parameter values
  string filename = cl.follow("parameters.pot","-p");
  cout<<"Reading parameters from "<<filename<<std::endl;
  // read parameters
  const parameters param=readParameters(filename,verbose);
  // Transfer parameters to local variables
  // I use references to save memory (not really an issue here, it is just
  // to show a possible  use of references)
  const int&    itermax= param.itermax;   //max number of iteration for Gauss-Siedel
  const double& toler=param.toler;   // Tolerance for stopping criterion
  // Here I use auto (remember that you need const and & if you want constant references)
  const auto& L= param.L;  // Bar length
  const auto& a1=param.a1; // First longitudinal dimension
  const auto& a2=param.a2; //  Second longitudinal dimension
  const auto& To=param.To; // Dirichlet condition
  const auto& Te=param.Te; // External temperature (Centigrades)
  const auto& k=param.k;  // Thermal conductivity
  const auto& hc=param.hc; // Convection coefficient
  const auto& M=param.M; // Number of grid elements
  const auto& name=param.name; //Name of output file
  const int& see=param.see; //See the results on screen or file
  const int& norm=param.norm; //Norm parameter
  
  //! Precomputed coefficient for adimensional form of equation
  const auto act=2.*(a1+a2)*hc*L*L/(k*a1*a2);

  // mesh size
  const auto h=1./M;
  
  // Solution vector
  std::vector<double> theta(M+1);


  /*
  //Challenge 1.2 (solve the system with different norms)
  //initialization for GS
  for(unsigned int m=0;m <= M;++m)
   theta[m]=(1.-m*h)*(To-Te)/Te;

  //solve the system with GS
  status=solve(M,act,toler,itermax,theta,norm); 
  */

  //Challenge 1.3 (Thomas Algorithm)
  //construction of the matrix
  vector<double> a(M,2+h*h*act),b(M-1,-1.),c(M-1,-1.);
  vector<double> alpha(M),beta(M-1),gamma(M-1);
  vector<double> f(M); 
  a[M-1]=1.;
  f[0]=(To-Te)/Te;

  //Factorization
  alpha[0]=a[0];
  for (unsigned int i=1; i<M; ++i){
  	beta[i-1]=b[i-1]/alpha[i-1];
  	alpha[i]=a[i]-beta[i-1]*c[i-1];
  	gamma[i-1]=c[i-1];
  }
  //solution of the system
  theta=thomas(alpha,beta,gamma,f);
  theta[0]=(To-Te)/Te;

 // Analitic solution

    vector<double> thetaa(M+1);
     for(int m=0;m <= M;m++)
       thetaa[m]=Te+(To-Te)*cosh(sqrt(act)*(1-m*h))/cosh(sqrt(act));

     // writing results with format
     // x_i u_h(x_i) u(x_i) and lauch gnuplot 

     //Gnuplot gp;
     //std::vector<double> coor(M+1);
     //std::vector<double> sol(M+1);
     //std::vector<double> exact(M+1);

    cout<<"You have requested to print the results in mode "<<see<<endl;
    ofstream file; 
    int mysee=see; //auxiliary variable that I can modify
    if (see%2==1) {file.open(name);  cout<<"Result file: "<<name<<endl;}
    //ostream &f = cout;
    do{
     //ofstream f(name);
      ostream &f = (mysee%2==1 ? file : cout);
     for(int m = 0; m<= M; m++)
       {
	 // \t writes a tab 
         f<<m*h*L<<"\t"<<Te*(1.+theta[m])<<"\t"<<thetaa[m]<<endl;
	 // An example of use of tie and tuples!
         
	 //std::tie(coor[m],sol[m],exact[m])=
	   //std::make_tuple(m*h*L,Te*(1.+theta[m]),thetaa[m]);
       }
       mysee++;
    } while (mysee==4);
       /*
     // Using temporary files (another nice use of tie)
     gp<<"plot"<<gp.file1d(std::tie(coor,sol))<<
       "w lp title 'uh',"<< gp.file1d(std::tie(coor,exact))<<
       "w l title 'uex'"<<std::endl;
       */
     if (see%2==1) file.close();
     return status;
}
