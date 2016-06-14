#include "rk45.hpp"
#include "rkimp.hpp"
#include "rksys.hpp"
#include <iostream>
#include <fstream>

#include <cmath>
int main()
{
  using namespace std;
  using namespace ODE;
  using namespace NonLinearSystems;
  // auto fun = [](double const & t, double const & y){return -10*y;};
  // auto fun = [](double const & t, double const & y){return -std::sin(t);};
   auto fun=[](double const &t, double const & y){return 1.0;};
  double t0=0;
  double y0=1;
  double T=100;
  double h_init=0.2;
  double errorDesired=1.e-7;
  int status;


//part 1
  //Template version
  rk<double> tmp(fun,t0,T,y0,h_init,(T-t0)/4.,errorDesired,status,10000);
  tmp.computeSolution();
  auto result=tmp.getSolution();

  ofstream file("result.dat");
  std::cout<<"Writing on file basic version"<<std::endl;
  for (auto v : result)
  file<<v.first<<" "<<v.second<<std::endl;
  file.close();
//end of part 1

 /* 
 //part 2
  //Diagonally implicit
  rkImp<double> tmp2(fun,t0,T,y0,h_init,(T-t0)/4.,errorDesired,status,10000);
  tmp2.computeSolution();
  auto result2=tmp2.getSolution();

  ofstream file("result.dat");
  std::cout<<"Writing on file diagonally implicit"<<std::endl;
  for (auto v : result2)
  file<<v.first<<" "<<v.second<<std::endl;
  file.close();
//end of part 2
  */


/*
//part 3
  //System
  //auto f1=[] (argumentType const & y){return y[1];};
  auto f1=[] (argumentType const & y){return -y[1];};
  auto f2=[] (argumentType const & y){return -y[2];};
  std::vector<std::function<double (argumentType const & y)> > Y;
  Y.push_back(f1); Y.push_back(f2);
  std::vector<double> Y0(2,1.0);
  rksys<double> sys(Y,t0,T,Y0,h_init,(T-t0)/4.,errorDesired,status,10000);
  sys.computeSolution();
  auto resultsys=sys.getSolution();

  ofstream file("result.dat");
  std::cout<<"Writing on file for systems"<<std::endl;
  for (auto v : resultsys)
    file<<v.first<<" "<<v.second<<std::endl;
  file.close();
//end of part 3
*/


/*
  // Only if I know the exact solution
  //auto exact=[](double const &t){return std::exp(-10.*t);}
  auto exact=[](double const &t){return std::cos(t);};
  double max_error(0);
  for (auto i : result2)
    max_error=std::max(max_error,std::abs(exact(i.first)-i.second));
  std::cout.setf(std::ios::scientific);
  std::cout<<"Max error "<<max_error<<" Desired max error "<<errorDesired;
  std::cout<<std::endl;
  */
}
