#ifndef HH_RKSYSTEM_HH
#define HH_RKSYSTEM_HH
#include "nonLinSys.hpp"
#include <functional>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <iostream>

using namespace NonLinearSystems;


argumentType convert(std::vector<double> & y){
argumentType Y(y.size());
for (unsigned int i=0; i<y.size(); i++) Y[i]=y[i];
 return Y;
};

std::vector<double> convert2(argumentType & y){
std::vector<double> Y(y.size());
for (unsigned int i=0; i<y.size(); i++) Y[i]=y[i];
 return Y;
};


namespace ODE
{

//constexpr std::size_t MAX_STEPS=10000; 
using namespace NonLinearSystems;

template<class prec>
class rksys {
private:
	std::vector<std::pair<prec,prec>> solution;
    //f is a vector of functions (the time is one of the elements of the vector)
    std::vector<std::function<prec (argumentType const & Y)> > dy;
	prec t0;
	prec T;
    //the initial data is also a vector
	std::vector<prec> y0; 
	prec h_initial; 
	prec h_max;
	prec final_error;
	int status;
	std::size_t maxSteps;
	std::vector<prec> rksys_step(std::vector<prec> const & y0, prec const & t0, prec const & h, prec & error);

	std::vector<prec> b1,b2,c;
	std::vector<prec> a;

public:

	rksys<prec> (std::vector<std::function<prec (argumentType const & Y)>> const & dy,prec const & t0,prec const & T,
		std::vector<prec> const & y0,prec const & h_initial,prec const & h_max,prec const & final_error,int & status,std::size_t const & maxSteps=MAX_STEPS):
		dy(dy), t0(t0), T(T), y0(y0), h_initial(h_initial), h_max(h_max), final_error(final_error), maxSteps(maxSteps)
		{};
	//rk45(dy): dy(dy) {};
	std::vector<std::pair<prec,prec>> getSolution(){return solution;};
	void computeSolution();
	void setButcher45();
	void setButcher23();
};

template<class prec>
void rksys<prec>::computeSolution(){ 
    status=0;
    const std::size_t maxReduction=maxSteps;
    // parameters for decreasing/increasing time step
    double const c1=1.0;
    // I need to have a sufficient decrease of the local error
    // to allow time step coarsening
    double const c2=1./64.;

    double length=T-t0;
    //! Make sure that h allows to reach T
    std::size_t initialNSteps=std::max(static_cast<size_t>(1),static_cast<size_t>(length/h_initial));
    double h=length/initialNSteps;
    // To avoid underflow we need in any case to limit the time step to a positive number
    // Here I allow h to become 128 time smaller than that giving the maximal number of steps
    double h_min = length/(128*maxSteps);
    // SOme counters
    std::size_t stepsCounter(0);
    // Initial data
    double time(t0);
    std::vector<prec> y(y0);
    double errorPerTimeStep=final_error/initialNSteps;
    if (initialNSteps>=maxSteps) throw std::runtime_error("RK: initial time step h too small!");
    //std::vector<std::pair<double,double>> solution;
    solution.emplace_back(std::make_pair(t0,y0[0]));
    double localError;
    std::vector<prec> newy;
    while (time<T && stepsCounter <maxSteps)
      {
	//Do a step
	//adjust h if needed for the last step
	if (time + h > T) h = T-time;
	newy = rksys_step(y,time,h,localError);
	while (h> h_min && localError > c1*errorPerTimeStep)
	  {
	    // half time step
	    h /=2;
	    errorPerTimeStep /=2;
	    newy = rksys_step(y,time,h,localError);
	  }
	if (localError>errorPerTimeStep)status=1;
	//! advance
	y = newy;
	time +=h;
	++stepsCounter;
	solution.emplace_back(std::make_pair(time,y[0]));
	//! check if we reached end
	if(localError<c2*errorPerTimeStep && h<h_max)
	  {
	    // Double step
	    h *=2;
	    errorPerTimeStep *=2;
	  }
      }
    //handle exceptions
    if(stepsCounter>=maxSteps && time < T)
      {
	status=2;
	throw std::runtime_error("RKsys: Max number of time steps exceeded");
      }
};


template<class prec>
void rksys<prec>::setButcher45(){
a.resize(15); b1.resize(6); b2.resize(6); c.resize(6);
a[0]=1./4.; a[1]=3./32.; a[2]=9./32.; a[3]=1932./2197.; a[4]=-7200./2197.; a[5]=7296./2197.;
a[6]=439./216.; a[7]=-8.; a[8]=3680./513.; a[9]=-845./4104.; 
a[10]=-8./27.; a[11]=2.; a[12]=-3544./2565.; a[13]=1859./4104.; a[14]=-11./40;

c[0]=0.; c[1]=a[0]; c[2]=a[1]+a[2]; c[3]=a[3]+a[4]+a[5]; c[4]=a[6]+a[7]+a[8]+a[9];
c[5]=a[10]+a[11]+a[12]+a[13]+a[14];

b1[0]=25./216.; b1[1]=0.; b1[2]=1408./2565.; b1[3]=2197./4104.; b1[4]=-1./5.; b1[5]=0.;
b2[0]=16./135.; b2[1]=0.; b2[2]=6656./12825.; b2[3]=28561./56430.; b2[4]=-9./50.; b2[5]=2./55;
};

template<class prec>
void rksys<prec>::setButcher23(){
a.resize(6); b1.resize(4); b2.resize(4); c.resize(4);
a[0]=1./2.; a[1]=0.0; a[2]=3./4.; a[3]=2./9.; a[4]=1./3.; a[5]=4./9.;

c[0]=0.; c[1]=a[0]; c[2]=a[1]+a[2]; c[3]=a[3]+a[4]+a[5];

b1[0]=2./9.; b1[1]=1./3.; b1[2]=4./9.; b1[3]=0.;
b2[0]=7./24.; b2[1]=1./4.; b2[2]=1./3.; b2[3]=1./8.;
};



template<class prec>
std::vector<prec> rksys<prec>::rksys_step(std::vector<prec> const & y0, prec const & t0, prec const & h, prec & error)
  {
    using namespace std;
    using namespace NonLinearSystems;
    NonLinSys F;
    setButcher45();
    //setButcher23();
    for (unsigned int i=0; i<dy.size(); i++) F.addToSystem(dy[i]);


    std::vector<prec> YY=y0; 
    argumentType Y0=convert(YY);
    YY.insert(YY.begin(),t0);
    argumentType Y=convert(YY);

    std::vector<argumentType> K(b1.size());
    K[0]=F(Y);

    auto it1=a.begin(),it2=it1;

    for (std::size_t i=1; i<b1.size(); i++){
    	it2=it1+i;
    	//K[i]=f(t0+c[i]*h,y0+h*std::inner_product(it1,it2,K.begin(),0.0));
        argumentType tmp(Y.size()); tmp[0]=t0+c[i]*h;
        argumentType NullV(Y0.size());
        argumentType aux=Y0+h*std::inner_product(it1,it2,K.begin(),NullV);
        for (unsigned int j=1; j<tmp.size(); j++) tmp[j]=aux[j-1];
        K[i]=F(tmp);
    	it1=it2;
    }
    
    argumentType Y1=Y0+h*std::inner_product(b1.begin(),b1.end(),K.begin(),argumentType(Y.size()-1));
    argumentType Y2=Y0+h*std::inner_product(b2.begin(),b2.end(),K.begin(),argumentType(Y.size()-1));
    error=(Y2-Y1).norm();
    //cout<<error<<endl;
    return convert2(Y2);
    
  };



  }// end namespace
#endif