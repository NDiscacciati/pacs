#ifndef HH_RKIMP_HH
#define HH_RKIMP_HH
#include <functional>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <numeric>
#include <iostream>
#include "nonLinSys.hpp"
#include "newton.hpp"

namespace ODE
{ 


template<class prec>
class rkImp {
private:
	std::vector<std::pair<prec,prec>> solution;
    std::function<prec (prec const &, prec const &)> dy;
	prec t0;
	prec T;
	prec y0; 
	prec h_initial; 
	prec h_max;
	prec final_error;
	int status;
	std::size_t maxSteps;
	prec rkImp_step(prec const & y0, prec const & t0, prec const & h, prec & error);

	std::vector<prec> b1,b2,c;
	std::vector<prec> a;

public:

	rkImp<prec> (std::function<prec (prec const &, prec const &)> const & dy,prec const & t0,prec const & T,
		prec const & y0,prec const & h_initial,prec const & h_max,prec const & final_error,int & status,std::size_t const & maxSteps=MAX_STEPS):
		dy(dy), t0(t0), T(T), y0(y0), h_initial(h_initial), h_max(h_max), final_error(final_error), maxSteps(maxSteps)
		{};
	//rk45(dy): dy(dy) {};
	std::vector<std::pair<prec,prec>> getSolution(){return solution;};
	void computeSolution();
	void setButcher2();
};

template<class prec>
void rkImp<prec>::computeSolution(){ 
    status=0;
    const std::size_t maxReduction=maxSteps;
    // parameters for decreasing/increasing time step
    double const c1=1.0;
    // I need to have a sufficient decrease of the local error
    // to allow time step coarsening
    double const c2=1./64.;

    prec length=T-t0;
    //! Make sure that h allows to reach T
    std::size_t initialNSteps=std::max(static_cast<size_t>(1),static_cast<size_t>(length/h_initial));
    prec h=length/initialNSteps;
    // To avoid underflow we need in any case to limit the time step to a positive number
    // Here I allow h to become 128 time smaller than that giving the maximal number of steps
    prec h_min = length/(128*maxSteps);
    // SOme counters
    std::size_t stepsCounter(0);
    // Initial data
    prec time(t0);
    prec y(y0);
    prec errorPerTimeStep=final_error/initialNSteps;
    if (initialNSteps>=maxSteps) throw std::runtime_error("RK: initial time step h too small!");
    //std::vector<std::pair<double,double>> solution;
    solution.emplace_back(std::make_pair(t0,y0));
    prec localError;
    prec newy;
    while (time<T && stepsCounter <maxSteps)
      {
	//Do a step
	//adjust h if needed for the last step
	if (time + h > T) h = T-time;
	newy = rkImp_step(y,time,h,localError);
	while (h> h_min && localError > c1*errorPerTimeStep)
	  {
	    // half time step
	    h /=2;
	    errorPerTimeStep /=2;
	    newy = rkImp_step(y,time,h,localError);
	  }
	if (localError>errorPerTimeStep)status=1;
	//! advance
	y = newy;
	time +=h;
	++stepsCounter;
	solution.emplace_back(std::make_pair(time,y));
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
	throw std::runtime_error("RK: Max number of time steps exceeded");
      }
};


template<class prec>
void rkImp<prec>::setButcher2(){
a.resize(3); b1.resize(2); b2.resize(2); c.resize(2);
a[0]=0.0; a[1]=1./2.; a[2]=1./2.;

c[0]=0.; c[1]=a[1]+a[2];

b1[0]=1./2; b1[1]=1./2;
b2[0]=1.; b2[1]=0.;
};



template<class prec>
prec rkImp<prec>::rkImp_step(prec const & y0, prec const & t0, prec const & h, prec & error)
  {
    using namespace std;
    using namespace NonLinearSystems;
    auto f=dy;
    setButcher2();

    std::vector<prec> K(b1.size());
    K[0]=f(t0,y0);

    auto it1=a.begin(),it2=it1;



    for (std::size_t i=1; i<b1.size(); i++){
    	it2=it1+i;
        NonLinSys F;
        F.addToSystem(
            [=](argumentType const & KK)->prec
            {return KK[0]-f(t0+c[i]*h,y0+h*std::inner_product(it1,it2-1,K.begin(),0.0)+h*(*(it2-1))*KK[0]);}
            );
        DiscreteJacobian J(&F);
        NewtonOptions optNewton{1e-6,1.e-8,100};
        argumentType yy(1);
        yy[0]=y0;
        NewtonStatus result=Newton(F,J,yy,optNewton);
    	K[i]=yy[0];
    	it1=it2;
    }
    prec y1=y0+h*std::inner_product(b1.begin(),b1.end(),K.begin(),0.0);
    prec y2=y0+h*std::inner_product(b2.begin(),b2.end(),K.begin(),0.0);
    error=std::abs(y2-y1);
    //cout<<error<<endl;
    return y2;
    
  };

  }// end namespace
#endif