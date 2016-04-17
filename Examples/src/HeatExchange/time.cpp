#include "time.hpp"

std::vector<double> time(double const dt, double const T, int const M, double const act, double const init){
	using namespace std;
	double h=1./M;

	//I add one line to the matrix in order to handle Dirichlet condition in a safer way
	vector<double> a(M+1,1+dt/(h*h)*(2+h*h*act)),b(M,-dt/(h*h)*1.),c(M,-dt/(h*h)*1.);
	b[M-1]=-1.; a[M]=1.; a[0]=1.; c[0]=0.; //modification of terms
	vector<double> alpha(M+1,0.),beta(M,0.),gamma(M,0.);
	//vector<double> rhs(M+1);
	vector<double> theta(M+1,init);
	unsigned int steps=round(T/dt);

cout<<"dt= "<<dt<<" , T= "<<T<<" , steps= "<<steps<<endl;

//factorization for Thomas algorithm (constant in time)
 alpha[0]=a[0];
 for (int i=1; i<M+1; ++i){
 beta[i-1]=b[i-1]/alpha[i-1];
 alpha[i]=a[i]-beta[i-1]*c[i-1];
 gamma[i-1]=c[i-1];
  }


for (unsigned int t=1; t<steps+1; ++t){

	//Right hand side
	theta[M]=0; //to impose Neumann condition 
	//if I don't put the above condition, I don't have to change the last line in the matrix

	//Direct solution of linear system
	theta=thomas(alpha,beta,gamma,theta);
	theta[0]=init; //Dirichlet condition (optional: comes from the solver, but I can have approximation errors)

	//Print on screen solution at each step
	for (auto m: theta) cout<<m<<"\t";
	cout<<endl<<endl;
}

return theta;
}

