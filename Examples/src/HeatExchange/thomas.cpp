#include "thomas.hpp"

std::vector<double> thomas(std::vector<double> const & alpha, std::vector<double> const & beta, 
	std::vector<double> const & gamma, std::vector<double> const & f){

	using namespace std;
	const unsigned int M=alpha.size();
	vector<double> y(M,0),x(M+1,0);

	//Risolvo il sistema Ly=f
	y[0]=f[0];
	for (unsigned int i=1; i<M; ++i)
		y[i]=f[i]-beta[i-1]*y[i-1];
	//Risolvo il sistema Ux=y
	x[M]=y[M-1]/alpha[M-1];
	for (unsigned int i=M-1; i>0; --i)
		x[i]=1/alpha[i-1]*(y[i-1]-gamma[i-1]*x[i+1]);

return x;
}