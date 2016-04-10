#include "time.hpp"

std::vector<double> time(){
	using namespace std;
	double dt=0.1,h=0.1, T=100,act=1;
	int M=100;
	vector<double> Md(M,h*h/dt*2/3),Mu(M-1,h*h/dt/6),Ml(M-1,h*h/dt/6);
	//vector<double> MAd(M,2+(act*h*h+h*h/dt)*2/3),MAu(M-1,(act*h*h+h*h/dt)*1/6),MAl(M-1,(act*h*h+h*h/dt)*1/6);
	vector<double> a(M,2+(act*h*h+h*h/dt)*2/3),c(M-1,(act*h*h+h*h/dt)*1/6),b(M-1,(act*h*h+h*h/dt)*1/6);
	vector<double> alpha(M),beta(M),gamma(M);
	vector<double> rhs(M);
	vector<double> theta(M+1);

//fattorizzazione alla Thomas
 alpha[0]=a[0];
 for (unsigned int i=1; i<M; ++i){
 beta[i-1]=b[i-1]/alpha[i-1];
 alpha[i]=a[i]-beta[i-1]*c[i-1];
 gamma[i-1]=c[i-1];
  }

for (double t=0; t<T; t+=dt){

	//termine noto a destra
	rhs[0]=Md[0]*theta[1]+Md[1]*theta[2]+theta[0];
	for (unsigned int i=1; i<M-1; i++)
		rhs[i]=Ml[i-1]*theta[i]+Md[i]*theta[i+1]+Mu[i+1]*theta[i+2];
	rhs[M-1]=Ml[M-2]*theta[M-1]+Md[M-1]*theta[M];

	//risoluzione del sistema lineare
	theta=thomas(alpha,beta,gamma,rhs);
}

return theta;
}

