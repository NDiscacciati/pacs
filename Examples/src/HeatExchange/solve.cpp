#include "solve.hpp"

int solve(const int M, const double act, const double toler, const int itermax, std::vector<double> &theta, const int norm){
	
	using namespace std;
  const auto h=1./M;

  int iter=0, status=0,aux=(norm==1 ? 0 : 1);
  double diff=0.,xnew, epsilon;
     do
       { epsilon=0.;

     //First row
	   xnew  = (theta[0]+theta[2])/(2.+h*h*act);
	   epsilon+=h/3*(xnew-theta[1])*(xnew-theta[1])+aux/h*((xnew-theta[1])*(xnew-theta[1]));
	   diff= xnew-theta[1];
	   theta[1] = xnew;

	 //Central elements
         for(int m=2;m < M;m++)
         {   
	   xnew  = (theta[m-1]+theta[m+1])/(2.+h*h*act);
	   epsilon += h/3*(diff*diff+(xnew-theta[m])*(xnew-theta[m])+diff*(xnew-theta[m]))
	   			+aux/h*((xnew-theta[m]-diff)*(xnew-theta[m]-diff));
	   diff=xnew-theta[m];
	   theta[m] = xnew;
         }

	 //Last row
	 xnew = theta[M-1]; 
	 epsilon += h/3*(diff*diff+(xnew-theta[M])*(xnew-theta[M])+diff*(xnew-theta[M]))
	 			+aux/h*((xnew-theta[M]-diff)*(xnew-theta[M]-diff));
	 theta[M]=  xnew; 

	 iter++;     
       }while((sqrt(epsilon) > toler) && (iter < itermax) );


    cout<<"The final "<<(norm==1 ? "L2": "H1")<<" norm is "<<sqrt(epsilon)<<endl;
    if(iter<itermax)
      cout << "M="<<M<<"  Convergence in "<<iter<<" iterations"<<endl;
    else
      {
	cerr << "NOT CONVERGING in "<<itermax<<" iterations "<<
	  "||dx||="<<sqrt(epsilon)<<endl;
	status=1;
      }

	return status;
}