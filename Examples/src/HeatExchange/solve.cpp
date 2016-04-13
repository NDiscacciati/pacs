#include "solve.hpp"

int solve(const int M, const double act, const double toler, const int itermax, std::vector<double> &theta, const int norm){
	
	using namespace std;
  const auto h=1./M;

  int iter=0, status=0;
  double xnew, epsilon;
  vector<double> diff(M+1,0.);
     do
     { epsilon=0.;

	 //Compute the solution
	 for (int m=1; m<M; m++){
	 	xnew = (theta[m-1]+theta[m+1])/(2.+h*h*act);
	 	diff[m]=xnew-theta[m];
	 	theta[m]=xnew;
	 }
	 xnew=theta[M-1];
	 diff[M]=xnew-theta[M];
	 theta[M]=xnew;

	 //Compute the Rn norm
	 //for (int m=1; m<M+1; m++) epsilon+=(diff[m]*diff[m]);
	 
	 //Compute the L2 norm
	 epsilon+=h/3*(diff[1]*diff[1]);
	 for (int m=2; m<M+1; m++)
	 	//epsilon+=h/3*(diff[m-1]*diff[m-1]+diff[m]*diff[m]+diff[m-1]*diff[m]);
	 	epsilon+=h/6*(diff[m-1]*diff[m-1]+diff[m]*diff[m]+(diff[m-1]+diff[m])*(diff[m-1]+diff[m]));
	

	 //Compute the H1 norm, if necessary
	 if (norm==1){
	 	epsilon+=1/h*(diff[1]*diff[1]);
	 	for (int m=2; m<M+1; m++)
	 		epsilon+=1/h*(diff[m]-diff[m-1])*(diff[m]-diff[m-1]);
	 }

	 iter++;     
       }while((sqrt(epsilon) > toler) && (iter < itermax) );


    cout<<"The final "<<(norm==0 ? "L2": "H1")<<" norm is "<<sqrt(epsilon)<<endl;
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