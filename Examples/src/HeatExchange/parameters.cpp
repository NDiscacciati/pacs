#include <iostream>
#include "parameters.hpp"
std::ostream & operator << (std::ostream & out,const parameters & p)
{
  out<<"PARAMETER VALUES:"<<"\n";
  out<<"itermax= "<<p.itermax<<"\n";
  out<<"toler= "<<p.toler<<"\n";
  out<<"L= "<<p.L<<"\n";
  out<<"a1= "<<p.a1<<"\n";
  out<<"a2= "<<p.a2<<"\n";
  out<<"To= "<<p.To<<"\n";
  out<<"Te= "<<p.Te<<"\n";
  out<<"k= "<<p.k<<"\n";
  out<<"hc= "<<p.hc<<"\n";
  out<<"M= "<<p.M<<"\n\n";
  out<<"Output file= "<<p.name<<"\n";
  out<<"Output mode= "<<p.see<<"\n";
  out<<"Stopping criterion "<<(p.norm==1 ? "L2": "H1")<<" norm"<<"\n";
  return out;
}
