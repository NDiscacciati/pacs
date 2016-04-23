#ifndef HH_EDGE_HH
#define HH_EDGE_HH
#include "Polygon.hpp"
#include <iostream>
#include <vector>
#include <array>
#include <utility>
#include <memory>
#include <fstream>
#include <sstream>
#include <iterator>

using namespace Geometry;
using namespace std;

class Edge {
public:
	//Edge ():a(0),b(0){};
	Edge()=default;
    Edge (const Edge & )=default;
    Edge & operator= (const Edge & )= default;
    Edge(unsigned int aa, unsigned int bb): a(aa),b(bb){};

    void set(unsigned int const &aa, unsigned int const &bb){
      a=aa; b=bb;};
    unsigned int x(){return a;};
    unsigned int y(){return b;};
    void print(){cout<<a<<" "<<b<<endl;};
    friend bool operator <(Edge const &f, Edge const &s)
    {
    	unsigned int m1,m2,M1,M2;
    	m1=min(f.a,f.b); m2=min(s.a,s.b);
    	M1=max(f.a,f.b); M2=max(s.a,s.b);
    	if (m1==m2) return M1<M2;
    	return m1<m2;
    };
    //friend ostream & operator << (ofstream & ost, Edge const & e){
    //	ost<<e.a<<" "<<e.b<<endl;
    //	return ost;
    //};



private:
	unsigned int a;
	unsigned int b;
};


#endif