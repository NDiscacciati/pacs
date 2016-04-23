#include "grid.hpp"
#include <string>

using namespace Geometry;
using namespace std;

Grid::Grid(ifstream & file){
string str,aux,aux1,aux2;
unsigned int Npoly=0,Nvert=0,tipo=0;
set<Edge> insieme, interni;
Edge e;
double i1,i2;
vector<unsigned int> lati;
vector<Point2D> vvv;

//read the first line
getline(file,str);
stringstream ss(str);
/*
getline(ss,aux,' ');
Nvert=stoi(aux);
getline(ss,aux,' ');
Npoly=stoi(aux);
*/
ss>>Nvert>>Npoly;
cout<<Nvert<<"   "<<Npoly<<endl;

vect.resize(Nvert);

//leggo le altre righe: coordinate dei vertici
for (unsigned int i=0; i<Nvert; i++){
	getline(file,str);
	stringstream ss(str);
	ss>>i1>>i1>>i2;
	//cout<<i1<<"space"<<i2<<endl;
	//vect[i].set(stod(aux1),stod(aux2));
	vect[i].set(i1,i2);
}

//stampa
//for (auto i = vect.begin(); i != vect.end(); i++) cout<<i->x()<<" "<<i->y()<<endl;
pair<set<Edge>::iterator,bool> check;


//leggo le altre righe: matrice di connettività
for (unsigned int i=0; i<Npoly; i++){
	getline(file,str);
	stringstream ss(str);
	getline(ss,aux,' ');
	getline(ss,aux,' ');
	tipo=stoi(aux);
	vvv.clear(); lati.clear();
	while (getline(ss,aux,' ')){
		vvv.emplace_back(vect[stoi(aux)]);
		lati.emplace_back(stoi(aux));
		cout<<stoi(aux)<<endl;
	}

	//create a set containing the edges and an auxiliary one to store internal edges
	for(unsigned int j=0; j<lati.size()-1; j++){
		e.set(lati[j],lati[j+1]);
		check=insieme.insert(e);
		if (!check.second) interni.insert(e);
	}
	e.set(lati[lati.size()-1],lati[0]);
	check=insieme.insert(e);
	if (!check.second) interni.insert(e);

	//build the vector
	switch(tipo){
		case(0): {Triangle t(vvv); abspol.emplace_back(make_shared<Triangle> (t));} break;
		case(1): {Square s(vvv); abspol.emplace_back(make_shared<Square> (s));} break;
		default: {Polygon p(vvv); abspol.emplace_back(make_shared<Polygon> (p));}
	}
}

//copy all edges in the vector
//copy(insieme.begin(),insieme.end(),back_inserter(AllEdges));
AllEdges.resize(insieme.size());
copy(insieme.begin(),insieme.end(),AllEdges.begin());

//create the vector containing Boundary edges
set_difference(insieme.begin(), insieme.end(),interni.begin(),interni.end(),back_inserter(Boundary));
}


double Grid::area(){
	unsigned int Npoly=abspol.size();
	double area{0.0};
	for (unsigned int i=0; i<Npoly; i++){
		//cout<<"Area del poligono "<<i<<": "<<abspol[i]->area()<<endl;
		area+=abspol[i]->area();
	}
return area;
}

void Grid::printedges(){
	cout<<"All Edges:"<<endl;
	for(auto i=AllEdges.begin(); i!=AllEdges.end(); i++)
		i->print();
	cout<<"Lati di bordo:"<<endl;
	for(auto i=Boundary.begin(); i!=Boundary.end(); i++)
		i->print();
	
	vector<Edge> aux;
	set_difference(AllEdges.begin(),AllEdges.end(),Boundary.begin(),Boundary.end(),back_inserter(aux));
	cout<<"Lati interni:"<<endl;
	for(auto i=aux.begin(); i!=aux.end(); i++)
		i->print();
	
	//ostream_iterator<Edge> ost(cout,"");
	//set_difference(AllEdges.begin(),AllEdges.end(),Boundary.begin(),Boundary.end(),ost);
}