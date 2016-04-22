#include "mesh.hpp"
#include <string>

using namespace Geometry;
using namespace std;

Grid::Grid(ifstream & file){
string str,aux,aux1,aux2;
unsigned int Npoly=0,Nvert=0;

//legge la prima riga
getline(file,str);
stringstream ss(str);
getline(ss,aux,' ');
Nvert=stoi(aux);
getline(ss,aux,' ');
Npoly=stoi(aux);

vect.resize(Nvert);

//leggo le altre righe
for (unsigned int i=0; i<Nvert; i++){
	getline(file,str);
	//stringstream ss(str);
	stringstream ss(str);
	getline(ss,aux,' ');
	getline(ss,aux1,' ');
	getline(ss,aux2,' ');
	vect[i].set(stod(aux1),stod(aux2));
}

//stampa
//for (auto i = vect.begin(); i != vect.end(); i++) cout<<i->x()<<" "<<i->y()<<endl;

//abspol.resize(Npoly);
int tipo;

//leggo le altre righe
for (unsigned int i=0; i<Npoly; i++){
	getline(file,str);
	stringstream ss(str);
	getline(ss,aux,' ');
	getline(ss,aux1,' ');
	tipo=stoi(aux1);
	vector<Point2D> vvv;
	while (getline(ss,aux,' ')){
		vvv.push_back(vect[stoi(aux)]);
	}

if (tipo==0) {Triangle t(vvv); abspol.push_back(make_shared<Triangle> (t));}
if (tipo==1) {Square s(vvv); abspol.push_back(make_shared<Square> (s));}
if (tipo>=2) {Polygon p(vvv); abspol.push_back(make_shared<Polygon> (p));}
}

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
