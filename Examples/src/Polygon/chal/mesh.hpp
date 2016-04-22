#ifndef HH_GRID_HH
#define HH_GRID_HH
#include "Polygon.hpp"
#include <iostream>
#include <vector>
#include <array>
#include <utility>
#include <memory>
#include <fstream>
#include <sstream>

using namespace Geometry;
using namespace std;

class Grid
  {
  public:
    Grid ()=default;
    Grid (const Grid & )=default;
    Grid & operator= (Grid & )= default;

    Grid (ifstream & file);
  private:
    std::vector<Point2D> vect;
    std::vector<std::shared_ptr<AbstractPolygon> > abspol;
  };

  #endif