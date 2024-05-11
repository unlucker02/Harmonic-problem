#include <iostream>
#include <vector>
#include <set>
#include <fstream>
#include <functional>
#include <chrono>
#include <iomanip>
#include <cmath>

using std::vector;

struct Mesh
{
   vector<double> xM{ };
   vector<double> yM{ };
   vector<double> zM{ };
   vector<vector<int>> areasMesh{ };
};

struct SplittingMesh
{
   vector<double> x{ };
   vector<double> y{ };
   vector<double> z{ };
};

struct ShiftsArrays
{
   vector<int> Ix{ };
   vector<int> Iy{ };
   vector<int> Iz{ };
};

struct BoundaryConds
{
   vector<vector<int>> bC{ };
};

struct SparseMatrix
{
   vector<int> ig{ }, jg{ };
   vector<double> ggl{ };
   vector<double> ggu{ };
   vector<double> di{ };
};

struct ProfileMatrix
{
   vector<int> ig{ };
   vector<double> ggl{ };
   vector<double> ggu{ };
   vector<double> di{ };

   void calcLDU();
   void calcZ(vector<double> &z, vector<double> &b);
   void calcY(vector<double> &y);
   void calcX(vector<double> &x);
   void solveSlae(vector<double> &z, vector<double> &b);
};

struct SLAE
{
   SparseMatrix A{ };
   vector<double> q;
   vector<double> b;
};

struct QuadratureNode
{
   double x = 0;
   double y = 0;
   double z = 0;
   double weight = 0;
};

struct LOS
{
   vector <double> r1{ }, rk{ }, z1{ }, p1{ }, Ar{ }, p{ }, mult{ };
};