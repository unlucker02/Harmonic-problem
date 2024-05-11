#include "functions.h"

void getQuadratures3D(vector<QuadratureNode> &quadratures)
{
//   vector<double> quadratures1D = { -0.77459666924148, 0, 0.77459666924148 };
//   vector<double> weights1D = { 0.55555555555556, 0.88888888888889, 0.55555555555556 };

   double p1 = sqrt((3. - 2. * sqrt(1.2)) / 7.);
   double p2 = sqrt((3. + 2. * sqrt(1.2)) / 7.);

   double w1 = (18. + sqrt(30)) / 36.;
   double w2 = (18. - sqrt(30)) / 36.;

   vector<double> quadratures1D = { -p1, -p2, p1, p2 };
   vector<double> weights1D = { w1, w2, w1, w2 };


   QuadratureNode node{ };

   for (int r = 0; r < 4; r++)
      for (int s = 0; s < 4; s++)
         for (int p = 0; p < 4; p++)
         {
            node.x = quadratures1D[p];
            node.y = quadratures1D[s];
            node.z = quadratures1D[r];
            node.weight = weights1D[p] * weights1D[s] * weights1D[r];

            quadratures.push_back(node);
         }
}

double normL2(ShiftsArrays &I, vector<vector<int>> &areasMesh, ParametresMesh &coefs, vector<double> &q, SplittingMesh &sMesh)
{
   auto &xCoords = sMesh.x;
   auto &yCoords = sMesh.y;
   auto &zCoords = sMesh.z;

   const int xSize = xCoords.size();
   const int ySize = yCoords.size();
   const int zSize = zCoords.size();

   vector<QuadratureNode> nodes{ };
   getQuadratures3D(nodes);

   double normL2 = 0;
   int l = 0;

   for (int r = 0; r < zSize - 1; r++)
   {
      for (int s = 0; s < ySize - 1; s++)
      {
         for (int p = 0; p < xSize - 1; p++)
         {
            if (numberOfArea(I, areasMesh, p, s, r, l) != -1)
            {
               double hx = xCoords[p + 1] - xCoords[p];
               double hy = yCoords[s + 1] - yCoords[s];
               double hz = zCoords[r + 1] - zCoords[r];

               double valueOfIntegral = 0;

               for (auto &node : nodes)
               {
                  double x = hx / 2. * node.x + (xCoords[p] + xCoords[p + 1]) / 2.;
                  double y = hy / 2. * node.y + (yCoords[s] + yCoords[s + 1]) / 2.;
                  double z = hz / 2. * node.z + (zCoords[r] + zCoords[r + 1]) / 2.;

                  valueOfIntegral += node.weight * uNumMinusUReal(I, areasMesh, coefs, q, sMesh, x, y, z, 2.);
               }

               normL2 += hx * hy * hz * valueOfIntegral / 8.;
            }
         }
      }
   }

   return sqrt(normL2);
}

