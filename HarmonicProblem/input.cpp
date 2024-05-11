#include "structs.h"

// —читывание сетки
void inputMesh(Mesh &mesh)
{
   std::ifstream m("mesh.txt");

   int xMSize = 0;
   int yMSize = 0;
   int zMSize = 0;
   int numElements = 0;

   // —читываем сетку по х
   m >> xMSize;
   mesh.xM.resize(xMSize);

   for (int i = 0; i < xMSize; i++)
      m >> mesh.xM[i];

   // —читываем сетку по у
   m >> yMSize;
   mesh.yM.resize(yMSize);

   for (int i = 0; i < yMSize; i++)
      m >> mesh.yM[i];

   // —читываем сетку по z
   m >> zMSize;
   mesh.zM.resize(zMSize);

   for (int i = 0; i < zMSize; i++)
      m >> mesh.zM[i];

   // —читываем элементы сетки (области)
   m >> numElements;
   mesh.areasMesh.resize(numElements);

   for (int i = 0; i < numElements; i++)
   {
      mesh.areasMesh[i].resize(7);

      int numArea = 0;
      int x0 = 0;
      int x1 = 0;
      int y0 = 0;
      int y1 = 0;
      int z0 = 0;
      int z1 = 0;

      m >> numArea;
      m >> x0 >> x1;
      m >> y0 >> y1;
      m >> z0 >> z1;

      mesh.areasMesh[i][0] = numArea - 1;
      mesh.areasMesh[i][1] = x0 - 1;
      mesh.areasMesh[i][2] = x1 - 1;
      mesh.areasMesh[i][3] = y0 - 1;
      mesh.areasMesh[i][4] = y1 - 1;
      mesh.areasMesh[i][5] = z0 - 1;
      mesh.areasMesh[i][6] = z1 - 1;
   }

   m.close();
}

// —читываем разбиение сетки (считывать будем дл€ каждой оси просто передава€ нужные массивы)
void inputSplittingMesh(vector<double> &coordAxisOfMesh, vector<double> &coordAxisOfSplitMesh, vector<int> &shiftArray, std::ifstream &split)
{
   const int axisCoordSize = coordAxisOfMesh.size();

   shiftArray.resize(axisCoordSize);

   int nk = 0;
   coordAxisOfSplitMesh.resize(1, coordAxisOfMesh[0]);

   for (int i = 0, j = 1; i < axisCoordSize - 1; i++, j++)
   {
      int countIntervals = 0;

      double coef = 0;
      double step = 0;

      split >> countIntervals >> coef;
      nk += countIntervals;
      coordAxisOfSplitMesh.resize(nk + 1);

      if (coef != 1)
      {
         double sumProgression = (pow(coef, countIntervals) - 1.) / (coef - 1.);
         step = (coordAxisOfMesh[i + 1] - coordAxisOfMesh[i]) / sumProgression;

         int jk = 1;
         for (j; j < nk; j++, jk++)
            coordAxisOfSplitMesh[j] = coordAxisOfMesh[i] + step * (pow(coef, jk) - 1.) / (coef - 1.);
      }
      else
      {
         step = (coordAxisOfMesh[i + 1] - coordAxisOfMesh[i]) / countIntervals;

         int jk = 1;
         for (j; j < nk; j++, jk++)
            coordAxisOfSplitMesh[j] = coordAxisOfMesh[i] + step * jk;
      }

      coordAxisOfSplitMesh[j] = coordAxisOfMesh[i + 1];
      shiftArray[i + 1] = j;
   }
}

void inputBoundaryConditions(BoundaryConds &conditions)
{
   auto &bC = conditions.bC;
   int countConditions = 0;

   std::ifstream conds("conditions.txt");

   conds >> countConditions;
   bC.resize(countConditions);

   for (int i = 0; i < countConditions; i++)
   {
      bC[i].resize(8);

      int typeCond = 0;
      int numFunc = 0;

      int x0 = 0;
      int x1 = 0;
      int y0 = 0;
      int y1 = 0;
      int z0 = 0;
      int z1 = 0;

      conds >> typeCond >> numFunc;
      conds >> x0 >> x1;
      conds >> y0 >> y1;
      conds >> z0 >> z1;

      bC[i][0] = typeCond;
      bC[i][1] = numFunc - 1;
      bC[i][2] = x0 - 1;
      bC[i][3] = x1 - 1;
      bC[i][4] = y0 - 1;
      bC[i][5] = y1 - 1;
      bC[i][6] = z0 - 1;
      bC[i][7] = z1 - 1;
   }

   conds.close();
}