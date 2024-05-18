#include "functions.h"

double uNumMinusUReal(ShiftsArrays &I, vector<vector<int>> &areasMesh, FunctionsProblem &funcs, vector<double> &q, SplittingMesh &sMesh, double x, double y, double z, double t)
{
   double u = uNum(I, areasMesh, funcs, q, sMesh, x, y, z, t);
   double uR = uReal(funcs, x, y, z, t);

   return (u - uR) * (u - uR);
}

double uReal(FunctionsProblem &funcs, double x, double y, double z, double t)
{
   auto &omega = funcs.omega;

   // (x * x + y * y + z * z) * sin(omega[0]() * t) + (x * x * y + z * z) * cos(omega[0]() * t)
   // sin(x + z) * sin(omega[0] * t) - cos(y + z) * cos(omega[0] * t)

   return sin(x + z) * sin(omega[0] * t) - cos(y + z) * cos(omega[0] * t);
}

double uNum(ShiftsArrays &I, vector<vector<int>> &areasMesh, FunctionsProblem &funcs, vector<double> &q, SplittingMesh &sMesh, double x, double y, double z, double t)
{
   auto &xCoord = sMesh.x;
   auto &yCoord = sMesh.y;
   auto &zCoord = sMesh.z;

   auto &omega = funcs.omega;

   const int xSize = xCoord.size();
   const int ySize = yCoord.size();
   const int zSize = zCoord.size();

   int begX = 0;
   int endX = xSize - 1;
   int begY = 0;
   int endY = ySize - 1;
   int begZ = 0;
   int endZ = zSize - 1;

   while (!(xCoord[begX] <= x && x <= xCoord[begX + 1]))
   {
      int indX = (begX + endX) / 2;

      if (xCoord[indX] < x)
         begX = indX;
      else
         endX = indX;
   }

   while (!(yCoord[begY] <= y && y <= yCoord[begY + 1]))
   {
      int indY = (begY + endY) / 2;

      if (yCoord[indY] < y)
         begY = indY;
      else
         endY = indY;
   }

   while (!(zCoord[begZ] <= z && z <= zCoord[begZ + 1]))
   {
      int indZ = (begZ + endZ) / 2;

      if (zCoord[indZ] < z)
         begZ = indZ;
      else
         endZ = indZ;
   }

   vector<int> globalNums{ };
   getGlobalNums(globalNums, xSize, ySize, begX, begY, begZ);

   int l = 0;
   int numArea = numberOfArea(I, areasMesh, begX, begY, begZ, l);

   double x0 = xCoord[begX];
   double x1 = xCoord[begX + 1];
   double y0 = yCoord[begY];
   double y1 = yCoord[begY + 1];
   double z0 = zCoord[begZ];
   double z1 = zCoord[begZ + 1];

   double hx = x1 - x0;
   double hy = y1 - y0;
   double hz = z1 - z0;

   double psi1 = (x1 - x) / hx * (y1 - y) / hy * (z1 - z) / hz;
   double psi2 = (x - x0) / hx * (y1 - y) / hy * (z1 - z) / hz;
   double psi3 = (x1 - x) / hx * (y - y0) / hy * (z1 - z) / hz;
   double psi4 = (x - x0) / hx * (y - y0) / hy * (z1 - z) / hz;
   double psi5 = (x1 - x) / hx * (y1 - y) / hy * (z - z0) / hz;
   double psi6 = (x - x0) / hx * (y1 - y) / hy * (z - z0) / hz;
   double psi7 = (x1 - x) / hx * (y - y0) / hy * (z - z0) / hz;
   double psi8 = (x - x0) / hx * (y - y0) / hy * (z - z0) / hz;

   double us = q[2 * globalNums[0]] * psi1 +
               q[2 * globalNums[1]] * psi2 +
               q[2 * globalNums[2]] * psi3 +
               q[2 * globalNums[3]] * psi4 +
               q[2 * globalNums[4]] * psi5 +
               q[2 * globalNums[5]] * psi6 +
               q[2 * globalNums[6]] * psi7 +
               q[2 * globalNums[7]] * psi8;

   double uc = q[2 * globalNums[0] + 1] * psi1 +
               q[2 * globalNums[1] + 1] * psi2 +
               q[2 * globalNums[2] + 1] * psi3 +
               q[2 * globalNums[3] + 1] * psi4 +
               q[2 * globalNums[4] + 1] * psi5 +
               q[2 * globalNums[5] + 1] * psi6 +
               q[2 * globalNums[6] + 1] * psi7 +
               q[2 * globalNums[7] + 1] * psi8;

   double u = us * sin(omega[numArea] * t) + uc * cos(omega[numArea] * t);

   return u;
}

void getGlobalNums(vector<int> &globalNums, int xSize, int ySize, int p, int s, int r)
{
   globalNums = { r * ySize * xSize + s * xSize + p,
                  r * ySize * xSize + s * xSize + p + 1,
                  r * ySize * xSize + (s + 1) * xSize + p,
                  r * ySize * xSize + (s + 1) * xSize + p + 1,
                  (r + 1) * ySize * xSize + s * xSize + p,
                  (r + 1) * ySize * xSize + s * xSize + p + 1,
                  (r + 1) * ySize * xSize + (s + 1) * xSize + p,
                  (r + 1) * ySize * xSize + (s + 1) * xSize + p + 1, };
}

void getVectorNumpsr(vector<vector<int>> &result, int p, int s, int r)
{
   result = { { p, s, r }, { p + 1, s, r }, { p, s + 1, r }, { p + 1, s + 1, r },
              { p, s, r + 1 }, { p + 1, s, r + 1 }, { p, s + 1, r + 1 }, { p + 1, s + 1, r + 1 } };
}

void getValuesFunc(vector<double> &result, SplittingMesh &sMesh, std::function<double(double, double, double)> &f, int p, int s, int r)
{
   auto &x = sMesh.x;
   auto &y = sMesh.y;
   auto &z = sMesh.z;

   result = { f(x[p], y[s], z[r]),
              f(x[p + 1], y[s], z[r]),
              f(x[p], y[s + 1], z[r]), 
              f(x[p + 1], y[s + 1], z[r]), 
              f(x[p], y[s], z[r + 1]), 
              f(x[p + 1], y[s], z[r + 1]), 
              f(x[p], y[s + 1], z[r + 1]), 
              f(x[p + 1], y[s + 1], z[r + 1]), };
}

int getSizeSet(std::set<int> &list, int i)
{
   int size = 0;

   for (auto j : list)
      if (j < i) size++;

   return size;
}

void setStringMatrixInZero(SparseMatrix &matrix, int numString)
{
   int i0 = matrix.ig[numString];
   int i1 = matrix.ig[numString + 1];
   const int countNonZeroElems = matrix.ig[matrix.di.size()];

   // Зануляем всю строку i в нижнем треугольнике
   for (i0; i0 < i1; i0++)
      matrix.ggl[i0] = 0;

   int j0 = matrix.ig[numString + 1];

   // Проходим по столбцам, начиная с (i + 1)го 
   // до последнего ненулевого элемента в верхнем треугольнике
   // Зануляем элемент, у которого номер строки совпал с numString
   for (j0; j0 < countNonZeroElems; j0++)
      if (matrix.jg[j0] == numString)
         matrix.ggu[j0] = 0;

   matrix.di[numString] = 1;
}

void convertSparseToProfile(SparseMatrix &matrix, ProfileMatrix &resultMatrix)
{
   const int sizeMatrix = matrix.di.size();

   resultMatrix.ig.resize(sizeMatrix + 1);

   for (int i = 0; i < sizeMatrix; i++)
   {
      int i0 = matrix.ig[i];
      int i1 = matrix.ig[i + 1];

      if (i1 - i0 > 0)
         resultMatrix.ig[i + 1] = resultMatrix.ig[i] + i - matrix.jg[i0];
   }

   const int countNonZeroElems = resultMatrix.ig[sizeMatrix];

   resultMatrix.ggl.resize(countNonZeroElems);
   resultMatrix.ggu.resize(countNonZeroElems);
   resultMatrix.di = matrix.di;

   for (int i = 0; i < sizeMatrix; i++)
   {
      int i0 = matrix.ig[i];
      int i1 = matrix.ig[i + 1];

      int ip0 = resultMatrix.ig[i];
      int ip1 = resultMatrix.ig[i + 1];

      int m = ip1 - ip0;
      int j = i - m;

      for (int k = ip0; k < ip1; k++, j++)
      {
         int ind = binarySearch(matrix.jg, j, i0, i1 - 1);

         if (ind != -1)
         {
            resultMatrix.ggl[k] = matrix.ggl[ind];
            resultMatrix.ggu[k] = matrix.ggu[ind];
         }
      }
   }
}

int binarySearch(vector<int> &vector, int target, int low, int high)
{
   while (low <= high)
   {
      int mid = (low + high) / 2;
      int value = vector[mid];

      if (value == target) return mid;
      else if (value < target) low = mid + 1;
      else high = mid - 1;
   }

   return -1;
}

int mu(int i)
{
   return i % 2;
}

int nu(int i)
{
   return (i / 2) % 2;
}

int upsilon(int i)
{
   return i / 4;
}

void outputForTests(SplittingMesh &sMesh, vector<double> &q, double normFromLos)
{
   auto &x = sMesh.x;
   auto &y = sMesh.y;
   auto &z = sMesh.z;

   const int xSize = x.size();
   const int ySize = y.size();
   const int zSize = z.size();

   FILE *out = NULL;
   fopen_s(&out, "test.txt", "w");

   if (out == NULL)
      exit(EXIT_FAILURE);

   fprintf_s(out, "\t\tx\t\t\t\t  y\t\t\t\t\tz\t\t\t\t  qS\t\t\t\tqC\n");

   for (int r = 0; r < zSize; r += 1)
   {
      for (int s = 0; s < ySize; s += 1)
      {
         for (int p = 0; p < xSize; p += 1)
         {
            int globalNum = r * ySize * xSize + s * xSize + p;

            fprintf_s(out, "%.15lf %.15lf %.15lf %.15lf %.15lf\n",
                             x[p], y[s], z[r], q[2 * globalNum], q[2 * globalNum + 1]);
         }
      }
   }

   fprintf_s(out, "\n||u - u*|| = %.15lf", normFromLos);
   fclose(out);
}

void outputForReseach(FunctionsProblem &funcs, double timeLOS, double normFromLOS, double timeLDU, double normFromLDU, double timeGMRES, double normFromGMRES)
{
   std::ofstream out("research.txt", std::ios::app);

   out << std::setprecision(3) << funcs.omega[0] << " "
       << std::setprecision(3) << funcs.lambda[0] << " "
       << std::setprecision(3) << funcs.khi[0] << " "
       << std::setprecision(3) << funcs.sigma[0] << " "
       << timeLOS << " " << normFromLOS << " "
       << timeLDU << " " << normFromLDU << " "
       << timeGMRES << " " << normFromGMRES << " "
       << "\n";

   out.close();
}

void clearVector(vector<double> &vec)
{
   const int sizeVector = vec.size();

   for (int i = 0; i < sizeVector; i++)
      vec[i] = 0;
}

void clearSLAE(SLAE &slae)
{
   const int sizeMatrix = slae.A.di.size();
   const int countNonZeroElems = slae.A.ig[sizeMatrix];

   for (int i = 0; i < sizeMatrix; i++)
   {
      slae.A.di[i] = 0;
      slae.q[i] = 0;
      slae.b[i] = 0;
   }

   for (int i = 0; i < countNonZeroElems; i++)
   {
      slae.A.ggl[i] = 0;
      slae.A.ggu[i] = 0;
   }
}