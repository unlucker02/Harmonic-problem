#include "matrices.h"

// Полная сборка глобальной матрицы и вектора правой части
void calcGlobalMatrixAndVector(SLAE &slae, Mesh &mesh, SplittingMesh &sMesh, ShiftsArrays &I, ParametresMesh &coefs, BoundaryConds &conds, FunctionsProblem &funcs)
{
   auto &lambda = coefs.lambda;
   auto &sigma = coefs.sigma;
   auto &khi = coefs.khi;
   auto &omega = coefs.omega;

   auto &x = sMesh.x;
   auto &y = sMesh.y;
   auto &z = sMesh.z;

   const int xSize = x.size();
   const int ySize = y.size();
   const int zSize = z.size();

   int areas = 0;

   vector<int> globalNum{ };
   vector<double> fs{ };
   vector<double> fc{ };
   vector<vector<int>> numPSR{ };

   for (int r = 0; r < zSize - 1; r++)
   {
      for (int s = 0; s < ySize - 1; s++)
      {
         for (int p = 0; p < xSize - 1; p++)
         {
            int numArea = numberOfArea(I, mesh.areasMesh, p, s, r, areas);
            getGlobalNums(globalNum, xSize, ySize, p, s, r);

            if (numArea != -1)
            {
               double hx = x[p + 1] - x[p];
               double hy = y[s + 1] - y[s];
               double hz = z[r + 1] - z[r];

               for (int i = 0; i < 8; i++)
               {
                  getValuesFunc(fs, sMesh, funcs.fs[numArea], p, s, r);
                  getValuesFunc(fc, sMesh, funcs.fc[numArea], p, s, r);

                  double sumbiS = 0;
                  double sumbiC = 0;

                  int mui = mu(i);
                  int nui = nu(i);
                  int upsiloni = upsilon(i);

                  for (int j = 0; j < 8; j++)
                  {
                     int muj = mu(j);
                     int nuj = nu(j);
                     int upsilonj = upsilon(j);

                     double stiffnessij = lambda[numArea]() *
                        (1. / hx * G[mui][muj] * hy * M[nui][nuj] * hz * M[upsiloni][upsilonj] +
                         hx * M[mui][muj] * 1. / hy * G[nui][nuj] * hz * M[upsiloni][upsilonj] +
                         hx * M[mui][muj] * hy * M[nui][nuj] * 1. / hz * G[upsiloni][upsilonj]);

                     double massij = hx * hy * hz * M[mui][muj] * M[nui][nuj] * M[upsiloni][upsilonj];

                     double pij = stiffnessij - omega[numArea]() * omega[numArea]() * khi[numArea]() * massij;

                     addElemLocalMatrixInGlobalMatrix(slae.A, pij, 2 * globalNum[i], 2 * globalNum[j]);
                     addElemLocalMatrixInGlobalMatrix(slae.A, pij, 2 * globalNum[i] + 1, 2 * globalNum[j] + 1);

                     double cij = omega[numArea]() * sigma[numArea]() * massij;

                     addElemLocalMatrixInGlobalMatrix(slae.A, -cij, 2 * globalNum[i], 2 * globalNum[j] + 1);
                     addElemLocalMatrixInGlobalMatrix(slae.A, cij, 2 * globalNum[i] + 1, 2 * globalNum[j]);

                     sumbiS += fs[j] * massij;
                     sumbiC += fc[j] * massij;
                  }

                  slae.b[2 * globalNum[i]] += sumbiS;
                  slae.b[2 * globalNum[i] + 1] += sumbiC;
               }
            }
            else
            {
               getVectorNumpsr(numPSR, p, s, r);
               int l1 = 0;

               for (int i = 0; i < 8; i++)
               {
                  if (isFictitiousNode(I, mesh.areasMesh, numPSR[i][0], numPSR[i][1], numPSR[i][2], l1))
                  {
                     slae.A.di[2 * globalNum[i]] = 1;
                     slae.A.di[2 * globalNum[i] + 1] = 1;
                  }
               }
            }
         }
      }
   }

   std::cout << "The adding conditions started.\n";

   auto start = std::chrono::steady_clock::now();
   addBoundConditionTypeFirst(slae, I, sMesh, conds, funcs);
   auto finish = std::chrono::steady_clock::now();

   double timeAddCond = std::chrono::duration<double, std::milli>(finish - start).count();

   std::cout << "The adding conditions finished.\n"
             << "Time for adding conditions: " << timeAddCond << "\n\n";
}

void addBoundConditionTypeFirst(SLAE &slae, ShiftsArrays &I, SplittingMesh &sMesh, BoundaryConds &conds, FunctionsProblem &funcs)
{
   auto &matrix = slae.A;
   auto &ugS = funcs.ugS;
   auto &ugC = funcs.ugC;
   auto &conditions = conds.bC;

   auto &x = sMesh.x;
   auto &y = sMesh.y;
   auto &z = sMesh.z;

   const int countConditions = conditions.size();
   const int xSize = x.size();
   const int ySize = y.size();
   const int zSize = z.size();
   const int matrixSize = 2 * xSize * ySize * zSize;

   for (int cond = 0; cond < countConditions; cond++)
   {
      const int typeCond = conditions[cond][0];

      if (typeCond == 1)
      {
         int numFunc = conditions[cond][1];
         int p0 = conditions[cond][2];
         int p1 = conditions[cond][3];
         int s0 = conditions[cond][4];
         int s1 = conditions[cond][5];
         int r0 = conditions[cond][6];
         int r1 = conditions[cond][7];

         if (p0 == p1)
         {
            int p = I.Ix[p0];
            
            s0 = I.Iy[s0];
            s1 = I.Iy[s1];
            r0 = I.Iz[r0];
            r1 = I.Iz[r1];

            for (int r = r0; r <= r1; r++)
            {
               for (int s = s0; s <= s1; s++)
               {
                  int globalNum = r * ySize * xSize + s * ySize + p;
                  int globalNum1 = 2 * globalNum;
                  int globalNum2 = 2 * globalNum + 1;

                  setStringMatrixInZero(matrix, globalNum1);
                  setStringMatrixInZero(matrix, globalNum2);

                  slae.b[globalNum1] = ugS[numFunc](x[p], y[s], z[r]);
                  slae.b[globalNum2] = ugC[numFunc](x[p], y[s], z[r]);
               }
            }
         }
         else if (s0 == s1)
         {
            int s = I.Iy[s0];

            p0 = I.Ix[p0];
            p1 = I.Ix[p1];
            r0 = I.Iz[r0];
            r1 = I.Iz[r1];

            for (int r = r0; r <= r1; r++)
            {
               for (int p = p0; p <= p1; p++)
               {
                  int globalNum = r * ySize * xSize + s * ySize + p;
                  int globalNum1 = 2 * globalNum;
                  int globalNum2 = 2 * globalNum + 1;

                  setStringMatrixInZero(matrix, globalNum1);
                  setStringMatrixInZero(matrix, globalNum2);

                  slae.b[globalNum1] = ugS[numFunc](x[p], y[s], z[r]);
                  slae.b[globalNum2] = ugC[numFunc](x[p], y[s], z[r]);
               }
            }
         }
         else
         {
            int r = I.Iz[r0];

            p0 = I.Ix[p0];
            p1 = I.Ix[p1];
            s0 = I.Iy[s0];
            s1 = I.Iy[s1];

            for (int s = s0; s <= s1; s++)
            {
               for (int p = p0; p <= p1; p++)
               {
                  int globalNum = r * ySize * xSize + s * ySize + p;
                  int globalNum1 = 2 * globalNum;
                  int globalNum2 = 2 * globalNum + 1;

                  setStringMatrixInZero(matrix, globalNum1);
                  setStringMatrixInZero(matrix, globalNum2);

                  slae.b[globalNum1] = ugS[numFunc](x[p], y[s], z[r]);
                  slae.b[globalNum2] = ugC[numFunc](x[p], y[s], z[r]);
               }
            }
         }
      }
   }
}

// Генерация портрета матрицы в разреженном формате
void generatePortraitSparseMatrix(SLAE &slae, int xSize, int ySize, int zSize)
{
   auto &matrix = slae.A;

   const int sizeMatrix = 2 * xSize * ySize * zSize;

   matrix.di.resize(sizeMatrix);
   slae.q.resize(sizeMatrix);
   slae.b.resize(sizeMatrix);

   vector<std::set<int>> list(sizeMatrix);

   for (int r = 0; r < zSize - 1; r++)
      for (int s = 0; s < ySize - 1; s++)
         for (int p = 0; p < xSize - 1; p++)
         {
            vector<int> globalNum{ };
            getGlobalNums(globalNum, xSize, ySize, p, s, r);
            
            for (int i = 0; i < 8; i++)
            {
               int ind11 = 2 * globalNum[i];
               int ind12 = 2 * globalNum[i] + 1;

               for (int j = 0; j < 8; j++)
               {
                  int ind21 = 2 * globalNum[j];
                  int ind22 = 2 * globalNum[j] + 1;

                  list[ind11].insert(ind21);
                  list[ind11].insert(ind22);
                  list[ind12].insert(ind21);
                  list[ind12].insert(ind22);
               }
            }
         }

   matrix.ig.resize(sizeMatrix + 1);
   matrix.ig[0] = 0;
   matrix.ig[1] = 0;

   for (int i = 0; i < sizeMatrix; i++)
      matrix.ig[i + 1] = matrix.ig[i] + getSizeSet(list[i], i);

   int numNonZeroElems = matrix.ig[sizeMatrix];
   matrix.ggl.resize(numNonZeroElems);
   matrix.ggu.resize(numNonZeroElems);
   matrix.jg.resize(numNonZeroElems);

   for (int i = 0, k = 0; i < sizeMatrix; i++)
      for (auto j : list[i])
      {
         if (j >= i) break;
         matrix.jg[k] = j;
         k++;
      }
}

// Добавление элемента в матрицу в разреженном формате
void addElemLocalMatrixInGlobalMatrix(SparseMatrix &matrix, double elem, int i, int j)
{
   auto &ig = matrix.ig, &jg = matrix.jg;

   if (i == j)
      matrix.di[i] += elem;
   else
   {
      // В нижний треугольник через бинарный поиск
      if (i > j)
      {
         int beg = ig[i];
         int end = ig[i + 1] - 1;

         while (jg[beg] != j)
         {
            int ind = (beg + end) / 2;

            if (jg[ind] < j)
               beg = ind + 1;
            else
               end = ind;
         }

         matrix.ggl[beg] += elem;
      }
      else // В верхний треугольник через бинарный поиск
      {
         int beg = ig[j];
         int end = ig[j + 1] - 1;

         while (jg[beg] != i)
         {
            int ind = (beg + end) / 2;

            if (jg[ind] < i)
               beg = ind + 1;
            else
               end = ind;
         }

         matrix.ggu[beg] += elem;
      }
   }
}