#include "functions.h"

void calcLU(SLAE &slae, SLAE &LU)
{
   auto &ig = slae.A.ig;
   auto &jg = slae.A.jg;

   auto &L = LU.A.ggl;
   auto &U = LU.A.ggu;

   LU.A.ig = slae.A.ig;
   LU.A.jg = slae.A.jg;
   LU.b = slae.b;

   const int matrixSize = slae.A.di.size();
   const int countNonZeroElems = slae.A.ig[matrixSize];

   LU.A.di.resize(matrixSize);
   L.resize(countNonZeroElems);
   U.resize(countNonZeroElems);

   for (int i = 0; i < matrixSize; i++)
   {
      double sumDi = 0;

      int i0 = ig[i];
      int i1 = ig[i + 1];

      for (int k = i0; k < i1; k++)
      {
         double sumLow = 0;
         double sumUpper = 0;

         int j = jg[k];
         int j0 = ig[j];
         int j1 = ig[j + 1];

         for (int ik = i0, kj = i0; ik < i1 && kj < j1;)
         {
            if (jg[ik] > jg[kj]) kj++;
            else if (jg[ik] < jg[kj]) ik++;
            else
            {
               sumLow += L[ik] * U[kj];
               sumUpper += L[kj] * U[ik];
               ik++;
               kj++;
            }
         }

         L[k] = (slae.A.ggl[k] - sumLow);
         U[k] = (slae.A.ggu[k] - sumUpper) / LU.A.di[j];

         sumDi += L[k] * U[k];
      }

      LU.A.di[i] = slae.A.di[i] - sumDi;
   }
}

void calcY(SLAE &LU, vector<double> &b, vector<double> &y)
{
   const int sizeMatrix = LU.A.di.size();

   for (int i = 0; i < sizeMatrix; i++)
   {
      double sum = 0;

      int i0 = LU.A.ig[i];
      int i1 = LU.A.ig[i + 1];

      for (i0; i0 < i1; i0++)
      {
         int j = LU.A.jg[i0];
         sum += LU.A.ggl[i0] * y[j];
      }

      y[i] = (b[i] - sum) / LU.A.di[i];
   }
}

void calcX(SLAE &LU, vector<double> &y, vector<double> &x)
{
   const int sizeMatrix = LU.A.di.size();

   vector <double> vectors = y;

   for (int i = sizeMatrix - 1; i >= 0; i--)
   {
      x[i] = vectors[i];

      int i0 = LU.A.ig[i];
      int i1 = LU.A.ig[i + 1];

      for (i0; i0 < i1; i0++)
      {
         int j = LU.A.jg[i0];
         vectors[j] -= x[i] * LU.A.ggu[i0];
      }
   }
}