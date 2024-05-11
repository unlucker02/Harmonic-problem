#include "functions.h"

void localOptimalSchemeLU(SLAE &slae, SLAE &LU, LOS &vectors, int maxIter, double eps)
{
   auto &r1 = vectors.r1;
   auto &z1 = vectors.z1;
   auto &p1 = vectors.p1;
   auto &mult = vectors.mult;
   auto &rk = vectors.rk;
   auto &Ar = vectors.Ar;
   auto &p = vectors.p;

   const int sizeMatrix = slae.A.di.size();
   double normb = 0;

   r1.resize(sizeMatrix);
   z1.resize(sizeMatrix);
   p1.resize(sizeMatrix);
   mult.resize(sizeMatrix);
   rk.resize(sizeMatrix);
   Ar.resize(sizeMatrix);
   p.resize(sizeMatrix);

   calcDiscrepancy(slae, vectors, normb);

   calcY(LU, r1, r1);
   calcX(LU, r1, z1);

   multMatrixToVector(slae.A, z1, p1);
   calcY(LU, p1, p1);

   double scalarr = scalarMultiply(r1, r1);
   double discrepancy = sqrt(scalarr / normb);

   int k = 1;

   for (k; k < maxIter && discrepancy > eps; k++)
   {
      std::cout << discrepancy << "\n";
      double scalarp = scalarMultiply(p1, p1);
      double alpha = scalarMultiply(p1, r1) / scalarp;

      calcVectorMultToCoef(z1, alpha, mult);
      calcVectorPlusVector(slae.q, mult, slae.q);

      calcVectorMultToCoef(p1, -alpha, mult);
      calcVectorPlusVector(r1, mult, r1);

      calcX(LU, r1, rk);
      multMatrixToVector(slae.A, rk, Ar);
      calcY(LU, Ar, p);

      double betta = -scalarMultiply(p1, p) / scalarp;

      calcVectorMultToCoef(z1, betta, mult);
      calcVectorPlusVector(rk, mult, z1);

      calcVectorMultToCoef(p1, betta, mult);
      calcVectorPlusVector(p, mult, p1);

      discrepancy = sqrt(scalarMultiply(r1, r1) / normb);
   }
   normb = 0;
   calcDiscrepancy(slae, vectors, normb);

   discrepancy = sqrt(scalarMultiply(r1, r1) / normb);

   std::cout << discrepancy << "\n";
}

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

void multMatrixToVector(SparseMatrix &A, vector<double> &vec, vector<double> &result)
{
   const int sizeMatrix = A.di.size();

   for (int i = 0; i < sizeMatrix; i++)
   {
      result[i] = A.di[i] * vec[i];

      int i0 = A.ig[i];
      int i1 = A.ig[i + 1];

      for (i0; i0 < i1; i0++)
      {
         int j = A.jg[i0];

         result[i] += A.ggl[i0] * vec[j];
         result[j] += A.ggu[i0] * vec[i];
      }
   }
}

void calcDiscrepancy(SLAE &slae, LOS &vectors, double &normb)
{
   auto &matrix = slae.A;
   auto &b = slae.b;
   auto &q = slae.q;

   const int sizeMatrix = slae.A.di.size();

   for (int i = 0; i < sizeMatrix; i++)
   {
      normb += b[i] * b[i];
      vectors.r1[i] = b[i] - matrix.di[i] * q[i];

      int i0 = matrix.ig[i];
      int i1 = matrix.ig[i + 1];

      for (i0; i0 < i1; i0++)
      {
         int j = matrix.jg[i0];
         vectors.r1[i] -= matrix.ggl[i0] * q[j];
         vectors.r1[j] -= matrix.ggu[i0] * q[i];
      }
   }
}

void calcVectorMultToCoef(vector <double> &vec, double coef, vector<double> &result)
{
   const int size = vec.size();

   for (int i = 0; i < size; i++)
      result[i] = vec[i] * coef;
}

void calcVectorPlusVector(vector<double> &vector1, vector<double> &vector2, vector<double> &result)
{
   const int size = vector1.size();

   for (int i = 0; i < size; i++)
      result[i] = vector1[i] + vector2[i];
}

double scalarMultiply(vector<double> &vector1, vector<double> &vector2)
{
   const int size = vector1.size();
   double result = 0;

   for (int i = 0; i < size; i++)
      result += vector1[i] * vector2[i];

   return result;
}