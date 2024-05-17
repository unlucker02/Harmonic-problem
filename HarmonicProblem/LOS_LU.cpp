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

   multSparseMatrixToVector(slae.A, z1, p1);
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
      multSparseMatrixToVector(slae.A, rk, Ar);
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