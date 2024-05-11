#include "functions.h"

void ProfileMatrix::calcLDU()
{
   const int sizeMatrix = di.size();

   for (int i = 0; i < sizeMatrix; i++)
   {
      int i0 = ig[i];
      int i1 = ig[i + 1];

      int m = i1 - i0;

      int k = i0;
      int j = i - m;
      double sumDi = 0;

      for (k, j; k < i1; k++, j++)
      {
         int j0 = ig[j];
         int j1 = ig[j + 1];

         int ik = i0;
         int kj = j0;

         int kol_i = k - i0;
         int kol_j = j1 - j0;
         int kol_r = kol_i - kol_j;

         if (kol_r < 0)
            kj -= kol_r;
         else
            ik += kol_r;

         int kd = 0;
         if (kol_j >= kol_i)
            kd = i - m;
         else
            kd = j - kol_j;

         double sumLow = 0;
         double sumUpper = 0;

         for (ik, kj, kd; ik < k; ik++, kj++, kd++)
         {
            sumLow += ggl[ik] * ggu[kj] * di[kd];
            sumUpper += ggu[ik] * ggl[kj] * di[kd];
         }

         ggl[k] = (ggl[k] - sumLow) / di[j];
         ggu[k] = (ggu[k] - sumUpper) / di[j];
         sumDi += ggu[k] * di[j] * ggl[k];
      }

      di[i] -= sumDi;

      if (abs(di[i]) < 1e-14)
      {
         std::cout << "LDU decomposition is impossible.";
         std::exit(-1);
      }
   }
}

void ProfileMatrix::calcZ(vector<double> &z, vector<double> &b)
{
   const int sizeMatrix = z.size();
   z = b;

   for (int i = 0; i < sizeMatrix; i++)
   {
      int i0 = ig[i];
      int i1 = ig[i + 1];

      int m = i1 - i0;

      int k = i0;
      int j = i - m;
      double sum = 0;

      for (k, j; k < i1; k++, j++)
         sum += ggl[k] * z[j];

      z[i] -= sum;
   }
}

void ProfileMatrix::calcY(vector<double> &y)
{
   const int sizeMatrix = y.size();

   for (int i = 0; i < sizeMatrix; i++)
      y[i] /= di[i];
}

void ProfileMatrix::calcX(vector<double> &x)
{
   const int sizeMatrix = x.size();

   for (int i = sizeMatrix - 1; i >= 0; i--)
   {
      int i0 = ig[i];
      int i1 = ig[i + 1];

      int k = i1 - 1;
      int j = i - 1;
      double xi = x[i];

      for (k, j; k >= i0; k--, j--)
         x[j] -= ggu[k] * xi;

      x[i] = xi;
   }
}

void ProfileMatrix::solveSlae(vector<double> &z, vector<double> &b)
{
   calcLDU();
   calcZ(z, b);
   calcY(z);
   calcX(z);
}