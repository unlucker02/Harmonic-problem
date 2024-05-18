#include "functions.h"

int GMRES(SLAE &slae, int m, int maxIter, double eps)
{
   int countIter = 1;

   auto &A = slae.A;
   auto &q = slae.q;

   const int sizeMatrix = A.di.size();

   vector<vector<double>> V{ };
   vector<vector<double>> H{ };
   vector<double> r(sizeMatrix);
   vector<double> xTemp(sizeMatrix);
   vector<double> bTemp(sizeMatrix);
   vector<double> y(sizeMatrix);
   vector<double> Ay(sizeMatrix);
   vector<double> vAdditional(sizeMatrix);

   SLAE LU{ };

   initSizeOfMatrix(V, sizeMatrix, m); // V и H буду хранить
   initSizeOfMatrix(H, m + 1, m);      // по столбцам

   calcLU(slae, LU);

   calcVectorDiscrepancy(slae, r);
   calcY(LU, r, r);

   multLowTriangleToVector(LU.A, q, xTemp);
   calcY(LU, slae.b, bTemp);

   double normDiscrepancy = calcNormVector(r);
   double normRightPart = calcNormVector(bTemp);

   double relativeDiscrepancy = normDiscrepancy / normRightPart;

   for (int k = 1; k < maxIter && relativeDiscrepancy > eps; k++)
   {
      countIter++;
//      std::cout << k << " " << relativeDiscrepancy << "\n";

      bool elemHisZero = false;
      int mTemp = m;

      initSizeOfMatrix(V, sizeMatrix, mTemp);
      initSizeOfMatrix(H, mTemp + 1, mTemp);

      initFirstColumn(V, r);

      for (int mu = 0; mu < mTemp && !elemHisZero; mu++)
      {

         vector<double> w(sizeMatrix);
         vector<double> column(sizeMatrix);

         getColumnOfMatrix(V, column, mu);

         calcX(LU, column, y);
         multSparseMatrixToVector(A, y, Ay);
         calcY(LU, Ay, w);

         calcColumnH(H, V, w, mu);
         vAdditional = w;

         calcAdditionalVector(H, V, vAdditional, mu);

         H[mu + 1][mu] = calcNormVector(vAdditional);

         if (H[mu + 1][mu] < 1e-13)
         {
            elemHisZero = true;
            mTemp = mu + 1;
         }
         else
            if (mu < mTemp - 1)
               calcColumnV(V, H[mu + 1][mu], vAdditional, mu + 1);
      }

      vector<double> d(mTemp + 1);
      d[0] = normDiscrepancy;

      initSizeOfMatrix(V, sizeMatrix, mTemp);
      initSizeOfMatrix(H, mTemp + 1, mTemp);

      doUpperTriangle(H, d, mTemp);

      vector<double> z(mTemp);
      vector<double> Vz(sizeMatrix);

      calcSlaeWithUpTriangle(H, z, d, mTemp);

      multMatrixToVector(V, z, Vz);
      vectorPlusVectorMultyToCoef(xTemp, xTemp, 1., Vz);

      calcX(LU, xTemp, q);
      calcVectorDiscrepancy(slae, r);
      calcY(LU, r, r);

      normDiscrepancy = calcNormVector(r);
      relativeDiscrepancy = normDiscrepancy / normRightPart;

      clearMatrix(V);
      clearMatrix(H);
   }

//   std::cout << relativeDiscrepancy << "\n";

   return countIter;
}

void calcVectorDiscrepancy(SLAE &slae, vector<double> &discrepancy)
{
   auto &A = slae.A;
   auto &q = slae.q;
   auto &b = slae.b;

   const int sizeMatrix = A.di.size();

   for (int i = 0; i < sizeMatrix; i++)
   {
      discrepancy[i] = b[i] - A.di[i] * q[i];

      int i0 = A.ig[i];
      int i1 = A.ig[i + 1];

      for (i0; i0 < i1; i0++)
      {
         int j = A.jg[i0];
         discrepancy[i] -= A.ggl[i0] * q[j];
         discrepancy[j] -= A.ggu[i0] * q[i];
      }
   }
}

void initFirstColumn(vector<vector<double>> &V, vector<double> &r)
{
   const int countRows = V.size();

   double normR = calcNormVector(r);

   for (int i = 0; i < countRows; i++)
      V[i][0] = r[i] / normR;
}

void initSizeOfMatrix(vector<vector<double>> &matrix, int rows, int columns)
{
   matrix.resize(rows);

   for (int i = 0; i < rows; i++)
      matrix[i].resize(columns);
}

void initRotationMatrix(vector<vector<double>> &R, vector<vector<double>> &H, vector<vector<double>> &H0, int numColumn)
{
   const int matrixRSize = R.size();

   for (int i = 0; i < matrixRSize; i++)
      R[i][i] = 1;

   double temp = sqrt(H[numColumn][numColumn] * H[numColumn][numColumn] + H0[numColumn + 1][numColumn] * H0[numColumn + 1][numColumn]);

   double si = H0[numColumn + 1][numColumn] / temp;
   double ci = H[numColumn][numColumn] / temp;

   R[numColumn][numColumn] = ci;
   R[numColumn + 1][numColumn + 1] = ci;
   R[numColumn][numColumn + 1] = si;
   R[numColumn + 1][numColumn] = -si;
}

void doUpperTriangle(vector<vector<double>> &H, vector<double> &d, int m)
{
   vector<double> dPrev(m + 1);
   vector<vector<double>> HPrev{ };
   vector<vector<double>> H0{ };
   vector<vector<double>> R{ };

   equalMatrixToMatrix(H, H0);
   initSizeOfMatrix(R, m + 1, m + 1);

   for (int mi = 0; mi < m; mi++)
   {
      std::swap(dPrev, d);
      equalMatrixToMatrix(H, HPrev);

      initRotationMatrix(R, HPrev, H0, mi);
      multMatrixToMatrix(R, HPrev, H);
      multMatrixToVector(R, dPrev, d);

      clearMatrix(R);
   }
}

void calcColumnH(vector<vector<double>> &H, vector<vector<double>> &V, vector<double> &w, int mu)
{
   const int countRows = V.size();

   for (int lambda = 0; lambda <= mu; lambda++)
   {
      double scalarMultiplication = 0;

      for (int i = 0; i < countRows; i++)
         scalarMultiplication += V[i][lambda] * w[i];

      H[lambda][mu] = scalarMultiplication;
   }
}

void calcColumnV(vector<vector<double>> &V, double elemH, vector<double> &vAdditional, int mu)
{
   const int countRows = V.size();

   for (int i = 0; i < countRows; i++)
      V[i][mu] = vAdditional[i] / elemH;
}

void getColumnOfMatrix(vector<vector<double>> &matrix, vector<double> &column, int numColumn)
{
   const int countRows = matrix.size();

   for (int i = 0; i < countRows; i++)
      column[i] = matrix[i][numColumn];
}

void calcAdditionalVector(vector<vector<double>> &H, vector<vector<double>> &V, vector<double> &vAdditional, int mu)
{
   const int countRows = V.size();

   for (int lambda = 0; lambda <= mu; lambda++)
   {
      vector<double> temp(countRows);

      for (int i = 0; i < countRows; i++)
         temp[i] = H[lambda][mu] * V[i][lambda];

      vectorPlusVectorMultyToCoef(vAdditional, vAdditional, -1., temp);
   }
}