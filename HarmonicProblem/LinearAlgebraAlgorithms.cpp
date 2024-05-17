#include "functions.h"

void multMatrixToVector(vector<vector<double>> &matrix, vector<double> &vec, vector<double> &vectorResult)
{
   const int countRows = matrix.size();
   const int countColumns = matrix[0].size();

   for (int i = 0; i < countRows; i++)
   {
      double sum = 0;

      for (int j = 0; j < countColumns; j++)
         sum += matrix[i][j] * vec[j];

      vectorResult[i] = sum;
   }
}

void multSparseMatrixToVector(SparseMatrix &A, vector<double> &vec, vector<double> &result)
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

void multLowTriangleToVector(SparseMatrix &matrix, vector<double> &vec, vector<double> &result)
{
   const int sizeMatrix = matrix.di.size();

   for (int i = 0; i < sizeMatrix; i++)
   {
      double sum = vec[i] * matrix.di[i];

      int i0 = matrix.ig[i];
      int i1 = matrix.ig[i + 1];

      for (i0; i0 < i1; i0++)
      {
         int j = matrix.jg[i0];
         sum += matrix.ggl[i0] * vec[j];
      }

      result[i] = sum;
   }
}

void multMatrixToMatrix(vector<vector<double>> &matrix1, vector<vector<double>> &matrix2, vector<vector<double>> &matrixResult)
{
   const int countRows = matrix1.size();
   const int countColumns = matrix2[0].size();

   for (int i = 0; i < countRows; i++)
   {
      for (int j = 0; j < countColumns; j++)
      {
         double sum = 0;

         for (int k = 0; k < countRows; k++)
            sum += matrix1[i][k] * matrix2[k][j];

         matrixResult[i][j] = sum;
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

double calcNormVector(vector<double> &vec)
{
   double norm = 0;
   const int sizeVector = vec.size();

   for (int i = 0; i < sizeVector; i++)
      norm += vec[i] * vec[i];

   return sqrt(norm);
}

void vectorPlusVectorMultyToCoef(vector<double> &result, vector<double> &vector1, double coef, vector<double> &vector2)
{
   const int vectorSize = vector1.size();

   for (int i = 0; i < vectorSize; i++)
      result[i] = vector1[i] + coef * vector2[i];
}

void calcSlaeWithUpTriangle(vector<vector<double>> &H, vector<double> &x, vector<double> &b, int sizeSlae)
{
   for (int i = sizeSlae - 1; i >= 0; i--)
   {
      double sum = 0;

      for (int j = i + 1; j < sizeSlae; j++)
         sum += H[i][j] * x[j];

      x[i] = (b[i] - sum) / H[i][i];
   }
}

void equalMatrixToMatrix(vector<vector<double>> &matrix1, vector<vector<double>> &matrix2)
{
   const int countRows = matrix1.size();
   const int countColumns = matrix1[0].size();

   matrix2.resize(countRows);

   for (int i = 0; i < countRows; i++)
   {
      matrix2[i].resize(countColumns);

      for (int j = 0; j < countColumns; j++)
         matrix2[i][j] = matrix1[i][j];
   }
}

void clearMatrix(vector<vector<double>> &matrix)
{
   const int countRows = matrix.size();
   const int countColumns = matrix[0].size();

   for (int i = 0; i < countRows; i++)
      for (int j = 0; j < countColumns; j++)
         matrix[i][j] = 0;
}