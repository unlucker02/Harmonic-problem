#include "structs.h"

// Правая часть и функции краевых
struct FunctionsProblem
{
   // Коэффициенты задачи (всегда константы)
   vector<double> lambda
   {
      100
   };

   vector<double> sigma
   {
      1
   };

   vector<double> omega
   {
      1
   };

   vector<double> khi
   {
      1
   };

   vector<std::function<double(double, double, double)>> fs
   {
      [this](double x, double y, double z) { return 2. * lambda[0] * sin(x + z) - omega[0] * sigma[0] * (-cos(y + z)) - omega[0] * omega[0] * khi[0] * sin(x + z); }, // 1
      [](double x, double y, double z) { return x * x * x + 216; }, // 2
      [](double x, double y, double z) { return 216 + y * y * y; } // 3
   }; // -3 * x - 2 * y - 2 * z
   // -12 - 20 * x * x * y - 48 * x * x - 48 * y * y - 68 * z * z
   // 2. * lambda[0] * sin(x + z) - omega[0] * sigma[0] * (-cos(y + z)) - omega[0] * omega[0] * khi[0] * sin(x + z)

   vector<std::function<double(double, double, double)>> fc
   {
      [this](double x, double y, double z) { return -2. * lambda[0] * cos(y + z) + omega[0] * sigma[0] * sin(x + z) - omega[0] * omega[0] * khi[0] * (-cos(y + z)); }, // 1
      [](double x, double y, double z) { return x * x * x + 216; }, // 2
      [](double x, double y, double z) { return 216 + y * y * y; } // 3
   }; // -4 - 4 * y + 20 * x * x - 48 * x * x * y + 20 * y * y - 28 * z * z
   // -2. * lambda[0] * cos(y + z) + omega[0] * sigma[0] * sin(x + z) - omega[0] * omega[0] * khi[0] *  (- cos(y + z))

   vector<std::function<double(double, double, double)>> ugS
   {
      [](double x, double y, double z) { return sin(x); }, // 1
      [](double x, double y, double z) { return sin(4 + z); }, // 2
      [](double x, double y, double z) { return sin(x + 4); }, // 3
      [](double x, double y, double z) { return sin(z); }, // 4
      [](double x, double y, double z) { return sin(x + z); }, // 5
      [](double x, double y, double z) { return sin(x + z); }, // 6
      [](double x, double y, double z) { return 27 + y * y * y; }, // 7
      [](double x, double y, double z) { return x * x * x + 1; }, // 8
      [](double x, double y, double z) { return 3 + y; }, // 9
      [](double x, double y, double z) { return x + 2; }, // 10
      [](double x, double y, double z) { return 2 + y; }, // 11
      [](double x, double y, double z) { return x + 1; } // 12
   };

   vector<std::function<double(double, double, double)>> ugC
   {
      [](double x, double y, double z) { return -cos(y); }, // 1
      [](double x, double y, double z) { return -cos(y + z); }, // 2
      [](double x, double y, double z) { return -cos(y + 4); }, // 3
      [](double x, double y, double z) { return -cos(y + z); }, // 4
      [](double x, double y, double z) { return -cos(z); }, // 5
      [](double x, double y, double z) { return -cos(4 + z); }, // 6
      [](double x, double y, double z) { return 27 + y * y * y; }, // 7
      [](double x, double y, double z) { return x * x * x + 1; }, // 8
      [](double x, double y, double z) { return 3 + y; }, // 9
      [](double x, double y, double z) { return x + 2; }, // 10
      [](double x, double y, double z) { return 2 + y; }, // 11
      [](double x, double y, double z) { return x + 1; } // 12
   };
};

// Считываение сетки, разбиения и КУ
void inputMesh(Mesh &mesh);
void inputSplittingMesh(vector<double> &coordAxisOfMesh, vector<double> &coordAxisOfSplitMesh, vector<int> &shiftArray, std::ifstream &split);
void inputBoundaryConditions(BoundaryConds &conditions);

// Генерация портретра для матрицы размера 2n
void generatePortraitSparseMatrix(SLAE &slae, int xSize, int ySize, int zSize);

// Добавление элемента в глобальную матрицу
void addElemLocalMatrixInGlobalMatrix(SparseMatrix &matrix, double elem, int i, int j);

// Полностью собираем глоб матрицу и вектор правой части (краевые в ней же учитываются)
void calcGlobalMatrixAndVector(SLAE &slae, Mesh &mesh, SplittingMesh &sMesh, ShiftsArrays &I, BoundaryConds &conds, FunctionsProblem &funcs);
void addBoundConditionTypeFirst(SLAE &slae, ShiftsArrays &I, SplittingMesh &sMesh, BoundaryConds &conds, FunctionsProblem &funcs);

// Функции для нахождения области, в которой находится элемент
// и для определения фиктивных узлов сетки
int numberOfArea(ShiftsArrays &I, vector<vector<int>> &areasMesh, int p, int s, int r, int &l);
bool isFictitiousNode(ShiftsArrays &I, vector<vector<int>> &areasMesh, int p, int s, int r, int &l);



// ЛОС(LU) и функции для него
void localOptimalSchemeLU(SLAE &slae, SLAE &LU, LOS &vectors, int maxIter, double eps);
void calcLU(SLAE &slae, SLAE &LU);
void calcY(SLAE &LU, vector<double> &b, vector<double> &y);
void calcX(SLAE &LU, vector<double> &y, vector<double> &x);
void calcDiscrepancy(SLAE &slae, LOS &vectors, double &normb);


// GMRES_LU и функции для него
int GMRES(SLAE &slae, int m, int maxIter, double eps);
void calcVectorDiscrepancy(SLAE &slae, vector<double> &discrepancy);
void initFirstColumn(vector<vector<double>> &V, vector<double> &r);
void initSizeOfMatrix(vector<vector<double>> &matrix, int rows, int columns);
void initRotationMatrix(vector<vector<double>> &R, vector<vector<double>> &H, vector<vector<double>> &H0, int numColumn);
void doUpperTriangle(vector<vector<double>> &H, vector<double> &d, int m);
void calcColumnH(vector<vector<double>> &H, vector<vector<double>> &V, vector<double> &w, int mu);
void calcColumnV(vector<vector<double>> &V, double elemH, vector<double> &vAdditional, int mu);
void getColumnOfMatrix(vector<vector<double>> &matrix, vector<double> &column, int numColumn);
void calcAdditionalVector(vector<vector<double>> &H, vector<vector<double>> &V, vector<double> &vAdditional, int mu);





// Функции линейной алгебры
void multMatrixToVector(vector<vector<double>> &matrix, vector<double> &vec, vector<double> &vectorResult);
void multSparseMatrixToVector(SparseMatrix &A, vector<double> &vec, vector<double> &result);
void multLowTriangleToVector(SparseMatrix &matrix, vector<double> &vec, vector<double> &result);
void multMatrixToMatrix(vector<vector<double>> &matrix1, vector<vector<double>> &matrix2, vector<vector<double>> &matrixResult);
void calcVectorMultToCoef(vector <double> &vec, double coef, vector<double> &result);
void calcVectorPlusVector(vector<double> &vector1, vector<double> &vector2, vector<double> &result);
void vectorPlusVectorMultyToCoef(vector<double> &result, vector<double> &vector1, double coef, vector<double> &vector2);
void calcSlaeWithUpTriangle(vector<vector<double>> &H, vector<double> &x, vector<double> &b, int sizeSlae);
void equalMatrixToMatrix(vector<vector<double>> &matrix1, vector<vector<double>> &matrix2);
double scalarMultiply(vector<double> &vector1, vector<double> &vector2);
double calcNormVector(vector<double> &vec);
void clearMatrix(vector<vector<double>> &matrix);


// Функции для определения номеров в одномерных лок матрицах
int mu(int i);
int nu(int i);
int upsilon(int i);


// Численное интегрирование для нормы разности в L2
void getQuadratures3D(vector<QuadratureNode> &quadratures);
double normL2(ShiftsArrays &I, vector<vector<int>> &areasMesh, FunctionsProblem &funcs, vector<double> &q, SplittingMesh &sMesh);


// Вспомогательные функции
double uNumMinusUReal(ShiftsArrays &I, vector<vector<int>> &areasMesh, FunctionsProblem &funcs, vector<double> &q, SplittingMesh &sMesh, double x, double y, double z, double t);
double uReal(FunctionsProblem &funcs, double x, double y, double z, double t);
double uNum(ShiftsArrays &I, vector<vector<int>> &areasMesh, FunctionsProblem &funcs, vector<double> &q, SplittingMesh &sMesh, double x, double y, double z, double t);
void getGlobalNums(vector<int> &globalNums, int xSize, int ySize, int p, int s, int r);
void getVectorNumpsr(vector<vector<int>> &result, int p, int s, int r);
void getValuesFunc(vector<double> &result, SplittingMesh &sMesh, std::function<double(double, double, double)> &f, int p, int s, int r);
int getSizeSet(std::set<int> &list, int i);
void setStringMatrixInZero(SparseMatrix &matrix, int numString);
int binarySearch(vector<int> &vector, int target, int low, int high);
void convertSparseToProfile(SparseMatrix &matrix, ProfileMatrix &resultMatrix);
void clearVector(vector<double> &vec);
void clearSLAE(SLAE &slae);
void outputForTests(SplittingMesh &sMesh, vector<double> &q, double normFromLos);
void outputForReseach(FunctionsProblem &funcs, double timeLOS, double normFromLOS, double timeLDU, double normFromLDU, double timeGMRES, double normFromGMRES);
