#include "functions.h"

int main()
{
   Mesh mesh{ };
   SplittingMesh sMesh{ };
   ShiftsArrays I{ };
   BoundaryConds conds{ };
   SLAE slae{ }, LU{ };
   ParametresMesh coefs;
   FunctionsProblem funcs;
   LOS vectors{ };
   ProfileMatrix pMatrix{ };

   // Считываем изначальную сетку
   inputMesh(mesh);

   // Считываем разбиение сетки
   std::ifstream split("splittingMesh.txt");
   inputSplittingMesh(mesh.xM, sMesh.x, I.Ix, split);
   inputSplittingMesh(mesh.yM, sMesh.y, I.Iy, split);
   inputSplittingMesh(mesh.zM, sMesh.z, I.Iz, split);
   split.close();

   // Считываем краевые условия
   inputBoundaryConditions(conds);

   // Строим портрет матрицы размера 2n
   generatePortraitSparseMatrix(slae, sMesh.x.size(), sMesh.y.size(), sMesh.z.size());

   // Полностью собираем глобальную матрицу
   calcGlobalMatrixAndVector(slae, mesh, sMesh, I, coefs, conds, funcs);
   
   std::cout << "LOS started to solve.\n";
   auto startLOS = std::chrono::steady_clock::now();
   
   // Получаем решение
   calcLU(slae, LU);
   localOptimalSchemeLU(slae, LU, vectors, 10000, 1e-14);
   
   auto finishLOS = std::chrono::steady_clock::now();
   double timeLOS = std::chrono::duration<double, std::milli>(finishLOS - startLOS).count();

   std::cout << "LOS finished to solve.\n"
      << "Time for LOS solver: " << timeLOS << "\n\n";

   // Перегенерация разреженной матрицы в профильную
   convertSparseToProfile(slae.A, pMatrix);

   vector<double> q(slae.A.di.size());

   std::cout << "LDU started to solve.\n";
   auto startLDU = std::chrono::steady_clock::now();

   // Решение через LDU разложение
   pMatrix.solveSlae(q, slae.b);

   auto finishLDU = std::chrono::steady_clock::now();
   double timeLDU = std::chrono::duration<double, std::milli>(finishLDU - startLDU).count();

   std::cout << "LDU finished to solve.\n"
             << "Time for LDU solver: " << timeLDU << "\n\n";

   // Считается нормы ||u_числ - u_анал|| для u_числ
   // от ЛОСа и от LDU 
   double normFromLOS = normL2(I, mesh.areasMesh, coefs, slae.q, sMesh);
   double normFromLDU = normL2(I, mesh.areasMesh, coefs, q, sMesh);

//   outputForTests(sMesh, q, normFromLDU);

   outputForReseach(coefs, timeLOS, normFromLOS, timeLDU, normFromLDU);

   return 0;
}