#include "functions.h"

int main()
{
   Mesh mesh{ };
   SplittingMesh sMesh{ };
   ShiftsArrays I{ };
   BoundaryConds conds{ };
   SLAE slae{ }, LU{ };
   FunctionsProblem funcs;
   LOS vectors{ };

   vector<double> omega = { 1e-4, 1e-1, 1e+2, 1e+5, 1e+9 };
   vector<double> lambda = { 1e+2, 1e+3, 1e+4, 1e+5, 8 * 1e+5 };
   vector<double> khi = { 8.81 * 1e-12, 1e-12, 1e-11, 1e-10 };
   vector<double> sigma = { 0., 1e+2, 1e+4, 1e+6, 1e+8 };

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

   
   //calcGlobalMatrixAndVector(slae, mesh, sMesh, I, conds, funcs);

   //funcs.omega[0] = 100.;
   //funcs.lambda[0] = 100.;
   //funcs.khi[0] = 1e-12;
   //funcs.sigma[0] = 1e+8;

   //printf_s("m    countIter     timeGMRES     normL2\n");

   //for (int m = 5; m <= 100; m += 5)
   //{
   //   auto startGMRES = std::chrono::steady_clock::now();
   //   int countIter = GMRES(slae, m, 10000, 1e-14);
   //   auto finishGMRES = std::chrono::steady_clock::now();

   //   double timeGMRES = std::chrono::duration<double, std::milli>(finishGMRES - startGMRES).count();

   //   double normFromGMRES = normL2(I, mesh.areasMesh, funcs, slae.q, sMesh);

   //   printf_s("%d %d %5lf %.15lf\n", m, countIter, timeGMRES, normFromGMRES);

   //   clearVector(slae.q);
   //}

   funcs.lambda[0] = 100.;
   funcs.khi[0] = 1e-12;
   funcs.sigma[0] = 1.;

   for (int i = 0; i < 5; i++)
   {
      funcs.omega[0] = omega[i];

      // Полностью собираем глобальную матрицу
      calcGlobalMatrixAndVector(slae, mesh, sMesh, I, conds, funcs);

      std::cout << "LOS started to solve.\n";
      auto startLOS = std::chrono::steady_clock::now();

      // Получаем решение
      calcLU(slae, LU);
      localOptimalSchemeLU(slae, LU, vectors, 10000, 1e-14);

      auto finishLOS = std::chrono::steady_clock::now();
      double timeLOS = std::chrono::duration<double, std::milli>(finishLOS - startLOS).count();

      std::cout << "LOS finished to solve.\n"
         << "Time for LOS solver: " << timeLOS << "\n\n";

      ProfileMatrix pMatrix{ };
      vector<double> qLOS = slae.q;

      clearVector(slae.q);

      // Перегенерация разреженной матрицы в профильную
      convertSparseToProfile(slae.A, pMatrix);

      vector<double> qLDU(slae.A.di.size());

      std::cout << "LDU started to solve.\n";
      auto startLDU = std::chrono::steady_clock::now();

      // Решение через LDU разложение
      pMatrix.solveSlae(qLDU, slae.b);

      auto finishLDU = std::chrono::steady_clock::now();
      double timeLDU = std::chrono::duration<double, std::milli>(finishLDU - startLDU).count();

      std::cout << "LDU finished to solve.\n"
         << "Time for LDU solver: " << timeLDU << "\n\n";

      std::cout << "GMRES started to solve.\n";
      auto startGMRES = std::chrono::steady_clock::now();

      // Получаем решение
      GMRES(slae, 10, 10000, 1e-14);

      auto finishGMRES = std::chrono::steady_clock::now();
      double timeGMRES = std::chrono::duration<double, std::milli>(finishGMRES - startGMRES).count();

      std::cout << "GMRES finished to solve.\n"
         << "Time for GMRES solver: " << timeGMRES << "\n\n";

      vector<double> qGMRES = slae.q;

      // Считается нормы ||u_числ - u_анал|| для u_числ
      // от ЛОСа, LDU и GMRES
      double normFromLOS = normL2(I, mesh.areasMesh, funcs, qLOS, sMesh);
      double normFromLDU = normL2(I, mesh.areasMesh, funcs, qLDU, sMesh);
      double normFromGMRES = normL2(I, mesh.areasMesh, funcs, qGMRES, sMesh);

      //outputForTests(sMesh, slae.q, normFromGMRES);

      outputForReseach(funcs, timeLOS, normFromLOS, timeLDU, normFromLDU, timeGMRES, normFromGMRES);

      clearSLAE(slae);
   }

   funcs.omega[0] = 100.;
   funcs.khi[0] = 1e-12;
   funcs.sigma[0] = 1.;

   for (int i = 0; i < 5; i++)
   {
      funcs.lambda[0] = lambda[i];

      // Полностью собираем глобальную матрицу
      calcGlobalMatrixAndVector(slae, mesh, sMesh, I, conds, funcs);

      std::cout << "LOS started to solve.\n";
      auto startLOS = std::chrono::steady_clock::now();

      // Получаем решение
      calcLU(slae, LU);
      localOptimalSchemeLU(slae, LU, vectors, 10000, 1e-14);

      auto finishLOS = std::chrono::steady_clock::now();
      double timeLOS = std::chrono::duration<double, std::milli>(finishLOS - startLOS).count();

      std::cout << "LOS finished to solve.\n"
         << "Time for LOS solver: " << timeLOS << "\n\n";

      ProfileMatrix pMatrix{ };
      vector<double> qLOS = slae.q;

      clearVector(slae.q);

      // Перегенерация разреженной матрицы в профильную
      convertSparseToProfile(slae.A, pMatrix);

      vector<double> qLDU(slae.A.di.size());

      std::cout << "LDU started to solve.\n";
      auto startLDU = std::chrono::steady_clock::now();

      // Решение через LDU разложение
      pMatrix.solveSlae(qLDU, slae.b);

      auto finishLDU = std::chrono::steady_clock::now();
      double timeLDU = std::chrono::duration<double, std::milli>(finishLDU - startLDU).count();

      std::cout << "LDU finished to solve.\n"
         << "Time for LDU solver: " << timeLDU << "\n\n";

      std::cout << "GMRES started to solve.\n";
      auto startGMRES = std::chrono::steady_clock::now();

      // Получаем решение
      GMRES(slae, 10, 10000, 1e-14);

      auto finishGMRES = std::chrono::steady_clock::now();
      double timeGMRES = std::chrono::duration<double, std::milli>(finishGMRES - startGMRES).count();

      std::cout << "GMRES finished to solve.\n"
         << "Time for GMRES solver: " << timeGMRES << "\n\n";

      vector<double> qGMRES = slae.q;

      // Считается нормы ||u_числ - u_анал|| для u_числ
      // от ЛОСа, LDU и GMRES
      double normFromLOS = normL2(I, mesh.areasMesh, funcs, qLOS, sMesh);
      double normFromLDU = normL2(I, mesh.areasMesh, funcs, qLDU, sMesh);
      double normFromGMRES = normL2(I, mesh.areasMesh, funcs, qGMRES, sMesh);

      //outputForTests(sMesh, slae.q, normFromGMRES);

      outputForReseach(funcs, timeLOS, normFromLOS, timeLDU, normFromLDU, timeGMRES, normFromGMRES);

      clearSLAE(slae);
   }

   funcs.lambda[0] = 100.;
   funcs.omega[0] = 100.;
   funcs.sigma[0] = 1.;

   for (int i = 0; i < 4; i++)
   {
      funcs.khi[0] = khi[i];

      // Полностью собираем глобальную матрицу
      calcGlobalMatrixAndVector(slae, mesh, sMesh, I, conds, funcs);

      std::cout << "LOS started to solve.\n";
      auto startLOS = std::chrono::steady_clock::now();

      // Получаем решение
      calcLU(slae, LU);
      localOptimalSchemeLU(slae, LU, vectors, 10000, 1e-14);

      auto finishLOS = std::chrono::steady_clock::now();
      double timeLOS = std::chrono::duration<double, std::milli>(finishLOS - startLOS).count();

      std::cout << "LOS finished to solve.\n"
         << "Time for LOS solver: " << timeLOS << "\n\n";

      ProfileMatrix pMatrix{ };
      vector<double> qLOS = slae.q;

      clearVector(slae.q);

      // Перегенерация разреженной матрицы в профильную
      convertSparseToProfile(slae.A, pMatrix);

      vector<double> qLDU(slae.A.di.size());

      std::cout << "LDU started to solve.\n";
      auto startLDU = std::chrono::steady_clock::now();

      // Решение через LDU разложение
      pMatrix.solveSlae(qLDU, slae.b);

      auto finishLDU = std::chrono::steady_clock::now();
      double timeLDU = std::chrono::duration<double, std::milli>(finishLDU - startLDU).count();

      std::cout << "LDU finished to solve.\n"
         << "Time for LDU solver: " << timeLDU << "\n\n";

      std::cout << "GMRES started to solve.\n";
      auto startGMRES = std::chrono::steady_clock::now();

      // Получаем решение
      GMRES(slae, 10, 10000, 1e-14);

      auto finishGMRES = std::chrono::steady_clock::now();
      double timeGMRES = std::chrono::duration<double, std::milli>(finishGMRES - startGMRES).count();

      std::cout << "GMRES finished to solve.\n"
         << "Time for GMRES solver: " << timeGMRES << "\n\n";

      vector<double> qGMRES = slae.q;

      // Считается нормы ||u_числ - u_анал|| для u_числ
      // от ЛОСа, LDU и GMRES
      double normFromLOS = normL2(I, mesh.areasMesh, funcs, qLOS, sMesh);
      double normFromLDU = normL2(I, mesh.areasMesh, funcs, qLDU, sMesh);
      double normFromGMRES = normL2(I, mesh.areasMesh, funcs, qGMRES, sMesh);

      //outputForTests(sMesh, slae.q, normFromGMRES);

      outputForReseach(funcs, timeLOS, normFromLOS, timeLDU, normFromLDU, timeGMRES, normFromGMRES);

      clearSLAE(slae);
   }

   funcs.lambda[0] = 100.;
   funcs.khi[0] = 1e-12;
   funcs.omega[0] = 100.;

   for (int i = 0; i < 5; i++)
   {
      funcs.sigma[0] = sigma[i];

      // Полностью собираем глобальную матрицу
      calcGlobalMatrixAndVector(slae, mesh, sMesh, I, conds, funcs);

      std::cout << "LOS started to solve.\n";
      auto startLOS = std::chrono::steady_clock::now();

      // Получаем решение
      calcLU(slae, LU);
      localOptimalSchemeLU(slae, LU, vectors, 10000, 1e-14);

      auto finishLOS = std::chrono::steady_clock::now();
      double timeLOS = std::chrono::duration<double, std::milli>(finishLOS - startLOS).count();

      std::cout << "LOS finished to solve.\n"
         << "Time for LOS solver: " << timeLOS << "\n\n";

      ProfileMatrix pMatrix{ };
      vector<double> qLOS = slae.q;

      clearVector(slae.q);

      // Перегенерация разреженной матрицы в профильную
      convertSparseToProfile(slae.A, pMatrix);

      vector<double> qLDU(slae.A.di.size());

      std::cout << "LDU started to solve.\n";
      auto startLDU = std::chrono::steady_clock::now();

      // Решение через LDU разложение
      pMatrix.solveSlae(qLDU, slae.b);

      auto finishLDU = std::chrono::steady_clock::now();
      double timeLDU = std::chrono::duration<double, std::milli>(finishLDU - startLDU).count();

      std::cout << "LDU finished to solve.\n"
         << "Time for LDU solver: " << timeLDU << "\n\n";

      std::cout << "GMRES started to solve.\n";
      auto startGMRES = std::chrono::steady_clock::now();

      // Получаем решение
      GMRES(slae, 10, 10000, 1e-14);

      auto finishGMRES = std::chrono::steady_clock::now();
      double timeGMRES = std::chrono::duration<double, std::milli>(finishGMRES - startGMRES).count();

      std::cout << "GMRES finished to solve.\n"
         << "Time for GMRES solver: " << timeGMRES << "\n\n";

      vector<double> qGMRES = slae.q;

      // Считается нормы ||u_числ - u_анал|| для u_числ
      // от ЛОСа, LDU и GMRES
      double normFromLOS = normL2(I, mesh.areasMesh, funcs, qLOS, sMesh);
      double normFromLDU = normL2(I, mesh.areasMesh, funcs, qLDU, sMesh);
      double normFromGMRES = normL2(I, mesh.areasMesh, funcs, qGMRES, sMesh);

      //outputForTests(sMesh, slae.q, normFromGMRES);

      outputForReseach(funcs, timeLOS, normFromLOS, timeLDU, normFromLDU, timeGMRES, normFromGMRES);

      clearSLAE(slae);
   }

   return 0;
}