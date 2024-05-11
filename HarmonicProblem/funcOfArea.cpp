#include "functions.h"

int numberOfArea(ShiftsArrays &I, vector<vector<int>> &areasMesh, int p, int s, int r, int &l)
{
   auto &Ix = I.Ix, &Iy = I.Iy, &Iz = I.Iz;

   const int countAreas = areasMesh.size();
   const int L1 = l;

   for (l; l < countAreas; l++)
   {
      const int numArea = areasMesh[l][0];
      const int p0 = Ix[areasMesh[l][1]];
      const int p1 = Ix[areasMesh[l][2]];
      const int s0 = Iy[areasMesh[l][3]];
      const int s1 = Iy[areasMesh[l][4]];
      const int r0 = Iz[areasMesh[l][5]];
      const int r1 = Iz[areasMesh[l][6]];

      if (p0 <= p && p <= p1 && p0 <= (p + 1) && (p + 1) <= p1 &&
          s0 <= s && s <= s1 && s0 <= (s + 1) && (s + 1) <= s1 &&
          r0 <= r && r <= r1 && r0 <= (r + 1) && (r + 1) <= r1)
         return numArea;
   }

   for (l = 0; l < L1; l++)
   {
      const int numArea = areasMesh[l][0];
      const int p0 = Ix[areasMesh[l][1]];
      const int p1 = Ix[areasMesh[l][2]];
      const int s0 = Iy[areasMesh[l][3]];
      const int s1 = Iy[areasMesh[l][4]];
      const int r0 = Iz[areasMesh[l][5]];
      const int r1 = Iz[areasMesh[l][6]];

      if (p0 <= p && p <= p1 && p0 <= (p + 1) && (p + 1) <= p1 &&
          s0 <= s && s <= s1 && s0 <= (s + 1) && (s + 1) <= s1 &&
          r0 <= r && r <= r1 && r0 <= (r + 1) && (r + 1) <= r1)
         return numArea;
   }

   return -1;
}

bool isFictitiousNode(ShiftsArrays &I, vector<vector<int>> &areasMesh, int p, int s, int r, int &l)
{
   auto &Ix = I.Ix, &Iy = I.Iy, &Iz = I.Iz;

   const int countAreas = areasMesh.size();
   const int L1 = l;

   for (l; l < countAreas; l++)
   {
      const int p0 = Ix[areasMesh[l][1]];
      const int p1 = Ix[areasMesh[l][2]];
      const int s0 = Iy[areasMesh[l][3]];
      const int s1 = Iy[areasMesh[l][4]];
      const int r0 = Iz[areasMesh[l][5]];
      const int r1 = Iz[areasMesh[l][6]];

      if (p0 <= p && p <= p1 &&
          s0 <= s && s <= s1 &&
          r0 <= r && r <= r1)
         return false;
   }

   for (l = 0; l < L1; l++)
   {
      const int p0 = Ix[areasMesh[l][1]];
      const int p1 = Ix[areasMesh[l][2]];
      const int s0 = Iy[areasMesh[l][3]];
      const int s1 = Iy[areasMesh[l][4]];
      const int r0 = Iz[areasMesh[l][5]];
      const int r1 = Iz[areasMesh[l][6]];

      if (p0 <= p && p <= p1 &&
          s0 <= s && s <= s1 &&
          r0 <= r && r <= r1)
         return false;
   }

   return true;
}