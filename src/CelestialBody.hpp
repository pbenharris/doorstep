#pragma once

#include <vector>

namespace doorstep
{
   struct CelestialBody
   {
      std::vector<double> position;
      std::vector<double> velocity;
      double radius;
      double mass;
      std::string name;
   };
}
