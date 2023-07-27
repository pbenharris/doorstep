#pragma once

#include <spdlog/spdlog.h>

namespace doorstep
{
   class BodyDistribution {

   public:

      BodyDistribution(const RunConfiguration rc, spdlog::logger& mlogger) :
         logger (mlogger)
      {
         numberBodies = rc.celestialBody.size();

         // Add count for random groups
         for (auto it = rc.randomConfig.begin();
              it!=rc.randomConfig.end();
              it++)
            numberBodies += it->number_bodies;

         // Add count of bodies for grids
         for (auto it = rc.gridConfig.begin();
              it!=rc.gridConfig.end();
              it++)
         {
            gridSize = it->nx * it->ny * it->nz;
            numberBodies += gridSize;
         }

         // Next use private functions to set up distributions in p, q
         
      }

      size_t getNumberBodies(void) const
      {
         return numberBodies;
      }
      
   private:
      
      spdlog::logger& logger;
      size_t numberBodies;
      
   };
} // end namespace doorstep
