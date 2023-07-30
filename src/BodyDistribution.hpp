#pragma once

#include <spdlog/spdlog.h>

#include "System.hpp"

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
            size_t gridSize = it->nx * it->ny * it->nz;
            numberBodies += gridSize;
         }

         // Next set up distributions for p, q
         mass = scalar_type(numberBodies, 0.0);
         radius = scalar_type(numberBodies, 0.0);
         metric = scalar_type(numberBodies, 0.0);
         p = container_type(numberBodies, 0.0);
         q = container_type(numberBodies, 0.0);
      }

      size_t bodyCount(void) const
      {
         return numberBodies;
      }

      scalar_type getMass(void) const
      {
         return mass;
      }

   private:

      scalar_type mass;
      scalar_type radius;
      scalar_type metric;

      container_type p;
      container_type q;
      
      spdlog::logger& logger;
      size_t numberBodies;
      
   };
} // end namespace doorstep
