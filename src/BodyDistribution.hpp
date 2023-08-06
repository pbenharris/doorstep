#pragma once

#include <random>
#include <cmath>


#include <spdlog/spdlog.h>

#include "System.hpp"

namespace doorstep
{
   class BodyDistribution {

   public:

      BodyDistribution(const RunConfiguration& rc, spdlog::logger& mlogger) :
         logger (mlogger)
      {
         numberBodies = rc.celestialBody.size();

         // Add count for random groups,
         // Set their properties
         for (auto it = rc.randomConfig.begin();
              it!=rc.randomConfig.end();
              it++)
         { 
            numberBodies += it->number_bodies;
         }

         // Add count of bodies for grids
         for (auto it = rc.gridConfig.begin();
              it!=rc.gridConfig.end();
              it++)
         {
            size_t gridSize = it->nx * it->ny * it->nz;
            numberBodies += gridSize;
         }

         mass = scalar_type(numberBodies, 1e-10);
         radius = scalar_type(numberBodies, 1e-10);

         // Next set up distributions for p, q
         p = container_type(numberBodies, 0.0);
         q = container_type(numberBodies, 0.0);

         // Set positions and velocities - set q to positions, q to velocities (for now)
         // Note there's a different call for each class of bodies
         size_t lastIdx = setCelestialBodies(rc, 0);
         lastIdx = setRandomClouds(rc, lastIdx);
         lastIdx = setGrids(rc, lastIdx);

         // Set the mean q and p to zero
         point_type pmean = center_of_mass(p, mass);
         point_type qmean = center_of_mass(q , mass );
         for( size_t i=0 ; i<numberBodies ; ++i )
         {
            q[i] -= qmean;
            p[i] -= pmean;
         }

         // Finally adjust p to be momentum
         for( size_t i=0 ; i<numberBodies ; ++i ) p[i] *= mass[i];
      }

      size_t bodyCount(void) const
      {
         return numberBodies;
      }

      scalar_type getMass(void) const
      {
         return mass;
      }

      scalar_type getRadius(void) const
      {
         return radius;
      }

      container_type getP(void) const
      {
	 return p;
      }

      container_type getQ(void) const
      {
	 return q;
      }

   private:

      // Set up each CelestialBody
      size_t setCelestialBodies(const RunConfiguration& rc, size_t idx)
      {  
         for (auto cfg = rc.celestialBody.begin();
              cfg!=rc.celestialBody.end(); cfg++, idx++)
         {
            mass[idx]   = cfg->mass;
            radius[idx] = cfg->radius;
   
            for (size_t j=0;j<cfg->position.size();j++)
            {
               q[idx][j] = cfg->position[j];
               p[idx][j] = cfg->velocity[j];
            }

           // p[idx] *= cfg->mass;
         }

         return idx;
      }

      // Set up each random configuration
      // Returns an index into the p and q vectors at the end of what was written
      // Note only 2D, uniform distribution is implemented
      size_t setRandomClouds(const RunConfiguration& rc, size_t idx)
      {
         const double pi = std::acos(-1.);

         std::default_random_engine rEngine;

         for (auto cfg = rc.randomConfig.begin();
              cfg!=rc.randomConfig.end();
              cfg++)
         {
                 
            std::uniform_real_distribution<double> xdist(cfg->min_x, cfg->max_x);
            std::uniform_real_distribution<double> ydist(cfg->min_y, cfg->max_y);
            std::uniform_real_distribution<double> thetaDist (-pi, pi);
            std::normal_distribution<double> speedDist (0, cfg->speed_sd);

            for (size_t n=0; n< cfg->number_bodies; n++, idx++)
            {
               double thisMass = cfg->particle_mass;
               mass[idx] = thisMass;

               // Random distribution for velocity
               double theta = thetaDist(rEngine);
               double vmag = speedDist(rEngine);

               q[idx] = cfg->initial_position;
               p[idx] = cfg->initial_velocity;

               q[idx][0] += xdist(rEngine);
               q[idx][1] += ydist(rEngine);

               p[idx][0] += vmag * std::cos(theta);
               p[idx][1] += vmag * std::sin(theta);

             //  p[idx] *= thisMass;
            }
         }

         return idx;
      }

      size_t setGrids(const RunConfiguration& rc, size_t idx)
      {
         for (auto cfg = rc.gridConfig.begin();
              cfg != rc.gridConfig.end();
              cfg++)
         {
            // Assume initial position is at the center of mass of the grid
            double d=cfg->distance;
            double minx = cfg->ic_x - 0.5*d*(cfg->nx - 1);
            double miny = cfg->ic_y - 0.5*d*(cfg->ny - 1);
            double minz = cfg->ic_z - 0.5*d*(cfg->nz - 1);

            for (int ix=0; ix!=cfg->nx; ix++)
               for (int iy=0; iy!=cfg->ny; iy++)
                  for (int iz=0; iz!=cfg->nz; iz++, idx++)
                  {
                     mass[idx] = cfg->mass;
                     
                     point_type pos;
                     pos[0] = minx + d*ix;
                     pos[1] = miny + d*iy;
                     pos[2] = minz + d*iz;
                     q[idx]=pos;
                     
                     point_type vel;
                     vel[0] = cfg->ic_vx;
                     vel[1] = cfg->ic_vy;
                     vel[2] = cfg->ic_vz;
                     p[idx] = vel;
                  } 
         }        
         
         return idx;
      }

      scalar_type mass;
      scalar_type radius;

      container_type p;
      container_type q;
      
      spdlog::logger& logger;
      size_t numberBodies;
      
   };
} // end namespace doorstep
