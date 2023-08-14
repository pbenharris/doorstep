#pragma once

#include <cmath>
#include <limits>
#include <set>
#include <string>
#include <exception>

#include "System.hpp"
#include "Configuration.hpp"

namespace doorstep
{
   class WindTunnel
   {
      public:

         enum Exit
         {
            DISTANCE_FROM_ORIGIN,
            BOX,
            ENERGY_LEVEL,
            UNSPECIFIED
         };

         WindTunnel(container_type& inputP, container_type& inputQ,
                    scalar_type& inputMass, mask_type& inputActive,
                    spdlog::logger& mlogger) :
            p(inputP), 
            q(inputQ),
            mass(inputMass),
            active(inputActive),
            exitType(WindTunnel::Exit::UNSPECIFIED),
            exitDistance(0.),
            logger(mlogger)
         {
            for (size_t i=0; i<inputActive.size();i++)
               if (!active[i]) inactiveParticles.insert(i);
            size_t iCount = inactiveParticles.size();
            minInactive = iCount;
            maxInactive = iCount;
            logger.info("Total {} bodies inactive.",iCount);
         };

         size_t getMinInactive(void)
         {
            return minInactive;
         }

         size_t getMaxInactive(void)
         {
            return maxInactive;
         }
      
         void setExitAsDistance(double distance)
         {
            exitType = WindTunnel::Exit::DISTANCE_FROM_ORIGIN;
            exitDistance = distance;
         }

      void setEntry(std::vector<GridConfiguration>& gridConfig)
      {
         const double pi = std::acos(-1.);
         
         size_t g_count=0;
         for (auto gc = gridConfig.begin();
              gc != gridConfig.end();
              gc ++, g_count++)
         {

            double v_mag = std::sqrt(gc->ic_vx*gc->ic_vx + gc->ic_vy*gc->ic_vy + gc->ic_vz*gc->ic_vz);
            double delta_t_e = gc->distance / v_mag;
            logger.info("Grid {} windtunnel entry time rate t_e is {}", g_count, delta_t_e);
            
            double theta = std::atan2(gc->ic_vy, gc->ic_vx);
            logger.info("Entry angle is {}", theta*180/pi);

            double d = gc->distance;
            double minx = gc->ic_x - 0.5*d*(gc->nx - 1);
            double miny = gc->ic_y - 0.5*d*(gc->ny - 1);
            double minz = gc->ic_z - 0.5*d*(gc->nz - 1);
            double maxx = gc->ic_x + 0.5*d*(gc->nx - 1);
            double maxy = gc->ic_y + 0.5*d*(gc->ny - 1);
            double maxz = gc->ic_z + 0.5*d*(gc->nz - 1);

            point_type thisp;
            thisp[0] = gc->mass * gc->ic_vx;
            thisp[1]  = gc->mass * gc->ic_vy;
            thisp[2]  = gc->mass * gc->ic_vz;
            
            // East bound
            if (std::fabs(theta)<2*std::numeric_limits<double>::epsilon())
            {
               logger.info("Configuring grid {} with east bound particles.", g_count);
               for (int i=0; i<gc->ny; i++)
               {
                  point_type x;
                  x[0] = minx;
                  x[1] = miny + d*i;
                  x[2] = 0;

                  gridEntryQ.push_back(x);
                  gridEntryP.push_back(thisp);
                  gridEntryMass.push_back(gc->mass);
                  
                  nextEntryTime.push_back(delta_t_e);
                  entryRate.push_back(delta_t_e);
               }        
            }
         }
      }

      void exitParticles(void)
      {
         if (exitType == Exit::DISTANCE_FROM_ORIGIN)
         {
            for( size_t i=0; i<q.size() ; ++i )
            {
               if (active[i] && (abs(q[i])>exitDistance))
               {
                  active[i]=false;
                  inactiveParticles.insert(i);

                  updateInactiveStats();
               }
            }
         }
      }

      void enterParticles(double t)
      {

         bool entered = false;
         for (size_t i=0; i<gridEntryP.size(); i++)
         {
            if (t>=nextEntryTime[i])
            {
               entered = true;
               nextEntryTime[i] += entryRate[i];

               // Check if there is room on inactive list for reassignment / entry
               // Else throw

               if (inactiveParticles.size()==0)
               {
                  std::string msg("Not enough inactive particles to support entry conditions.");
                  logger.error(msg.c_str());
                  throw std::runtime_error(msg);
               }
               else
               {
                  auto iter = inactiveParticles.begin();
                  size_t j = *iter;
                  inactiveParticles.erase(iter);
                  active[j] = true;
                  p[j] = gridEntryP[i];
                  q[j] = gridEntryQ[i];
                  mass[j] = gridEntryMass[i];

                  updateInactiveStats();

               }
            }
         }

         if (entered)
            logger.info("Particles entering at time {}.",t);
      }
      
   private:

      void updateInactiveStats (void)
      {
         size_t sizeNow = inactiveParticles.size();
         if (sizeNow > maxInactive)
            maxInactive = sizeNow;
         if (sizeNow < minInactive)
            minInactive = sizeNow;         
      }
      
      WindTunnel::Exit exitType;     
      double exitDistance;

      container_type& p;
      container_type& q;
      scalar_type& mass;
      mask_type& active;

      // Single set/list/vector of points to enter for boundary conditions
      // derived from grids.
      container_type gridEntryP;
      container_type gridEntryQ;
      scalar_type gridEntryMass;

      // Book keeping for grid entries
      // The last time each point entered.
      scalar_type nextEntryTime;
      // The rate or cadence at which points entered (units of time)
      scalar_type entryRate;
      
      spdlog::logger& logger;

      std::set<size_t> inactiveParticles;
      size_t minInactive, maxInactive;
   };
} // end namespace doorstep
