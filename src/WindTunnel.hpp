#pragma once

#include "System.hpp"

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
                    scalar_type& inputMass, mask_type& inputActive) :
            p(inputP), 
            q(inputQ),
            mass(inputMass),
            active(inputActive),
            exitType(WindTunnel::Exit::UNSPECIFIED),
            exitDistance(0.)
         {};

         void setExitAsDistance(double distance)
         {
            exitType = WindTunnel::Exit::DISTANCE_FROM_ORIGIN;
            exitDistance = distance;
         }

         void detectExits(void)
         {
            if (exitType == Exit::DISTANCE_FROM_ORIGIN)
            {
               for( size_t i=0; i<q.size() ; ++i )
               {
                  if (abs(q[i])>exitDistance)
                     active[i]=false;
               }
            }
         }

   private:

      WindTunnel::Exit exitType;     
      double exitDistance;

      container_type& p;
      container_type& q;
      scalar_type& mass;
      mask_type& active;

   };
} // end namespace doorstep
