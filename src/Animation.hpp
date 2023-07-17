#pragma once

#include <string>

namespace doorstep
{
   struct Axis
   {
      bool applyLimits;
      double xScale;
      double yScale;
      double minX, minY, maxX, maxY;
   };
      
   class Animation
   {
      public:
         Animation(const Axis& inputAxis) :
            axis(inputAxis)
         {}
      
         void generateNextFrame(const std::vector<double>& x,
                                const std::vector<double>& y);
                    
      private:
         int frameNumber;
         Axis axis;
         double scale_reduction;
         std::string baseFilename;         
         int frameNumberDigits;
      
         std::filesystem::path outputDirectory;
         std::string outputFilename;
          
   };
}
