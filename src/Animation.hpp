#pragma once

#include <string>
#include <filesystem>
#include <exception>
#include <vector>

namespace doorstep
{
   struct FrameOutputConfiguration
   {
      double axis_min, axis_max;
      double length_scale;
   };

   struct AnimationConfiguration
   {
      AnimationConfiguration(): animationOutputPath(".")
      {
      }

      void setWorkingDirectory(const std::string& animationWorkingDirName)
      {
         std::filesystem::path homePath = getHome();
         animationOutputPath = homePath / animationWorkingDirName;
         if (!exists(animationOutputPath))
            throw std::runtime_error("Animation working directory " +
                                     animationOutputPath.generic_string() +
                                     "does not exist.");
      }

      std::vector<FrameOutputConfiguration> frameOutputConfig;
      std::filesystem::path animationOutputPath;
      std::string outputFilename;
   };

   class Animation
   {
      public:
      Animation()
      {}
                    
      private:

   };
}
