#pragma once

#include <filesystem>

namespace doorstep
{
   
   std::filesystem::path getHome(void)
   {
   #ifndef _WIN32
      // Tested in Windows 10
      return std::filesystem::path(getenv("HOME"));
   #else
      return std::filesystem::path(getenv("USERPROFILE"));
   #endif
      return std::filesystem::current_path();
   }
   
   void cleanupFiles(std::filesystem::path workingDir, std::string frameWildcard, spdlog::logger& logger)
   {
      std::filesystem::path inputFilespec = workingDir / frameWildcard;
      logger.info("Removing files with specification {} ", inputFilespec.c_str());

      int allRemoved = 0;
      for (auto& de : glob::glob(inputFilespec.generic_string()))
      {
         int thisRemoved = remove_all(de);
         allRemoved += thisRemoved;
      }
      logger.info("Removed {} files.",allRemoved);
      return;
   }

}
