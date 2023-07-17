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

}
