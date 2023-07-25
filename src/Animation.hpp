#pragma once

#include <string>
#include <filesystem>
#include <exception>
#include <vector>

#include "Utils.hpp"
#include "Configuration.hpp"

namespace doorstep
{
   class Animation
   {
      public:

         enum Format
         {
            GIF,
            AVIMPG4,
            WEBP,
            UNSPECIFIED
         };

         Animation(void) : animationOutputPath("."),
                           format(Animation::Format::UNSPECIFIED),
                           openViewer(false)
         {}

         Animation(const AnimationConfiguration& acfg,
                   const ImageStream& is)
         {
            configure(acfg, is);
         }

         void configure(const AnimationConfiguration& acfg,
                        const ImageStream& is)
         {
            animationOutputPath = acfg.animationOutputPath;
            outputFilename = acfg.outputFilename;
            openViewer = acfg.openViewer;

            std::filesystem::path outputfn = animationOutputPath / outputFilename;
            std::string thisExtension = outputfn.extension().generic_string();

            tolower(thisExtension);

            if (thisExtension.compare(".gif")==0)
               format = Animation::Format::GIF;
            else if (thisExtension.compare(".avi")==0)
               format = Animation::Format::AVIMPG4;
            else if (thisExtension.compare(".webp")==0)
               format = Animation::Format::WEBP;
            else
               format = Animation::Format::UNSPECIFIED;

            // Set frame output C specified (for output and input to ffmpet)
            // and wildcard (for globbing)
            frameCspec = is.c_filespec;
            frameWildcard = is.file_wildcard;
            frameInputDir = is.imageOutputPath;

         }

         std::string getFormatDescription(void)
         {
            std::string desc;

            switch (format)
            {
               case GIF :
                  desc = "Animated GIF";
                  break;

               case AVIMPG4 :
                  desc = "AVI file, MPEG-4 encoding";
                  break;

               case WEBP :
                  desc = "Lossless WEBP animation";
                  break;

               default :
                  desc = "Unspecified";
            }

            return desc;
         }

         int generate(void)
         {
            if (format == Format::UNSPECIFIED)
               throw std::runtime_error("Cannot generate animation, format unspecified");

            std::string command;
            std::filesystem::path inputPath = frameInputDir / frameCspec;
            std::filesystem::path outputPath = animationOutputPath / outputFilename;

            if (format == Format::WEBP)
            {
               command = "ffmpeg -i " + inputPath.generic_string() +
                         " -vcodec libwebp -lossless 1 -loop 0 " +
                         outputPath.generic_string();
            }

            if (format == Format::AVIMPG4)
            {
               command = "ffmpeg -i " + inputPath.generic_string() +
                         " -vcodec mpeg4 -b:v 2M " + outputPath.generic_string();
            }

            if (format == Format::GIF)
            {
               inputPath = frameInputDir / frameWildcard;
               command = "convert " + inputPath.generic_string() +
                         " " + outputPath.generic_string();
            }
            
            return system(command.c_str());
         }


        int show(void)
        {
           std::filesystem::path animationPath = animationOutputPath / outputFilename;

           std::string command;

           #ifdef _WIN32
              command = "start " + animationPath.generic_string();
              return system(command.c_str());
           #else
	      if (format == Format::GIF)
	      {
	       std::string command = "eog " + animationPath.generic_string(); // Eye of Gnome

	      }

	       if (format == Format::AVIMPG4)
		  command = "vlc " + animationPath.generic_string();

	       if (format == Format::WEBP)
		  command = "firefox " + animationPath.generic_string();

	       return system(command.c_str());
	    #endif          
        } // end show() function
           
      bool openViewer;
      
      private:
         std::filesystem::path animationOutputPath;
         std::string outputFilename;

         std::string frameCspec;
         std::string frameWildcard;
         std::filesystem::path frameInputDir;
      
         Animation::Format format;
   };
}
