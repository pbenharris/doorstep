#pragma once

#include <string>
#include <fstream>
#include <list>
#include <vector>
#include <map>
#include <filesystem>

#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>

#include "Utils.hpp"
#include "CelestialBody.hpp"

namespace doorstep {

   struct GridConfiguration
   {
      double mass;
      size_t nx, ny, nz;
      double distance;
      double ic_x, ic_y, ic_z;
      double ic_vx, ic_vy, ic_vz;
   };

   struct RandomConfiguration
   {
      size_t number_bodies;
   };
   
   struct ImageStream
   {
      double axis_min, axis_max;
      double length_scale;
      std::filesystem::path imageOutputPath;
      std::string c_filespec, file_wildcard;

      void setOutputDirectory(const std::string& workingDirName)
      {
         std::filesystem::path homePath = getHome();
         imageOutputPath = homePath / workingDirName;
         if (!exists(imageOutputPath))
            throw std::runtime_error("Image working directory " +
                                     imageOutputPath.generic_string() +
                                     "does not exist.");
      }
   };

   struct AnimationConfiguration
   {
      AnimationConfiguration(): animationOutputPath("."),
                                openViewer(false)
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
      std::filesystem::path animationOutputPath;
      std::string outputFilename;
      std::string imageStreamName;
      bool openViewer;
   };

   struct RunConfiguration
   {
      int numberBodies;
      double gravitationalConstant;
      double initialTime;
      double finalTime;
      double timeStep;
      bool useWindTunnel;
      double windTunnelExit;
      std::vector<CelestialBody> celestialBody;
      std::vector<GridConfiguration> gridConfig;
      std::vector<RandomConfiguration> randomConfig;
      std::map<std::string, ImageStream> imageStreamMap;
      std::vector<AnimationConfiguration> animationConfig;
   };

   class ConfigurationFile
   {
      public:
         ConfigurationFile(const std::string& configFilename, spdlog::logger& mlogger)
         {
            // Logic to find the configuration file
            // First place it is found is used. Which is found is logged.
            std::list<std::filesystem::path> pathList;
            pathList.push_back(configFilename); // represents the current working directory
            std::filesystem::path home = getHome();
            pathList.push_back(getHome() / configFilename); // represents home directory
            pathList.push_back(getHome() / "doorstep" / "src" / configFilename); // development directory off of home

            bool foundConfig = false;
        
            std::filesystem::path configPathUsing;

            for (auto it = pathList.begin(); it!=pathList.end(); it++)
            {
               std::ifstream testOpen(*it);
               if (!testOpen.fail())
               {
                  foundConfig = true;
                  configPathUsing = *it;
                  mlogger.debug("Found configuration file {}",(*it).generic_string().c_str());
               }
            }

            if (!foundConfig)
            {
               mlogger.error("Cannot find a configuration file {}",
                            configFilename.c_str() );
               throw(std::runtime_error("Cannot open configuration file"+
                                         configFilename));
            }  

            std::ifstream configFile(configPathUsing.generic_string());
            jsonData = nlohmann::json::parse(configFile);     
            mlogger.info("Read configuration from {}",
                         configPathUsing.generic_string().c_str());
         }

         RunConfiguration getRunConfiguration(void) const
         {
            RunConfiguration rc;

            // Required configurations defining simulation time frame and G
            rc.numberBodies = jsonData["number_bodies"];
            rc.gravitationalConstant = jsonData["gravitational_constant"];
            rc.initialTime = jsonData["initial_time"].template get<double>();
            rc.finalTime = jsonData["final_time"].template get<double>();
            rc.timeStep = jsonData["time_step"].template get<double>();

            // Optional use of wind tunnel
            if (jsonData.contains("wind_tunnel_exit"))
            {
               rc.useWindTunnel = true;
               rc.windTunnelExit = jsonData["wind_tunnel_exit"];
            }
            else
            {
               rc.useWindTunnel = false;
               rc.windTunnelExit = 0.0;
            }
                    
            // Optional configuration to define celestial bodies
            if (jsonData.contains("celestial_body"))
            {
               rc.celestialBody.resize(jsonData["celestial_body"].size());
               for (auto i=0; i<rc.celestialBody.size();i++) 
               {
                  nlohmann::json body_data=jsonData["celestial_body"][i];
                  CelestialBody thisBody;
                  thisBody.name = body_data["name"];
                  thisBody.radius = body_data["radius"];
                  thisBody.mass = body_data["mass"];

                  thisBody.position.resize(body_data["initial_position"].size());
                  size_t j;
                  for (j=0; j<thisBody.position.size(); j++)
                     thisBody.position[j]=body_data["initial_position"][j];

                  thisBody.velocity.resize(body_data["initial_velocity"].size());
                  for (j=0; j<thisBody.velocity.size(); j++)
                     thisBody.velocity[j]=body_data["initial_velocity"][j];

                  rc.celestialBody[i] = thisBody;
               }
            } // end check for celestial bodies in the simulation

            // Configuration for animation output - option here but not in main()
            if (jsonData.contains("animation"))
            {
               rc.animationConfig.resize(jsonData["animation"].size());
               for (auto i=0; i<rc.animationConfig.size();i++)
               {
                  AnimationConfiguration ac;
                  ac.setWorkingDirectory(jsonData["animation"][i]["working_directory"]);
                  ac.outputFilename = jsonData["animation"][i]["output_filename"];
                  ac.imageStreamName = jsonData["animation"][i]["image_stream"];
                  ac.openViewer = jsonData["animation"][i]["open_in_viewer"];
                  rc.animationConfig[i] = ac;
               }
            } // end check for animation output
 

            if (jsonData.contains("image_stream"))
            {
               for (auto i=0; i<jsonData["image_stream"].size();i++)
               {
                  ImageStream is;
                  is.axis_min = jsonData["image_stream"][i]["axis_min"];
                  is.axis_max = jsonData["image_stream"][i]["axis_max"];
                  is.length_scale = jsonData["image_stream"][i]["length_scale"];
                  is.c_filespec = jsonData["image_stream"][i]["c_filespec"];
                  is.file_wildcard = jsonData["image_stream"][i]["file_wildcard"];
                  is.setOutputDirectory(jsonData["image_stream"][i]["working_directory"]);

                  rc.imageStreamMap[jsonData["image_stream"][i]["name"]] = is;
               }
            } // end check for image streams

            // Use a grid to specify initial conditions
            if (jsonData.contains("grid"))
            {
               rc.gridConfig.resize(jsonData["grid"].size());
               for (size_t i=0; i<rc.gridConfig.size(); i++)
               {
                  nlohmann::json gcjo = jsonData["grid"][i];
                  GridConfiguration gc;

                  gc.mass     = gcjo["particle_mass"];
                  gc.nx       = gcjo["n_x"];
                  gc.ny       = gcjo["n_y"];
                  gc.nz       = gcjo["n_z"];
                  gc.distance = gcjo["distance"];
                  gc.ic_x     = gcjo["initial_position"][0];
                  gc.ic_y     = gcjo["initial_position"][1];
                  gc.ic_z     = gcjo["initial_position"][2];
                  gc.ic_vx    = gcjo["initial_velocity"][0];
                  gc.ic_vy    = gcjo["initial_velocity"][1];
                  gc.ic_vz    = gcjo["initial_velocity"][2];

                  rc.gridConfig[i]=gc;
               }
            }
            return rc;
         }

   private:
      nlohmann::json jsonData;
   }; // end class ConfigurationFile

} // end namespace doorstep
