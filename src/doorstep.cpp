// DOORSTEP
// C++ simulation of dark matter motion in a two body system
//
// Ben Harris
//
// Website: www.darkmatteratourdoostep.com
// About this code: see the README.md and LICENSE
// Plans for this code: see TODO.txt

// C++ standard includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <filesystem>
#include <vector>
#include <list>
#include <exception>

// C includes (where possible the C++ standard version)
#include <cstdlib>

// External libraries
#include <boost/numeric/odeint.hpp>       // odeint function definitions
#include <matplot/matplot.h>
#include <glob.hpp>
#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/ostream_sink.h>

// DOORSTEP includes
#include "point_type.hpp" // from ODEINT examples
#include "CelestialBody.hpp"
#include "Utils.hpp"
#include "Animation.hpp"

typedef point< double , 3 > point_type;
typedef std::vector< point_type > container_type;
typedef std::vector< double > scalar_type;

using namespace std;
using namespace boost::numeric::odeint;
using namespace matplot;
using json = nlohmann::json;
using namespace doorstep;

int generateWEBP(std::filesystem::path workingDir,
                 std::string frameCspec,
                 std::string outputName)
{
   std::string command = "ffmpeg ";
   
   std::filesystem::path inputPath = workingDir / frameCspec;
   std::filesystem::path outputPath = workingDir / outputName;

   command += " -i " + inputPath.generic_string() + " -vcodec libwebp -lossless 1 -loop 0 " + outputPath.generic_string();
   
   //Move echo out of this to configuration confirmation?
   cout << "Reading frames from dir: " << inputPath << endl;
   cout << "Writing movie to file: " << outputPath << endl;
   cout << "Command is: " << command << endl;
   
   int success=system(command.c_str());
   
   cout << "Done." << endl << flush;
   return success;
}

int generateAVI(std::filesystem::path workingDir,
                std::string frameCspec,
                std::string outputAVIname)
{
   std::string command = "ffmpeg ";
   
   std::filesystem::path inputPath = workingDir / frameCspec;
   std::filesystem::path outputPath = workingDir / outputAVIname;

   command += " -i " + inputPath.generic_string() + " -vcodec mpeg4 -b:v 2M " + outputPath.generic_string();
   
   //Move echo out of this to configuration confirmation?
   cout << "Reading frames from dir: " << inputPath << endl;
   cout << "Writing movie to file: " << outputPath << endl;
   cout << "Command is: " << command << endl;
   
   int success=system(command.c_str());
   
   cout << "Done." << endl << flush;
   return success;
}

int generateGIF(std::filesystem::path workingDir, std::string framewildcard,
                std::string outputGIFname)
{
   std::string command = "convert ";
   
   std::filesystem::path inputPath = workingDir / framewildcard;
   std::filesystem::path outputPath = workingDir / outputGIFname;
   command += inputPath.generic_string() + " " + outputPath.generic_string();
   
   //Move echo out of this to configuration confirmation?
   cout << "Reading frames from dir: " << inputPath << endl;
   cout << "Writing movie to file: " << outputPath << endl;
   
   int success=system(command.c_str());
   
   cout << "Done." << endl << flush;
   return success;
}


void viewMovie(std::filesystem::path workingDir, std::string movieFilename)
{
   std::filesystem::path inputFilespec = workingDir / movieFilename;
#ifdef _WIN32
   std::string command = "start " + inputFilespec.generic_string(); 
#else
   string thisExtension = inputFilespec.extension().generic_string();

   std::string command = "eog " + inputFilespec.generic_string();

   int cmp =thisExtension.compare(".avi"); 
   if (cmp==0)
      command = "vlc " +  inputFilespec.generic_string();

   cmp = thisExtension.compare(".webp");
   if (cmp==0)
      command = "firefox " +  inputFilespec.generic_string();
#endif

   cout << "Spawning view." << endl << flush;
   system(command.c_str());
   return;
}

struct GridConfiguration
{
   double mass;
   size_t nx, ny, nz;
   double distance;
   double ic_x, ic_y, ic_z;
   double ic_vx, ic_vy, ic_vz;
};

struct RunConfiguration
{
   int numberBodies;
   double gravitationalConstant;
   double initialTime;
   double finalTime;
   double timeStep;
   std::vector<CelestialBody> celestialBody;
   AnimationConfiguration animationConfig;
   std::vector<GridConfiguration> gridConfig;
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
         std:filesystem::path configPathUsing;

         for (auto it = pathList.begin(); it!=pathList.end(); it++)
         {
            ifstream testOpen(*it);
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

         ifstream configFile(configPathUsing.generic_string());
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
         }

         // Configuration for animation output - option here but not in main()
         if (jsonData.contains("animation"))
         {
            rc.animationConfig.setWorkingDirectory(jsonData["animation"]["working_directory"]);
            rc.animationConfig.outputFilename = jsonData["animation"]["output_filename"];
            rc.animationConfig.frameOutputConfig.resize(jsonData["animation"]["frame_outputs"].size());

            for (auto i=0; i<rc.animationConfig.frameOutputConfig.size();i++)
            {
               nlohmann::json frameJsonData = jsonData["animation"]["frame_outputs"][i];
               rc.animationConfig.frameOutputConfig[i].axis_min = frameJsonData["axis_min"];
               rc.animationConfig.frameOutputConfig[i].axis_max = frameJsonData["axis_max"];
               rc.animationConfig.frameOutputConfig[i].length_scale = frameJsonData["length_scale"];
             }
         }

         // Grid of dm to use for testing
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
};


// Coordinate function
struct nbody_system_coor
{
    const scalar_type &mass;

    nbody_system_coor( const scalar_type &inputMass ) : mass( inputMass ) { }

    void operator()( const container_type &p , container_type &dqdt ) const
    {
       const size_t n = p.size();
       for( size_t i=0 ; i<n ; ++i )
          dqdt[i] = p[i] / mass[i];
    }
};

// Momentum function
struct nbody_system_momentum
{
   const scalar_type &mass;
   const scalar_type &radius; // for calculating force within bodies
   scalar_type &metric; // to debug detecting conditions (e.g., DM within bodies)
   double gravitational_constant;

   nbody_system_momentum( const scalar_type &inputMass,
                          const scalar_type &inputRadius,
                          scalar_type &inputMetric,
                          double inputG ) :
      mass(inputMass), radius(inputRadius), metric(inputMetric),
      gravitational_constant (inputG)
   { }

   void operator()( const container_type &q , container_type &dpdt )
   {
      const size_t n = q.size();
      for( size_t i=0 ; i<n ; ++i )
      {
         metric[i] =0.0;
         dpdt[i] = 0.0; // Ah, so dpdt is force. So, p is momentum.  
         for( size_t j=0 ; j<i ; ++j ) // The canonical form doesn't have velocity as a state variable, just positions/coordinates
         {
            point_type diff = q[j] - q[i];
            double d = abs( diff );
            double acceleration = gravitational_constant * mass[i] * mass[j] / d / d;
            
            //if (d<radius[i])
            //   metric[i] = 1.0;
            if (i==1)
            {
               if (radius[j] > d)
                  metric[i] = 1.0;
            }

            if (radius[j] > d)
            {
               acceleration = gravitational_constant * mass[i] * mass[j] * d / radius[j] / radius[j] / radius[j];
            }
            
            diff *=  acceleration / d ; // Nice choice to reuse the difference variable for force
            dpdt[i] += diff;
            dpdt[j] -= diff;
         } // inner loop over combinations of points
      } // outer loop over combinatios of points
   } // end parenthesis operator
};

point_type center_of_mass( const container_type &x , const scalar_type &m )
{
    const size_t n = x.size();
    double overall_mass = 0.0;
    point_type mean( 0.0 );
    for( size_t i=0 ; i<n ; ++i )
    {
        overall_mass += m[i];
        mean += m[i] * x[i];
    }
    if( !x.empty() ) mean /= overall_mass;
    return mean;
}

// Defining a state type
//typedef std::vector< double > state_type;
 
// System of DEs to be solved
/*void my_system( const state_type &x , state_type &dxdt , const double t )
{
   size_t n = x.size()/2;
   size_t i;
   
   for (i=0; i<n; i++) dxdt[i] = x[i+n];

   // zero force
   for (i=n; i<2*n; i++) dxdt[i]=0;
   
}
*/

struct streaming_observer
{
    std::ostream& m_out;
    std::filesystem::path outputPath;
    std::string framePattern;
    FrameOutputConfiguration frameConfig;
    size_t numberCB;

    scalar_type &radius; // for graphing. Body radius
    scalar_type &metric; // for drawing conditions (e.g., inside another body)

    //mutable list<double> t_hist;
    //mutable list<state_type> x_hist;
   
   streaming_observer( std::ostream &out, std::filesystem::path iOutPath,
                       FrameOutputConfiguration& ifc, std::string iFramePattern,
                       scalar_type &iRadius, scalar_type &iMetric,
                       size_t iNumberCelestialBodies) :
      m_out( out ), outputPath(iOutPath),
      frameConfig(ifc), framePattern(iFramePattern),
      radius(iRadius), metric(iMetric), numberCB(iNumberCelestialBodies)
   {
   }

   template< class state_type >
   void operator()( const state_type &pq , double t ) const
   {
      container_type &x = pq.first;
      container_type &v = pq.second;

      m_out << t;
      for( size_t i=0 ; i<x.size() ; ++i ) m_out << "\t" << x[i];
      m_out << "\n";

      // Copy coordinates from state into separate vectors
      size_t np = x.size();
      vector<double> xc(np), yc(np);
      auto xc_iter = xc.begin();
      auto yc_iter = yc.begin();

      double length_scale = frameConfig.length_scale;

      matplot::cla();
      axis(square);
      axis({frameConfig.axis_min, frameConfig.axis_max,
            frameConfig.axis_min, frameConfig.axis_max});

      // Draw the celestial bodies.
      for (size_t i=0;i<numberCB;i++)
      {
         // These get compile errors
         ///   *xc_iter++ = x_iter++;
         //*yc_iter++ = x_iter++;

         // Scaled for doorstep.json.gps
         
         xc[i] = x[i][0]/length_scale;
         yc[i] = x[i][1]/length_scale;
         double thisRadius = radius[i]/length_scale;
         auto e=ellipse(xc[i]-thisRadius,
                        yc[i]-thisRadius,
                        2*thisRadius,
                        2*thisRadius);
         e->fill("green");
         e->color("green");
      }
            
      // Draw the dark matter particles
      hold(on);
      size_t n_dm = x.size() - numberCB;
      vector<double> x_dm(n_dm), y_dm(n_dm);
      for (size_t i=0; i<n_dm; i++)
      {
         x_dm[i] = x[i+numberCB][0] / length_scale;
         y_dm[i] = x[i+numberCB][1] / length_scale;
      }
      scatter(x_dm, y_dm, 2);

      // Kill and redistribute (aka, the wind tunnel)
      for (size_t i=numberCB; i<(n_dm+numberCB); i++)
      {
         double alimit = 4;
         
         // Check for distance from origin
         double distx = abs(x[i][0]/length_scale);
         double disty = abs(x[i][1]/length_scale);
         double distz = abs(x[i][2]/length_scale);

         if (distx > alimit) x[i][0]=-x[i][0];
         if (disty > alimit) x[i][1]=-x[i][1];
         if (distz > alimit) x[i][2]=-x[i][2];
         
      }
      
      
      // Time in title
      char timelabel[60];
      sprintf(timelabel,"t=%g",t);
      title(timelabel);

      // Output metric value in title
      //char metricLabel[60];
      //sprintf(metricLabel,"%g",metric[1]);
      //title(metricLabel);

      static int framecount = 0;
      char fname[60];
      sprintf(fname,framePattern.c_str(),framecount++);
      std::filesystem::path outputFramePath = outputPath / fname;
  
      
      save(outputFramePath.generic_string());
      
      // show();
   }
    
};

// Version info
const int DOORSTEP_MAJOR_VERSION = 0;
const int DOORSTEP_MINOR_VERSION = 1;

// ------  Main
int main(int argc, char* argv[])
{
   try {
      auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
      console_sink->set_level(spdlog::level::debug);

      auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("doorstep-log.txt", true);
      file_sink->set_level(spdlog::level::debug);
      
      spdlog::logger logger("multi_sink", {console_sink, file_sink});
      logger.set_level(spdlog::level::debug);
      
      logger.info("DOORSTEP version {}.{}",
                   DOORSTEP_MAJOR_VERSION, DOORSTEP_MINOR_VERSION);

      string configfileName = "doorstep.json";
      if (argc ==2)
         configfileName = argv[1];

      logger.info("Reading configuration file {}", configfileName.c_str());

      ConfigurationFile cf(configfileName, logger);
      RunConfiguration rc = cf.getRunConfiguration();

      
      size_t n_body = rc.numberBodies;
      scalar_type mass (n_body, 0.);
      scalar_type radius (n_body, 0.01);
      scalar_type metric (n_body, 0.);
      container_type p(n_body, 0.), q(n_body, 0.);

      // Clean up old frames that use the wildcard
      string frameWildcard = "frame*.png";
      string frameCspec = "frame%04d.png";

      //cleanupFiles(workingDir, outputAVIname);
      cleanupFiles(rc.animationConfig.animationOutputPath, rc.animationConfig.outputFilename, logger);
      cleanupFiles(rc.animationConfig.animationOutputPath, frameWildcard, logger);

      //state_type x0(n*4); // 2D, position and velocity --> *4

      
      // Populate the initial state: mass and position
      //default_random_engine defEngine(time(0));
      default_random_engine rEngine;
      uniform_real_distribution<double> initialDist(-4,4);
      size_t i;
      for (i=0; i<n_body; i++)
      {
         if (i<rc.celestialBody.size())
         {
            mass[i] = rc.celestialBody[i].mass;
            radius[i] = rc.celestialBody[i].radius;
            size_t j;
            for (j=0;j<3;j++)
            {
               q[i][j] = rc.celestialBody[i].position[j];
               p[i][j] = rc.celestialBody[i].velocity[j];
            }
         }
         else // random distribution
         {
            mass[i] = 1e-9;
            radius[i] = 0.01;
            q[i][0] = initialDist(rEngine);
            q[i][1] = initialDist(rEngine);

            double x_plus = initialDist(rEngine);
            double x_minus = initialDist(rEngine);
            double y_plus =  initialDist(rEngine);
            double y_minus = initialDist(rEngine);

            double scale = 1./100;
            p[i][0] = scale*(x_plus - x_minus);
            p[i][1] = scale*(y_plus - y_minus);
         }
      }

      point_type qmean = center_of_mass( q , mass );
      point_type pmean = center_of_mass( p , mass );
      for( size_t i=0 ; i<n_body ; ++i )
      {
         q[i] -= qmean ;
         p[i] -= pmean;
      }

      for( size_t i=0 ; i<n_body ; ++i ) p[i] *= mass[i];

      // Integration parameters
      double t0 = rc.initialTime;
      double t1 = rc.finalTime;
      double dt = rc.timeStep;

      ofstream stateHistFile("state_hist.txt");

      /*
      filesystem::path homePath = getHome();
      filesystem::path animationOutputPath = homePath / "animation";

      if (exists(animationOutputPath))
         cout << "Outputing animation to " << animationOutputPath << endl;
      else
      {
         cerr << "Animation output directory " << animationOutputPath
              << " does not exist." << endl;
         exit(1);
      }
      */

      // Run integrator
      streaming_observer so(stateHistFile, rc.animationConfig.animationOutputPath,
                            rc.animationConfig.frameOutputConfig[0],
                            frameCspec,
                            radius, metric, rc.celestialBody.size());   
      //integrate_const( runge_kutta4<state_type>(), my_system, x0, t0, t1, dt, so );
      typedef symplectic_rkn_sb3a_mclachlan< container_type > stepper_type;
      //typedef runge_kutta4< container_type > stepper_type;

      integrate_const (
          stepper_type(),
          make_pair(nbody_system_coor(mass),
                    nbody_system_momentum(mass,
                                          radius,
                                          metric,
                                          rc.gravitationalConstant) ), // System
          make_pair( boost::ref( q ), boost::ref( p ) ), // i.c.
          t0, t1, dt, boost::ref(so)
          );

      //generateGIF(animationOutputPath, frameWildcard, string("test.gif"));
      //generateAVI(animationOutputPath, frameCspec, animationFilename);
      generateWEBP(rc.animationConfig.animationOutputPath, frameCspec,
                   rc.animationConfig.outputFilename);
      viewMovie(rc.animationConfig.animationOutputPath, rc.animationConfig.outputFilename);
   }
   catch (const spdlog::spdlog_ex &ex)
   {
      cout << "Log init failed: " << ex.what() << endl;
   }
   catch (const std::runtime_error& ex)
   {
      cout << "Runtime exception: " << ex.what() <<  endl;
   }
   
   return(EXIT_SUCCESS);
}
