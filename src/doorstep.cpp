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

typedef point< double , 3 > point_type;
typedef std::vector< point_type > container_type;
typedef std::vector< double > mass_type;

using namespace std;
using namespace boost::numeric::odeint;
using namespace matplot;
using json = nlohmann::json;

// Coordinate function
struct nbody_system_coor
{
    const mass_type &m_masses;

    nbody_system_coor( const mass_type &masses ) : m_masses( masses ) { }

    void operator()( const container_type &p , container_type &dqdt ) const
    {
       const size_t n = p.size();
       for( size_t i=0 ; i<n ; ++i )
          dqdt[i] = p[i] / m_masses[i];
    }
};

// Momentum function
struct nbody_system_momentum
{
   const mass_type &m_masses;
   double gravitational_constant;

   nbody_system_momentum( const mass_type &masses, double inputG ) :
      m_masses( masses ), gravitational_constant (inputG)
   { }

   void operator()( const container_type &q , container_type &dpdt ) const
   {
      const size_t n = q.size();
      for( size_t i=0 ; i<n ; ++i )
      {
         dpdt[i] = 0.0; // Ah, so dpdt is force. So, p is momentum.  
         for( size_t j=0 ; j<i ; ++j ) // The canonical form doesn't have velocity as a state variable, just positions/coordinates
         {
            point_type diff = q[j] - q[i];
            double d = abs( diff );
            diff *= ( gravitational_constant * m_masses[i] * m_masses[j] / d / d / d ); // Nice choice to reuse the difference variable for force
            dpdt[i] += diff;
            dpdt[j] -= diff;
         } // inner loop over combinations of points
      } // outer loop over combinatios of points
   } // end parenthesis operator
};

point_type center_of_mass( const container_type &x , const mass_type &m )
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
    
    //mutable list<double> t_hist;
    //mutable list<state_type> x_hist;
   
   streaming_observer( std::ostream &out, std::filesystem::path iOutPath,
                       std::string iFramePattern) :
      m_out( out ), outputPath(iOutPath), framePattern(iFramePattern)
   {
   }

   template< class state_type >
   void operator()( const state_type &pq , double t ) const
   {
      container_type &x = pq.first;
      
      m_out << t;
      for( size_t i=0 ; i<x.size() ; ++i ) m_out << "\t" << x[i];
      m_out << "\n";

      // Copy coordinates from state into separate vectors
      size_t np = x.size();
      vector<double> xc(np), yc(np);
      auto xc_iter = xc.begin();
      auto yc_iter = yc.begin();

      for (size_t i=0;i<np;i++)
      {
         // These get compile errors
         ///   *xc_iter++ = x_iter++;
         //*yc_iter++ = x_iter++;
         xc[i] = x[i][0];
         yc[i] = x[i][1];
      }

      scatter(xc, yc, 1);

      //auto r1 = rectangle(0,0,1,1);
      //r1->color("red");

      char timelabel[60];
      sprintf(timelabel,"t=%g",t);
      //text(-4,-4,timelabel);
      title(timelabel);
      
      static int framecount = 0;
      char fname[60];
      sprintf(fname,framePattern.c_str(),framecount++);
      std::filesystem::path outputFramePath = outputPath / fname;
      save(outputFramePath.generic_string());
          
      axis({-5, 5, -5, 5});
      //axis({0, 1, 0 , 1});
      //axis(matplot::equal);

      // show();
   }
    
};

std::filesystem::path getHome(void)
{
#ifndef _WIN32
   return std::filesystem::path(getenv("HOME"));
#else
   return std::filesystem::path(getenv("USERPROFILE"));
#endif
   return std::filesystem::current_path();
}

// Read frames from working dir, generates movie there
// Assumes ImageMagick or similar is installed
// Assumes working directory exists
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

void cleanupFrames(std::filesystem::path workingDir, std::string frameWildcard)
{
   std::filesystem::path inputFilespec = workingDir / frameWildcard;
   cout << "Removing frame image files with specification: "
        << inputFilespec << endl;

   int allRemoved = 0;
   for (auto& de : glob::glob(inputFilespec))
   {
      int thisRemoved = remove_all(de);
      allRemoved += thisRemoved;
   }
   cout << "Removed " << allRemoved << " prior frame image files." << endl;
   return;
}

void viewMovie(std::filesystem::path workingDir, std::string movieFilename)
{
   std::filesystem::path inputFilespec = workingDir / movieFilename;
   std::string command = "eog " + inputFilespec.generic_string();
   cout << "Spawning view." << endl << flush;
   system(command.c_str());
   return;
}

struct RunConfiguration
{
   int numberBodies;
   double gravitationalConstant;
};

class ConfigurationFile
{
   public:
      ConfigurationFile(const std::string& configFilename, spdlog::logger& mlogger)
      {
         std::list<std::filesystem::path> pathList;
         pathList.push_back(configFilename);
         std::filesystem::path home = getHome();
         mlogger.debug("Configuration has detected home dir at {}",home.generic_string().c_str());
         pathList.push_back(getHome() / configFilename);
         pathList.push_back(getHome() / "doorstep" / "src" / configFilename);

         bool foundConfig = false;
         std:filesystem::path configPathUsing;

         for (auto it = pathList.begin(); it!=pathList.end(); it++)
         {
            ifstream testOpen(*it);
            if (!testOpen.fail())
            {
               foundConfig = true;
               configPathUsing = *it;
               mlogger.info("Found configuration file {}",(*it).generic_string().c_str());
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
         data = nlohmann::json::parse(configFile);     
         mlogger.info("Read configuration from {}",
                      configPathUsing.generic_string().c_str());
      }

      RunConfiguration getRunConfiguration(void) const
      {
         RunConfiguration rc;
         rc.numberBodies = data["number_bodies"].template get<int>();
         rc.gravitationalConstant = data["gravitational_constant"].template get<int>();
         return rc;
      }

private:
   nlohmann::json data;
};
// Version info
const int DOORSTEP_MAJOR_VERSION = 0;
const int DOORSTEP_MINOR_VERSION = 1;

// ------  Main
int main()
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

      ConfigurationFile cf("doorstep.json", logger);
      RunConfiguration rc = cf.getRunConfiguration();
      
      mass_type masses (rc.numberBodies, 0.);
      container_type p(rc.numberBodies, 0.), q(rc.numberBodies, 0.);

      //state_type x0(n*4); // 2D, position and velocity --> *4

      
      // Populate the initial state: mass and position
      //default_random_engine defEngine(time(0));
      default_random_engine rEngine;
      uniform_real_distribution<double> initialDist(0,1);
      size_t i;
      for (i=0; i<rc.numberBodies; i++)
      {
         masses[i] = 1.0;
         q[i][0] = initialDist(rEngine);
         q[i][1] = initialDist(rEngine);

         double x_plus = initialDist(rEngine);
         double x_minus = initialDist(rEngine);
         double y_plus =  initialDist(rEngine);
         double y_minus = initialDist(rEngine);

         double scale = 1./10;
         p[i][0] = scale*(x_plus - x_minus);
         p[i][1] = scale*(y_plus - y_minus);
      }

      point_type qmean = center_of_mass( q , masses );
      point_type pmean = center_of_mass( p , masses );
      for( size_t i=0 ; i<rc.numberBodies ; ++i )
      {
         q[i] -= qmean ;
         p[i] -= pmean;
      }

      for( size_t i=0 ; i<rc.numberBodies ; ++i ) p[i] *= masses[i];

      // Integration parameters
      double t0 = 0.0;
      double t1 = 40.0;
      double dt = 0.1;

      ofstream stateHistFile("state_hist.txt");
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

      // Clean up old frames that use the wildcard
      string frameWildcard = "frame*.png";
      string frameCspec = "frame%04d.png";
      string movieFilename = "movie.gif";

      cleanupFrames(animationOutputPath, frameWildcard);

      // Run integrator
      streaming_observer so(stateHistFile, animationOutputPath, frameCspec);   
      //integrate_const( runge_kutta4<state_type>(), my_system, x0, t0, t1, dt, so );
      typedef symplectic_rkn_sb3a_mclachlan< container_type > stepper_type;
      //typedef runge_kutta4< container_type > stepper_type;

      integrate_const (
          stepper_type(),
          make_pair(nbody_system_coor(masses),
                    nbody_system_momentum(masses,
                                          rc.gravitationalConstant) ), // System
          make_pair( boost::ref( q ), boost::ref( p ) ), // i.c.
          t0, t1, dt, boost::ref(so)
          );

      generateGIF(animationOutputPath, frameWildcard, movieFilename);

      viewMovie(animationOutputPath, movieFilename);
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
