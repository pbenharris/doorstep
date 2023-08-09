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
#include <string>
#include <limits>

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
#include "System.hpp"
#include "CelestialBody.hpp"
#include "Utils.hpp"
#include "Animation.hpp"
#include "Configuration.hpp"
#include "BodyDistribution.hpp"

using namespace std;
using namespace boost::numeric::odeint;
using namespace matplot;
using json = nlohmann::json;
using namespace doorstep;

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
   mask_type &active; // to debug detecting conditions (e.g., DM within bodies)
   double gravitational_constant;

   nbody_system_momentum( const scalar_type &inputMass,
                          const scalar_type &inputRadius,
                          mask_type &inputActive,
                          double inputG ) :
      mass(inputMass), radius(inputRadius), active(inputActive),
      gravitational_constant (inputG)
   { }

   void simple(const container_type &q , container_type &dpdt )
   {
      const size_t n = q.size();
      for( size_t i=0 ; i<n ; ++i )
      {
         active[i] = true;
         dpdt[i] = 0.0; // Ah, so dpdt is force. So, p is momentum.  
         for( size_t j=0 ; j<i ; ++j ) // The canonical form doesn't have velocity as a state variable, just positions/coordinates
         {
            point_type diff = q[j] - q[i];
            double d = abs( diff );
            double acceleration = gravitational_constant * mass[i] * mass[j] / d / d;

            // if active i and j...
            if (radius[j] > d)
            {
               acceleration = gravitational_constant * mass[i] * mass[j] * d / radius[j] / radius[j] / radius[j];
            }
            
            if (d > std::numeric_limits<double>::epsilon())
               diff *= acceleration / d; // Nice choice to reuse the difference variable for force
            else
               diff = 0;

            dpdt[i] += diff;
            dpdt[j] -= diff;
         } // inner loop over combinations of points
      } // outer loop over combinatios of points
   } // End simple function to loop through all points

   void operator()(const container_type &q , container_type &dpdt )
   {
      simple(q,dpdt);
   } // end operator()
};


// relic from integrating with first order rk
// Keep until first order rk is folded back in as an option
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
   // std::filesystem::path outputPath;
   // std::string framePattern;
   // FrameOutputConfiguration frameConfig;
   size_t numberCB;
   std::vector<int> frameCount;
   
   std::map<std::string, ImageStream> &imageStreamMap;

   scalar_type &radius; // for graphing. Body radius


   //mutable list<double> t_hist;
   //mutable list<state_type> x_hist;

   // Wind tunnel function
   bool windtunnel; // Are we enforcing the wind tunenl
   double windTunnelExitDist; // A fixed distance in the scale beyond which a dm point is reused.
   mask_type &active; // Outside windtunnel means inactive
   
   streaming_observer( std::ostream &out,
                       std::map<std::string,ImageStream> &inputImageStreamMap,
                       bool inputUseWindtunnel, double inputWindtunnelExit,
                       //std::filesystem::path iOutPath,
                       //FrameOutputConfiguration& ifc,
                       //std::string iFramePattern,
                       scalar_type &iRadius, mask_type &iActive,
                       size_t iNumberCelestialBodies) :
      m_out( out ), imageStreamMap(inputImageStreamMap), 
      radius(iRadius), active(iActive),
      numberCB(iNumberCelestialBodies),
      frameCount(inputImageStreamMap.size(),0),
      windtunnel(inputUseWindtunnel),
      windTunnelExitDist(inputWindtunnelExit)      
   {
   }

   template< class state_type >
   void operator()( const state_type &pq , double t )
   {
      container_type &x = pq.first;
      int n_dm = x.size() - numberCB;      

      // Process state output
      m_out << t;
      for( size_t i=0 ; i<x.size() ; ++i ) m_out << "\t" << x[i];
      m_out << "\n";

      // Wind tunnel processing
      // First tag those dm that left the windtunnel
      size_t n_dm_total = x.size() - numberCB; // Total number of dm particles
      size_t n_dm_active=0; // count of dm particles that are active
      size_t n_dm_inactive; // count of dm particles that are inactive
      if (windtunnel)
      {
         for( size_t i=0 ; i<x.size() ; ++i )
         {
            if (abs(x[i])>windTunnelExitDist)
               active[i]=false;
            else
            {
               if (i>=numberCB)
               {
                  active[i]=true;
                  n_dm_active++;
               }
            }
         }
      }
      n_dm_inactive = n_dm_total - n_dm_active;

      // Map coordinates from state into vectors for plotting
      size_t np = x.size();
      vector<double> xc(np), yc(np); // c is for celestial body
      auto xc_iter = xc.begin();
      auto yc_iter = yc.begin();

      // Output each stream
      int streamIndex = 0;
      for (auto i = imageStreamMap.begin();
           i != imageStreamMap.end();
           i++)
      {
         double length_scale = i->second.length_scale;

         matplot::cla();
         axis(square);
         axis({i->second.axis_min, i->second.axis_max,
               i->second.axis_min, i->second.axis_max});

         // Draw the celestial bodies.
         for (size_t i=0;i<numberCB;i++)
         {
         
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

         vector<double> x_dm_a(n_dm_active), y_dm_a(n_dm_active); // X,Y for actives
         vector<double> x_dm_i(n_dm_inactive), y_dm_i(n_dm_inactive); // X,Y for inactives
         size_t idx_i=0; // counter into inactives
         size_t idx_a=0; // counter into actives
         // Loop to sort through active from inactive into two vectors
         for (size_t i=0; i<n_dm_total; i++)
         {
            if (active[i+numberCB])
            {
               x_dm_a[idx_a] = x[i+numberCB][0] / length_scale;
               y_dm_a[idx_a++] = x[i+numberCB][1] / length_scale;
            }
            else
            {
               x_dm_i[idx_i] = x[i+numberCB][0] / length_scale;
               y_dm_i[idx_i++] = x[i+numberCB][1] / length_scale;
            }
         }

         // Plot actives
         auto pa = scatter(x_dm_a, y_dm_a, 2.0);
         pa->marker_color({0.0f, 0.5f, 0.5f});
         pa->marker_face_color({0.f, 0.5f, 0.5f});

         // Plot inactives
         auto pi = scatter(x_dm_i, y_dm_i, 2.0);
         pi->marker_color({1.0f, 0.0f, 0.0f});
         pi->marker_face_color({1.0f, 0.0f, 0.0f});
         
         // Time in title
         char timelabel[60];
         sprintf(timelabel,"t=%.2lf",t);
         title(timelabel);

         // Output active value in title
         //char activeLabel[60];
         //sprintf(activeLabel,"%g",active[1]);
         //title(activeLabel);

         char fname[255];
         sprintf(fname,i->second.c_filespec.c_str(),frameCount[streamIndex++]++);
         std::filesystem::path outputFramePath = i->second.imageOutputPath / fname;

         save(outputFramePath.generic_string());
      }

      // Control function to kill and redistribute, a "wind tunnel"
      if (windtunnel)
         for (size_t i=numberCB; i<(n_dm+numberCB); i++)
         {
            double alimit = windTunnelExitDist;

            // Check for distance from origin
            double distx = abs(x[i][0]);
            double disty = abs(x[i][1]);
            double distz = abs(x[i][2]);

            if (distx > alimit) active[i]=false;
            if (disty > alimit) active[i]=false;
            if (distz > alimit) active[i]=false;

         } // end control function
   } // end operator()
    
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
      
      spdlog::logger logger("multi_log", {console_sink, file_sink});
      logger.set_level(spdlog::level::debug);
      
      logger.info("DOORSTEP version {}.{}",
                   DOORSTEP_MAJOR_VERSION, DOORSTEP_MINOR_VERSION);

      string configfileName = "doorstep.json";
      if (argc ==2)
         configfileName = argv[1];

      logger.info("Reading configuration file {}", configfileName.c_str());

      ConfigurationFile cf(configfileName, logger);
      RunConfiguration rc = cf.getRunConfiguration();

      BodyDistribution bd(rc, logger);

      size_t n_body = bd.bodyCount();

      scalar_type mass2 (n_body, 0.);
      scalar_type radius2 (n_body, 0.01);
      container_type p2(n_body, 0.), q2(n_body, 0.);

      scalar_type mass = bd.getMass();
      scalar_type radius = bd.getRadius();      
      container_type p = bd.getP();
      container_type q = bd.getQ();

      mask_type active(n_body, true);

      logger.info("Simulation configured for {} bodies", n_body);
      // Clean up output files for animations
      for (auto it = rc.animationConfig.begin();
           it != rc.animationConfig.end();
           it++)
         cleanupFiles(it->animationOutputPath,
                      it->outputFilename, logger);
      
      // Clean up outputfiles for image streams (plots etc)
      for (auto it = rc.imageStreamMap.begin();
           it!=rc.imageStreamMap.end();
           it++)  
         cleanupFiles(it->second.imageOutputPath,
                      it->second.file_wildcard,
                      logger);

      //state_type x0(n*4); // 2D, position and velocity --> *4

      // Integration parameters
      double t0 = rc.initialTime;
      double t1 = rc.finalTime;
      double dt = rc.timeStep;

      ofstream stateHistFile("state_hist.txt");

      // Run integrator      
      streaming_observer so(stateHistFile,
                            rc.imageStreamMap,
                            rc.useWindTunnel,
                            rc.windTunnelExit,
                            //rc.animationConfig.animationOutputPath,
                            //rc.animationConfig.frameOutputConfig[0],
                            //animation.getFrameCspec(),
                            radius, active,
                            rc.celestialBody.size());
      
      //integrate_const( runge_kutta4<state_type>(), my_system, x0, t0, t1, dt, so );
      typedef symplectic_rkn_sb3a_mclachlan< container_type > stepper_type;
      //typedef runge_kutta4< container_type > stepper_type;

      integrate_const (
          stepper_type(),
          make_pair(nbody_system_coor(mass),
                    nbody_system_momentum(mass,
                                          radius,
                                          active,
                                          rc.gravitationalConstant) ), // System
          make_pair( boost::ref( q ), boost::ref( p ) ), // i.c.
          t0, t1, dt, boost::ref(so)
          );

      for (auto i=0; i<rc.animationConfig.size(); i++)
      {
         // Don't trust the stream is correct - check for input error
         auto iter = rc.imageStreamMap.find(rc.animationConfig[i].imageStreamName);
         if (iter==rc.imageStreamMap.end())
         {
            logger.error("Animation input stream {} not matched.",
                         rc.animationConfig[i].imageStreamName.c_str());
            throw std::runtime_error("Cannot locate animation input stream " +
                                     rc.animationConfig[i].imageStreamName);
         }
         Animation animation(rc.animationConfig[i],
                             iter->second);
         logger.info("Animation output format {}",animation.getFormatDescription().c_str());
         animation.generate();
         if (animation.openViewer)
            animation.show();
      }
   }
   catch (const spdlog::spdlog_ex &ex)
   {
      cerr << "Log init failed:" << ex.what() << endl;
   }
   catch (const std::runtime_error& ex)
   {
      cout << "Runtime exception: " << ex.what() <<  endl;
   }
   catch (const std::exception& ex)
   {
      cout << "Standard exception: " << ex.what() <<  endl;
   }
   catch(...)
   {
      cout << "Caught an unexpected error." << endl;
   }
   
   return(EXIT_SUCCESS);
}
