DOORSTEP
========

Dark Matter Object ORiented Simulation with TErrestrial Physics
*(c) Ben Harris, 2023*

DOORSTEP is a configurable N-body simulation tailored to physics
of dark matter interacting with the Earth. The code is modern C++, with
features as late as C++17. 
Resources

The project website is at darkmatteratourdoostep.com

## Build notes

DOORSTEP was written and tested in
Ubuntu 22.04 LTS. The first version of DOORSTEP was build using these
packages

 - gcc 11.3
 - cmake 3.22
 - ImageMagick 6
 - gnuplot 5.4 patchlevel 2
 - Matplot++
 - ffmpeg 4.4.2

DOORSTEP is built using CMake. Notes about building under different platforms
are below.

### Windows build

Windows has to be told where Boost is installed, because there is no
standard place to install libraries on that operating system. This sequence
of command will build DOORSTEP in Windows.

set Boost_ROOT=c:\Users\my_user_name\boost_1_82_0
set Matplot++_DIR="c:\Program Files (x86)\Matplot++"
cd doorstep\src
mkdir build
cd build
cmake .. -B .