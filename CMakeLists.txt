# CMakeLists.txt for the "kloe-simu" package. It creates a library with a dictionary and a main program.
# If ROOT is not installed in a default system location you need to tell CMake where to find it.
# Sourcing `thisroot.sh` already sets the required environment variables.
# Otherwise, you'll have to tell the build system where to look for ROOT,
# e.g. by passing `-DROOT_DIR="/path/to/root/installation` at cmake configuration time.

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(kloe-simu LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_PLATFORM_INDEPENDENT_CODE ON)

# Locate the ROOT package and define a number of useful targets and variables.
find_package(ROOT REQUIRED COMPONENTS RIO Net)

# Locate EDep-sim
find_package(EDepSim)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-variable -Wno-unused-parameter")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
ROOT_GENERATE_DICTIONARY(G__Struct struct.h LINKDEF include/Linkdef.h)

# set output directory
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Creates a libStruct shared library with a generated dictionary.
add_library(Struct SHARED G__Struct.cxx)
target_link_libraries(Struct PUBLIC ROOT::RIO ROOT::Net)

# Creates a libUtils shared library
add_library(Utils SHARED src/utils.cpp)
target_link_libraries(Utils PUBLIC ROOT::RIO ROOT::Net EDepSim::edepsim_io)

# Creates Digitize executable.
add_executable(Digitize src/digitization.cpp)
target_link_libraries(Digitize Struct Utils)

# Creates Reconstruct executable.
add_executable(Reconstruct src/reconstruction.cpp)
target_link_libraries(Reconstruct Struct Utils)

# Creates Analyze executable.
add_executable(Analyze src/analysis.cpp)
target_link_libraries(Analyze Struct Utils ROOT::EG)

# Creates Display executable.
add_executable(Display src/display.cpp)
target_link_libraries(Display Struct Utils)

# Creates FasctCheck executable.
add_executable(FasctCheck src/fastcheck.cpp)
target_link_libraries(FasctCheck Struct Utils)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "default install path" FORCE )
endif()

# Copy setup.sh configuration file
configure_file(setup.sh "${CMAKE_INSTALL_PREFIX}/setup.sh" COPYONLY)

install(TARGETS Utils Struct Digitize Reconstruct Analyze Display FasctCheck
        RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
        
install(
   DIRECTORY "${CMAKE_BINARY_DIR}/."
   DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
   FILES_MATCHING PATTERN "*.pcm"
   PATTERN "lib" EXCLUDE
   PATTERN "bin" EXCLUDE
   PATTERN "CMakeFiles" EXCLUDE)

install(
  DIRECTORY include 
  DESTINATION "${CMAKE_INSTALL_PREFIX}"
  PATTERN "Linkdef.h" EXCLUDE)