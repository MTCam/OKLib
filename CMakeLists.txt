cmake_minimum_required(VERSION 3.2)
project(OKLib CXX C Fortran)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/Modules")
include(SetColors)

# force full path to libraries
if(NOT CMAKE_VERSION VERSION_LESS "3.3")
  cmake_policy(SET CMP0060 NEW)
endif()

#
# options for using CS-Tools
#
option(CS-Tools         "Enable CS Tools" OFF)

include(FortranCInterface)
FortranCInterface_HEADER(${PROJECT_BINARY_DIR}/include/FC.h MACRO_NAMESPACE "FC_")
FortranCInterface_VERIFY(CXX)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/include)


# Project-specific Options
OPTION(BUILD_STATIC "Build static libraries")
#OPTION(BUILD_PROVENANCE "Record repository provenance" ON)
#OPTION(ENABLE_MPI "Build with MPI Support" ON)
#OPTION(ENABLE_MPIP "Enables MPIP MPI-profiling library" OFF)
#OPTION(ENABLE_HDF5 "Build with HDF5 Support" ON)
OPTION(ENABLE_OPENMP "Build with OpenMP parallelism" OFF)
#OPTION(AUTODETECT_MPI "Automatically detect MPI" ON)
OPTION(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ON)
OPTION(STATIC_KERNELS "Build static library for Fortran computation kernels" OFF)
OPTION(ENABLE_COVERAGE "Enable coverage testing with gcov (GCC only)" OFF)
OPTION(ENABLE_KNL "Build for Intel Knights Landing architecture" OFF)
OPTION(ENABLE_PAPI "Build with support for PAPI hardware performance counters." OFF)
OPTION(ENABLE_EFENCE "Build with ElectricFence memory fault detection." OFF)
#OPTION(ENABLE_OVERKIT "Use Overkit" ON)
#OPTION(BUILD_OVERKIT "Build Overkit with PlasCom2 (if Overkit is enabled)" ON)
#OPTION(ENABLE_INLINE_ASSEMBLY "Get and build Fortran Overkit if enabled." OFF)
#OPTION(ENABLE_PETSC "Use PETSc" OFF)
#OPTION(ENABLE_STELLA "Use stella" OFF)
OPTION(BUILD_TESTS "Build tests" ON) 
OPTION(SHOWTESTS "List individual tests at config time" OFF)
#OPTION(INTEGRATED_TESTING "Enable longer integrated tests in testing suite" OFF)
OPTION(DISABLE_WARNINGS "Disable all compiler warnings" OFF)
OPTION(AGGRESSIVE_WARNINGS "Enable extra Fortran warnings (GCC only, sets DISABLE_WARNINGS to OFF)" OFF)
#OPTION(ENABLE_MAXWELL "Enable Maxwell CEM solver"  OFF)
#OPTION(ENABLE_CANTERA "Enable Cantera support for chemistry and equations of state" OFF)
#OPTION(USE_AMPI "Let build system know you plan on building with AMPI (optional)" OFF)
OPTION(USE_HYDRA "Enable Hydra offloading support" OFF)
OPTION(ENABLE_HYDRA "Enable Hydra offloading support" OFF)
#OPTION(FORCE_HDF5 "Proceed with configuration despite not finding paralell HDF5" OFF)
OPTION(DISABLE_FAILING_TESTS "Disable running of tests that are known to fail." OFF)
OPTION(ENABLE_THREAD_SANITIZER "Enable thread sanitizer thread debugger." OFF)
#OPTION(ENABLE_CVODE "Enable CVODE support" ON)
#OPTION(BUILD_CVODE "Force PlasCom2 to build its own CVODE" ON)
#OPTION(ENABLE_PLATO "Enable PLATO support" OFF)
OPTION(ENABLE_CUDA "Enable CUDA kernels and tests." OFF)
#OPTION(ENABLE_EOS_TIMERS "Enable high resolution EOS timings" OFF)
#OPTION(ENABLE_PROMETHEUS "Enable building new prometheus mechanisms" OFF)
#OPTION(ENABLE_Y6_TESTS "Enable Y6 testing" ON)

if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "${Green}${PROJECT_NAME}: No build type selected, default to Release${ColorReset}")
  set(CMAKE_BUILD_TYPE "Release")
endif()

SET(CMAKE_CXX_FLAGS_RELEASE "-g1 -O3 -DNDEBUG")
SET(CMAKE_Fortran_FLAGS_RELEASE "-g1 -O2 -DNDEBUG")
SET(CMAKE_C_FLAGS_RELEASE "-g1 -O3 -DNDEBUG")


#IF(ENABLE_EOS_TIMERS)
#  ADD_DEFINITIONS(-DEOS_TIMERS)
#ENDIF()


IF(ENABLE_CUDA)
  ENABLE_LANGUAGE(CUDA)
  INCLUDE(CheckLanguage)
  CHECK_LANGUAGE(CUDA)
  ADD_DEFINITIONS(-DHAVE_CUDA)
  MESSAGE(STATUS "${Green}${PROJECT_NAME}: CUDA           = ON${ColorReset}")
  FIND_PACKAGE(CUDA)
  INCLUDE_DIRECTORIES(${CUDA_TOOLKIT_ROOT_DIR}/include)
ENDIF()
IF(ENABLE_OPENMP)
  MESSAGE(STATUS "${Green}${PROJECT_NAME}: OpenMP         = ON${ColorReset}")
ENDIF()
#IF(NOT ENABLE_HDF5)
#  MESSAGE(STATUS "${Green}${PROJECT_NAME}: HDF            = OFF${ColorReset}")
#ENDIF()
#IF(ENABLE_OVERKIT)
#  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Overkit        = ON${ColorReset}")
#  IF(ENABLE_INLINE_ASSEMBLY)
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: Fortran Overkit = ON${ColorReset}")
#  ENDIF()
#ENDIF()
#IF(NOT ENABLE_CVODE)
#  SET(BUILD_CVODE FALSE)
#ENDIF()
#IF(BUILD_CVODE)
#  SET(ENABLE_CVODE TRUE)
#  MESSAGE(STATUS "${Green}${PROJECT_NAME}: CVODE          = ON/BUILD${ColorReset}") 
#ELSEIF(ENABLE_CVODE)
#  MESSAGE(STATUS "${Green}${PROJECT_NAME}: CVODE          = ON/EXTERNAL${ColorReset}")
#ENDIF()
#IF(ENABLE_CANTERA)
#  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Cantera        = ON${ColorReset}")
#ENDIF()
#IF(ENABLE_PETSC)
#  MESSAGE(STATUS "${Green}${PROJECT_NAME}: PETSC          = ON${ColorReset}")
#ENDIF()
#IF(ENABLE_STELLA)
#  MESSAGE(STATUS "${Green}${PROJECT_NAME}: STELLA         = ON${ColorReset}")
#ENDIF()
#IF(USE_AMPI)
#  MESSAGE(STATUS "${Green}${PROJECT_NAME}: AMPI           = ON${ColorReset}")
#ENDIF()
IF(BUILD_STATIC)
  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Static Build   = ON${ColorReset}")
ENDIF()
#IF(NOT BUILD_PROVENANCE)
#  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Provenance     = OFF${ColorReset}")
#ENDIF()
IF(ENABLE_COVERAGE)
  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Coverage Tests = ON${ColorReset}")
ENDIF()
#IF(ENABLE_Y6_TESTS)
#  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Y6 Testing = ON${ColorReset}")
#ELSE()
#  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Y6 Testing = OFF${ColorReset}")
#  ADD_DEFINITIONS(-DDISABLE_Y6)
#ENDIF()

#OPTION(ENABLE_AGGRESSIVE "Aggressive opt" OFF)

# Dynamic kernels on mac os x automatically
IF("${CMAKE_SYSTEM_VERSION}" STREQUAL "Darwin")
  SET(STATIC_KERNELS OFF)
ENDIF()
IF(APPLE)
  SET(STATIC_KERNELS OFF)
ENDIF()
IF(BUILD_STATIC)
  SET(STATIC_KERNELS ON)  
#  SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") 
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "IBM")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXIT_LINKER_FLAGS} -qstatic")
  ENDIF()
  SET(BUILD_SHARED_LIBS OFF)
ELSE()
  SET(BUILD_SHARED_LIBS ON)
ENDIF()

IF(USE_HYDRA) 
  SET(ENABLE_HYDRA TRUE)
ENDIF()

if(ENABLE_HYDRA)
  cmake_host_system_information(RESULT myhost QUERY HOSTNAME)
  if(${myhost} MATCHES "^(lassen|p9)")
    message(STATUS "${Green}${PROJECT_NAME}: Building Hydra for Lassen/P9${ColorReset}")
    add_definitions(-DHYDRA_LASSEN)
  endif()
  set(BUILD_STATIC ON)
  set(STATIC_KERNELS ON)
  set(ENABLE_OPENMP ON CACHE BOOL "" FORCE)
  set(BUILD_PROVENANCE OFF)
  add_definitions(-DUSE_HYDRA)
  find_package(CUDA)
  include_directories(${CUDA_TOOLKIT_ROOT_DIR}/include)
endif()

#
# guard against in-source builds
#
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "${Red}${PROJECT_NAME}: In-source builds not allowed. Please make a new directory and run CMake from there. You may need to remove CMakeCache.txt. ${ColorReset}")
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpic")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpic")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpic -w")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic -w")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpic -fstack-arrays")
  IF(AGGRESSIVE_WARNINGS)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra -pedantic -fbounds-check -fcheck-array-temporaries -fimplicit-none -ffpe-trap=invalid,zero,overflow -fbacktrace -std=f95")
  ENDIF()
  # using GCC
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpic")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic -std=c++11")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpic")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "IBM")
  # using IBM
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "XL")
  cmake_minimum_required(VERSION 3.11)
  # using IBM
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qstrict -qfloat=nomaf -qsuppress=1500-030")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qstrict -qfloat=nomaf -qsuppress=1500-030")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qstrict -qfloat=nomaf -qsuppress=1501-510")
  if (NOT BUILD_STATIC)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qPIC")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qPIC")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qPIC")
  endif()
endif()
IF(ENABLE_THREAD_SANITIZER)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=thread")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fsanitize=thread")
ENDIF()  


IF(ENABLE_EFENCE)
  FIND_LIBRARY(LIBEFENCE efence)
  IF(NOT LIBEFENCE)
    MESSAGE(WARNING "${Yellow}${PROJECT_NAME}: ElectricFence library not found. ${ColorReset}")
  ELSE()
    MESSAGE(STATUS "${Green}${PROJECT_NAME}: Building with ElectricFence memory fault detection ON. ${ColorReset}")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g")
  ENDIF()
ENDIF()

if(ENABLE_COVERAGE)
  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Coverage testing enabled. ${ColorReset}")
  SET(COVERAGE TRUE)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g -fprofile-arcs -ftest-coverage")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g -fprofile-arcs -ftest-coverage")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0 -g -fprofile-arcs -ftest-coverage")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgcov --coverage")
  SET(CMAKE_LIB_LINKER_FLAGS "${CMAKE_LIB_LINKER_FLAGS} -lgcov --coverage")
  SET(CMAKE_BUILD_TYPE "Debug")
endif()

#if(ENABLE_AGGRESSIVE)
#  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -xHost")
#  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -xHost")
#  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -xHost")
#endif()

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
get_filename_component (CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)

MESSAGE(STATUS "${Green}${PROJECT_NAME}: CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}${ColorReset}")
MESSAGE(STATUS "${Green}${PROJECT_NAME}: CXX_COMPILER_NAME = ${CXX_COMPILER_NAME}${ColorReset}")
MESSAGE(STATUS "${Green}${PROJECT_NAME}: CMAKE_CXX_COMPILER = ${CMAKE_CXX_COMPILER}${ColorReset}")
MESSAGE(STATUS "${Green}${PROJECT_NAME}: Fortran_COMPILER_NAME = ${Fortran_COMPILER_NAME}${ColorReset}")
MESSAGE(STATUS "${Green}${PROJECT_NAME}: CMAKE_Fortran_COMPILER = ${Fortran_COMPILER_ID}${ColorReset}")
MESSAGE(STATUS "${Green}${PROJECT_NAME}: CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}${ColorReset}")
MESSAGE(STATUS "${Green}${PROJECT_NAME}: CMAKE_Fortran_COMPILER_VERSION = ${CMAKE_Fortran_COMPILER_VERSION}${ColorReset}")
#SET( buildinfo_in          ${CMAKE_SOURCE_DIR}/include/buildinfo.in )
#SET( buildinfo             ${CMAKE_BINARY_DIR}/include/buildinfo.h  )
#CONFIGURE_FILE( ${buildinfo_in} ${buildinfo} @ONLY )

IF(ENABLE_OPENMP) 
  find_package(OpenMP)
  if (OPENMP_FOUND)
    MESSAGE(STATUS "${Green}${PROJECT_NAME}: Building with OpenMP enabled. ${ColorReset}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    ADD_DEFINITIONS(-DUSE_OMP)
  else()
    MESSAGE(WARNING "${Yellow}${PROJECT_NAME}: OpenMP not found. ${ColorReset}")
  endif()
ENDIF(ENABLE_OPENMP)

IF(ENABLE_KNL)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xMIC-AVX512")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xMIC-AVX512")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xMIC-AVX512")
ENDIF(ENABLE_KNL)


if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "XL")
  if(ENABLE_HYDRA)
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -qoffload")
    SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -qoffload")
  endif()
endif()

IF(DISABLE_WARNINGS)
  IF(NOT AGGRESSIVE_WARNINGS)
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -w")
    SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -w")
    SET(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -w")
  ENDIF()
ENDIF(DISABLE_WARNINGS)

# Compiled objects output paths
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single directory for all executables.")
set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE PATH "Single directory for all libraries and archives.")
mark_as_advanced (LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)
# Set up RPATH
if(NOT BUILD_STATIC)
  SET(CMAKE_MACOSX_RPATH TRUE)
  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()
IF(APPLE)
  set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
ENDIF()

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-ffree-form -pipe -ffree-line-length-0 -fconvert=big-endian -funroll-all-loops -fno-f2c -cpp -O3 -g")
  set (CMAKE_Fortran_FLAGS_DEBUG "-fno-f2c -cpp -O0 -g")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-fpic -f77rtl -O3 -g")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fpic -f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32 -g")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "XL")
  if (ENABLE_HYDRA)
    # Use -qoptimize=2, not -O2, to avoid spurious warning by the XL wrapper script
    set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -g1 -qoptimize=2 -qoffload -qsmp=omp -qsuppress=1501-510 -Xptxas --disable-optimizer-constants")
    # Need -lxlsmp before the other flags due to https://gitlab.kitware.com/cmake/cmake/issues/18518
    set (CMAKE_CXX_FLAGS_RELEASE "-lxlsmp ${CMAKE_CXX_FLAGS} -g1 -qoptimize=2 -qoffload -qsmp=omp -qsuppress=1501-510 -Xptxas --disable-optimizer-constants")
  else()
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -g1")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g1")
#    set (CMAKE_C_FLAGS_RELEASE "-O2 -g")
  endif()
else ()
  #message ("${PROJECT_NAME}: CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  #message ("${PROJECT_NAME}: Fortran compiler: " ${Fortran_COMPILER_NAME})
  message (STATUS "${Green}${PROJECT_NAME}: No optimized Fortran compiler flags are known, we just try -O3 ... ${ColorReset}")
  set (CMAKE_Fortran_FLAGS_RELEASE "-fpic -O3 -g1")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g1")
endif()

# 
# ========= USE IN EXISTING PROJECT ============
#
# In order to use this XPACC Project Template
# in an existing project COPY this file from this point
# on and paste it into your project's CMakeLists.txt 
# file *before* your existing INSTALL directives.
#
INCLUDE(CTest)
IF(BUILD_TESTING)
  SET(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
  MARK_AS_ADVANCED(BUILDNAME)
ENDIF(BUILD_TESTING)

#IF(AUTODETECT_MPI)
#  find_package(MPI REQUIRED)
#  include_directories(${MPI_INCLUDE_PATH})
#ENDIF()

# Enable the CS-Tools
#if(CS-Tools)
#  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/CS-Tools/CMakeLists.txt")
#    message(STATUS "${Green}${PROJECT_NAME}: Enabling CS Tools${ColorReset}")
#    add_subdirectory(CS-Tools)
#  else()
#    message(WARNING "${Yellow}${PROJECT_NAME}: CS-Tools were enabled but no files were found in CS-Tools. ${ColorReset}" 
#   "Try running the following commands: 
#   git submodule init
#    git submodule update"
#    "\nDisabling CS-Tools")
#    set(CS-Tools OFF)
#  endif()
#endif()


#IF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/googletest/")
#  add_subdirectory(googletest)
#ELSE()
#  add_definitions(-DNOGOOGLETEST)
#ENDIF()

#if(ENABLE_HDF5)
#  if(BUILD_STATIC OR APPLE)
#    set(HDF5_USE_STATIC_LIBRARIES TRUE)
#  endif()

#  SET(HDF5_PREFER_PARALLEL TRUE)
#  #SET(HDF5_FIND_DEBUG TRUE)
#  FIND_PACKAGE(HDF5 COMPONENTS C)
#  IF(HDF5_FOUND)  
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: HDF5_IS_PARALLEL       = ${HDF5_IS_PARALLEL}${ColorReset}")
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: HDF5_C_COMPILER_EXECUTABLE       = ${HDF5_C_COMPILER_EXECUTABLE}${ColorReset}")
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: HDF5_LIBRARY_DIRS      = ${HDF5_LIBRARY_DIRS}${ColorReset}")
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: HDF5_INCLUDE_DIRS      = ${HDF5_INCLUDE_DIRS}${ColorReset}")
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: HDF5_INCLUDE_DIR       = ${HDF5_INCLUDE_DIR}${ColorReset}")
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: HDF5_DEFINITIONS       = ${HDF5_DEFINITIONS}${ColorReset}")
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: HDF5_LIBRARIES         = ${HDF5_LIBRARIES}${ColorReset}")                              
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: HDF5_C_LIBRARIES       = ${HDF5_C_LIBRARIES}${ColorReset}")                         
#    #MESSAGE(STATUS "${Green}${PROJECT_NAME}: HDF5_CXX_LIBRARIES     = ${HDF5_CXX_LIBRARIES}${ColorReset}")                            
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: HDF5_Fortran_LIBRARIES = ${HDF5_Fortran_LIBRARIES}${ColorReset}")   
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: HDF5_C_COMPILER_EXECUTABLE = ${HDF5_C_COMPILER_EXECUTABLE}${ColorReset}")   
#    #MESSAGE(STATUS "${Green}${PROJECT_NAME}: HDF5_CXX_COMPILER_EXECUTABLE = ${HDF5_CXX_COMPILER_EXECUTABLE}${ColorReset}")   
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: HDF5_Fortran_COMPILER_EXECUTABLE = ${HDF5_Fortran_COMPILER_EXECUTABLE}${ColorReset}")
#    IF(HDF5_IS_PARALLEL) 
#      MESSAGE(STATUS "${Green}${PROJECT_NAME}: Enabling HDF5. ${ColorReset}")
#      SET(USE_HDF TRUE)
#    ELSE()
#      MESSAGE(WARNING "${Yellow}${PROJECT_NAME}: Disabling HDF5 (not parallel). ${ColorReset}")
#      SET(USE_HDF FALSE)
#    ENDIF(HDF5_IS_PARALLEL)
#  ENDIF() 
#  IF(USE_HDF OR FORCE_HDF5)
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: Setting up build for HDF5. ${ColorReset}")
#    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})                                                    
#    ADD_DEFINITIONS(-DENABLE_HDF5)
#    IF(FORCE_HDF5)
#      MESSAGE(WARNING "${Yellow}${PROJECT_NAME}: FindHDF5 failed, but forcing the build to continue regardless. ${ColorReset}")
#      MESSAGE(WARNING "${Yellow}${PROJECT_NAME}: This can be the right thing to do on systems that wrap hdf5 with their own compilers${ColorReset}")
#    ENDIF()
#  ELSE()
#    MESSAGE(FATAL_ERROR "${Red}${PROJECT_NAME}: Missing parallel HDF5. ${ColorReset}")
#  ENDIF()   
#ENDIF()

#
# Add suport for Cantera for combustion
# 2 options:
#   1. Put the path to Cantera in CMAKE_PREFIX_PATH
#   2. Define CANTERA_DIR in the environment
#
#IF(ENABLE_CANTERA)
#  MESSAGE(STATUS "${PROJECT_NAME}: Looking for Cantera")
#  FIND_PACKAGE(Cantera)
#  if(CANTERA_FOUND)
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: Found Cantera Version ${CANTERA_VERSION}${ColorReset}")
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: CANTERA_INCLUDE_DIR ${CANTERA_INCLUDE_DIRS${ColorReset}}")
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: CANTERA_CORE_LIBRARY ${CANTERA_CORE_LIBRARY}${ColorReset}")
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: CANTERA_EXTRA_LIBRARIES ${CANTERA_EXTRA_LIBRARIES}${ColorReset}")
    #message(STATUS "${Green}${PROJECT_NAME}: CANTERA_Fortran_LIBRARY ${CANTERA_Fortran_LIBRARY}${ColorReset}")
#
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: Setting up build to use Cantera. ${ColorReset}")
#    INCLUDE_DIRECTORIES(${CANTERA_INCLUDE_DIRS})
#    SET(CANTERA_LIBRARIES ${plascomExtraLibs}
#        #${CANTERA_Fortran_LIBRARY}
#        ${CANTERA_CORE_LIBRARY} 
#        ${CANTERA_EXTRA_LIBRARIES}
#        )
#    # need c++-11 for version >= 2.3
#    #if(${CANTERA_VERSION} VERSION_GREATER_EQUAL 2.3)
#    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#    #endif()
#    SET(USE_CANTERA "TRUE")
#    ADD_DEFINITIONS(-DHAVE_CANTERA)
#    # Hydra precludes the use of cantera
#    # can still be used in non-offloaded code though
#    # in this case, use HAVE_CANTERA (i.e. tests and initializations) inside the code
#    IF(NOT ENABLE_HYDRA)
#      ADD_DEFINITIONS(-DUSE_CANTERA)
#    ENDIF()
#  ELSE(CANTERA_FOUND)
#    MESSAGE(WARNING "${Yellow}${PROJECT_NAME}: Cantera not found, disabling Cantera. ${ColorReset}")
#    SET(USE_CANTERA "FALSE")
#  ENDIF(CANTERA_FOUND)
#
#ENDIF(ENABLE_CANTERA)
#
#IF(ENABLE_PROMETHEUS)
#  IF(USE_CANTERA)
#    SET(USE_PROMETHEUS "TRUE")
#    add_subdirectory(Prometheus)
#    ADD_DEFINITIONS(-DHAVE_PROMETHEUS)
#    ADD_DEFINITIONS(-DUSE_PROMETHEUS)
#  ELSE()
#    MESSAGE(FATAL_ERROR "${Red}${PROJECT_NAME}: Enabling Prometheus requires Cantera support. ${ColorReset}")
#  ENDIF()
#ENDIF()
#
#
#IF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/IX/")
#  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Building IX from source. ${ColorReset}")
#  ADD_SUBDIRECTORY(IX)
#  SET(IX_INCPATH ${CMAKE_SOURCE_DIR}/IX/include)
#  SET(RUNTEST ${EXECUTABLE_OUTPUT_PATH}/runtest)
#  SET(TESTRESULT ${EXECUTABLE_OUTPUT_PATH}/testresults)
#  SET(IX_LIB IX)
#ELSE()
#  FIND_LIBRARY(IX_LIB IX)
#  FIND_FILE(IX_HDR COMM.H)
#  FIND_PROGRAM(RUNTEST runtest)
#  FIND_PROGRAM(TESTRESULT testresults)
#  GET_FILENAME_COMPONENT(IX_LIBPATH ${IX_LIB} PATH)
#  GET_FILENAME_COMPONENT(IX_INCPATH ${IX_HDR} PATH)
#  GET_FILENAME_COMPONENT(IX_EXEPATH ${RUNTEST} PATH)
#  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Using existing IX install (${IX_EXEPATH}). ${ColorReset}")
#ENDIF()

#IF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/PCPP/")
#  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Building PCPP from source. ${ColorReset}")
#  ADD_SUBDIRECTORY(PCPP)
#  SET(PCPP_INCPATH ${CMAKE_SOURCE_DIR}/PCPP/include)
#  SET(PCPP_LIB PCPP)
#ELSE()
#  FIND_LIBRARY(PCPP_LIB PCPP)
#  FIND_FILE(PCPP_HDR PCPPTypes.H)
#  GET_FILENAME_COMPONENT(PCPP_LIBPATH ${PCPP_LIB} PATH)
#  GET_FILENAME_COMPONENT(PCPP_INCPATH ${PCPP_HDR} PATH)
#  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Using existing PCPP install (${PCPP_LIBPATH}). ${ColorReset}")
#ENDIF()

#IF(ENABLE_OVERKIT)
#  IF(BUILD_OVERKIT)
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: Configuring Overkit ... ${ColorReset}")
#    SET(OVERKIT_BUILD_TYPE ${CMAKE_BUILD_TYPE})
#    IF(OVERKIT_BUILD_TYPE STREQUAL "Debug")
#      SET(OVERKIT_BUILD_TYPE "FastDebug")
#    ENDIF()
#    IF(BUILD_STATIC)
#      SET(OVERKIT_SHARED OFF)
#    ELSE()
#      SET(OVERKIT_SHARED ON)
#    ENDIF()
#    SET(OVERKIT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Overkit")
#    SET(OVERKIT_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/Overkit")
#    INCLUDE(${OVERKIT_SOURCE_DIR}/OverkitSubproject.cmake)
#    CONFIGURE_OVERKIT(
#      SOURCE_DIR "${OVERKIT_SOURCE_DIR}"
#      BINARY_DIR "${OVERKIT_BINARY_DIR}"
#      INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin"
#      INSTALL_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib"
#      INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include"
#      OPTIONS -DCMAKE_BUILD_TYPE=${OVERKIT_BUILD_TYPE} -DXPACC=ON -DBUILD_SHARED_LIBS=${OVERKIT_SHARED}
#        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
#      QUIET
#    )
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: Configured Overkit ${Overkit_VERSION}${ColorReset}")
#  ELSE()
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: Looking for Overkit${ColorReset}")
#    FIND_PACKAGE(Overkit QUIET CONFIG REQUIRED)
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: Found Overkit ${Overkit_VERSION} in ${Overkit_DIR}${ColorReset}")
#    IF(NOT Overkit_OPTION_XPACC)
#      MESSAGE(FATAL_ERROR "${Red}Detected Overkit installation was not built with XPACC functionality enabled. ${ColorReset}")
#    ENDIF()
#  ENDIF()
#  ADD_DEFINITIONS(-DUSE_OVERKIT)
#  INCLUDE_DIRECTORIES("${Overkit_INCLUDES}")
#ENDIF()

# Fortran Overkit
#MESSAGE(STATUS "${Green}${PROJECT_NAME}: Looking for Fortran Overkit${ColorReset}")
#FIND_PACKAGE(OverkitFortran 2 QUIET CONFIG)
#IF(OverkitFortran_FOUND)
#  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Found Fortran Overkit ${OverkitFortran_VERSION} in ${OverkitFortran_DIR}${ColorReset}")
#  ADD_DEFINITIONS(-DUSE_OVERKIT_FORTRAN)
#  SET(OverkitFortran_LIBRARIES overkit-fortran)
#  SET(USE_OVERKIT_FORTRAN TRUE)
#ELSE()
#  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Fortran Overkit not found${ColorReset}")
#  SET(USE_OVERKIT_FORTRAN FALSE)
#  IF(ENABLE_INLINE_ASSEMBLY) # Then try to get and build Fortran Overkit
#    IF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/overkit-fortran/")
#      ADD_SUBDIRECTORY(overkit-fortran)
#      ADD_DEFINITIONS(-DUSE_OVERKIT_FORTRAN)
#      SET(Overkit_LIBRARIES overkit)
#      SET(USE_OVERKIT_FORTRAN TRUE)
#      INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}/overkit-fortran/src/core/mod")
#      INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}/overkit-fortran/src/support/mod")
#    ELSE()
#      MESSAGE(STATUS "${Green}${PROJECT_NAME}: Attempting to download Fortran Overkit${ColorReset}")
#      FIND_PACKAGE(Git QUIET)
#      IF(GIT_FOUND)
#        EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} clone https://github.com/majosm/overkit-fortran.git 
#          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#          RESULT_VARIABLE GOT_OVK_FORTRAN)
#        IF(NOT GOT_OVK_FORTRAN)
#          MESSAGE(FATAL_ERROR "${Red}${PROJECT_NAME}: Failed to get/clone Fortran Overkit (Reqd for inline assembly)${ColorReset}")
#        ELSE()
#          MESSAGE(STATUS "${Green}${PROJECT_NAME}: Download success for Fortran Overkit${ColorReset}")
#          ADD_SUBDIRECTORY(overkit-fortran)
#          ADD_DEFINITIONS(-DUSE_OVERKIT_FORTRAN)
#          SET(OverkitFortran_LIBRARIES overkit-fortran)
#          SET(USE_OVERKIT_FORTRAN TRUE)
#          INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}/overkit-fortran/src/core/mod")
#          INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}/overkit-fortran/src/support/mod")
#        ENDIF()
#      ELSE()
#        MESSAGE(FATAL_ERROR "${Red}${PROJECT_NAME}: Git (NOT FOUND) is required for inline assembly with Fortran Overkit${ColorReset}")
#      ENDIF()
#    ENDIF()
#  ENDIF()
#ENDIF()
#
#IF(ENABLE_PETSC)
#  FIND_PACKAGE(PETSc REQUIRED)
#  IF(NOT ${PETSC_FOUND})
#    MESSAGE(FATAL_ERROR "${Red}${PROJECT_NAME}: Unable to find PETSc. ${ColorReset}")
#  ENDIF()
#  ADD_DEFINITIONS(-DUSE_PETSC)
#  INCLUDE_DIRECTORIES(${PETSC_INCLUDES})
#ENDIF()
#
#IF(ENABLE_STELLA)
#  IF(NOT ENABLE_PETSC)
#    MESSAGE(FATAL_ERROR "${Red}${PROJECT_NAME}: Unable to configure with stella; needs PETSc. ${ColorReset}")
#  ENDIF()
#  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Building stella from source. ${ColorReset}")
#  add_subdirectory(stella)
#  SET(STELLA_LIB stella)
#  SET(STELLA_INCPATH ${CMAKE_SOURCE_DIR}/stella/src)
#  ADD_DEFINITIONS(-DUSE_STELLA)
#  INCLUDE_DIRECTORIES("${STELLA_INCPATH}")
#  IF(ENABLE_CEDAR)
#    ADD_DEFINITIONS(-DUSE_CEDAR)
#  ENDIF()
#ENDIF()
#
#IF(ENABLE_CVODE OR BUILD_CVODE)
#  IF(NOT BUILD_CVODE)
#    FIND_LIBRARY(CVODE_LIB sundials_cvode)
#    IF(NOT CVODE_LIB)
#      FIND_LIBRARY(CVODE_LIB cvode)
#    ENDIF(NOT CVODE_LIB)
#    IF(CVODE_LIB)
#      #  FIND_LIBRARY(NVEC_PAR_LIB sundials_nvecparallel)
#      #  FIND_LIBRARY(NVEC_SER_LIB sundials_nvecserial)
#      #  FIND_LIBRARY(SUNLINSOL_LIB sundials_sunlinsoldense)
#      #  FIND_LIBRARY(SUNMATRIX_LIB sundials_sunmatrixdense)
#      MESSAGE(STATUS "${Green}${PROJECT_NAME}: Found CVODE library: ${CVODE_LIB}${ColorReset}")
#      #  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Found SERIAL NVECTOR library: ${NVEC_SER_LIB}${ColorReset}")
#      #  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Found PARALLEL NVECTOR library: ${NVEC_PAR_LIB}${ColorReset}")
#      #  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Found SUNLINSOL library: ${SUNLINSOL_LIB}${ColorReset}")
#      #  MESSAGE(STATUS "${Green}${PROJECT_NAME}: Found SUNMATRIX library: ${SUNMATRIX_LIB}${ColorReset}")
#      GET_FILENAME_COMPONENT(CVODE_LIBPATH ${CVODE_LIB} PATH)
#      ADD_DEFINITIONS(-DHAVE_CVODE)
#      MESSAGE(STATUS "${Green}${PROJECT_NAME}: CVODE library path: ${CVODE_LIBPATH}${ColorReset}")
#      #  FIND_FILE(CVODE_HEADER cvode/cvode_direct.h)
#      FIND_FILE(CVODE_HEADER cvdense.h)
#      IF(NOT CVODE_HEADER)
#        MESSAGE(ERROR "${Red}${PROJECT_NAME}: CVODE header files not found. ${ColorReset}")
#      ENDIF(NOT CVODE_HEADER)
#      GET_FILENAME_COMPONENT(CVODE_INCPATH ${CVODE_HEADER} PATH)
#      #  MESSAGE(STATUS "${Green}${PROJECT_NAME}: CVODE include path: ${CVODE_INCPATH}/../ ${ColorReset}")
#      #  INCLUDE_DIRECTORIES(${CVODE_INCPATH}/../)
#      MESSAGE(STATUS "${Green}${PROJECT_NAME}: CVODE include path: ${CVODE_INCPATH}/ ${ColorReset}")
#      INCLUDE_DIRECTORIES(${CVODE_INCPATH}/)
#    ENDIF()
#  ENDIF()
#  IF(NOT CVODE_LIB)
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: Building internal CVODE library. ${ColorReset}")
#    ADD_SUBDIRECTORY(cvode)
#    ADD_DEFINITIONS(-DHAVE_CVODE)
#    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/cvode/include)
#    SET(CVODE_LIB PC2CVODE)
#  ENDIF(NOT CVODE_LIB)
#ENDIF()
#
#IF(ENABLE_PLATO)
#  MESSAGE(STATUS "Looking for Plato")
#  FIND_LIBRARY(PLATO_LIB plato)
#  IF(PLATO_LIB)
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: Found Plato library: ${PLATO_LIB}${ColorReset}")
#    GET_FILENAME_COMPONENT(PLATO_LIBPATH ${PLATO_LIB} PATH)
#    MESSAGE(STATUS "${Green}Plato library path: ${PLATO_LIBPATH}${ColorReset}")
#    FIND_FILE(PLATO_HEADER plato_library_interface_Cpp.h)
#    IF(NOT PLATO_HEADER)
#      MESSAGE(ERROR "${Red}${PROJECT_NAME}: PLATO header files not found.")
#    ENDIF(NOT PLATO_HEADER)
#    GET_FILENAME_COMPONENT(PLATO_INCPATH ${PLATO_HEADER} PATH)
#    MESSAGE(STATUS "${Green}${PROJECT_NAME}: Plato include path: ${PLATO_INCPATH}${ColorReset}")
#    INCLUDE_DIRECTORIES(${PLATO_INCPATH})
#    SET(USE_PLATO "TRUE")
#    ADD_DEFINITIONS(-DHAVE_PLATO)
#  ELSE(PLATO_LIB)
#    MESSAGE(WARNING "${Yellow}${PROJECT_NAME}: Plato not found. ${ColorReset}")
#    SET(USE_PLATO "FALSE")
#  ENDIF(PLATO_LIB)
#ENDIF(ENABLE_PLATO)
#
#
if(ENABLE_PAPI)
  FIND_PROGRAM(PAPI_AVAIL papi_avail)
  IF(NOT PAPI_AVAIL) 
    MESSAGE(WARNING "${Yellow}${PROJECT_NAME}: PAPI performance counter access not found. ${ColorReset}")
  ELSE()
    MESSAGE(STATUS "${Green}Building with PAPI performance counters enabled. ${ColorReset}")
    GET_FILENAME_COMPONENT(PAPI_PATH ${PAPI_AVAIL} PATH)
    SET(PAPI_PATH ${PAPI_PATH}/..)
    MESSAGE(STATUS "${Green}${PROJECT_NAME}: Found PAPI Install: ${PAPI_PATH}${ColorReset}")
    FIND_LIBRARY(PAPI_LIB papi PATHS ${PAPI_PATH}/lib NO_DEFAULT_PATH REQUIRED)
    MESSAGE(STATUS "${Green}${PROJECT_NAME}: Found PAPI Library: ${PAPI_LIB}${ColorReset}")
    GET_FILENAME_COMPONENT(PAPI_LIBPATH ${PAPI_LIB} PATH)
    ADD_DEFINITIONS(-DUSE_PAPIX)
    LINK_DIRECTORIES(${PAPI_LIBPATH})
    INCLUDE_DIRECTORIES(${PAPI_LIBPATH}/../include)
  ENDIF()
endif()
#
#
# Enable utilities
#add_subdirectory(utils)
#
#SET(PlasCom2_LibSources src/Simulation.C src/PlasCom2.C src/Stencil.C src/BoundaryCondition.C src/Info.C
#  src/NavierStokesBC.C src/Grid.C src/PC2Initialize.C src/PC2Util.C src/Halo.C
#  src/PC2Configuration.C src/PC2ReadDomains.C src/PC2WriteDomains.C src/State.C
#  src/EulerUtil.C src/ViscidUtil.C src/RunPlasCom2.C src/PlasComCM.C src/WENO.C
#  src/EOS.C src/NavierStokesDictionary.C src/Chemistry.C src/TransportModel.C
#  src/ChemistrySourceRHS.C src/CudaUtils.C src/NavierStokesRHS.C src/NavierStokesDomain.C
#  src/PAPIX.C src/Timers.C)
#
#
# Enable pre-built prometheus mechanisms
#set(PlasCom2_LibSources ${PlasCom2_LibSources} 
#  Prometheus/mechanisms/src/ETHYLENE.cpp Prometheus/mechanisms/src/ETHYLENE_EXTENDED.cpp
#  Prometheus/mechanisms/src/ETHYLENE_INERT.cpp Prometheus/mechanisms/src/ETHYLENE_INERT_EXT9.cpp
#  Prometheus/mechanisms/src/N2O2_reacN7.cpp Prometheus/mechanisms/src/ETHYLENE_PLASMA.cpp
#  Prometheus/mechanisms/src/SANDIEGO.cpp Prometheus/mechanisms/src/SANDIEGO_PLASMA.cpp
#  Prometheus/mechanisms/src/PrometheusInterface.cpp)

#if(ENABLE_HDF5)
#  set(PlasCom2_LibSources ${PlasCom2_LibSources} src/PC2IO.C)
#endif()
#
#IF(USE_OVERKIT_FORTRAN)
#  LIST(APPEND PlasCom2_LibSources src/OverkitFortranInterface.F90)
#ENDIF()

#SET(PlasCom2_KernelSources kernels/Limits.f90 kernels/barebones.f90 kernels/Operators.f90 kernels/MetricOps.f90
  #kernels/Roe.f90 kernels/RungeKutta.f90 kernels/Grid.f90 kernels/Euler.f90 kernels/EOS.f90
  #kernels/Viscid.f90 kernels/SATUtil.f90 kernels/NavierStokesRHSWENO.f90 
  #kernels/NavierStokesUtil.f90)

if(ENABLE_HYDRA)
  set(kernel_prefix hydra)
else()
  set(kernel_prefix kernels)
endif()

# this is an order that satisfies build order dependence
SET(PlasCom2_KernelSources kernels/Limits.f90 kernels/barebones.f90 kernels/NavierStokesUtil.f90
kernels/EOS.f90 ${kernel_prefix}/Operators.f90 kernels/Grid.f90 kernels/Euler.f90 kernels/Viscid.f90
${kernel_prefix}/MetricOps.f90 ${kernel_prefix}/Roe.f90 kernels/SATUtil.f90 ${kernel_prefix}/NavierStokesRHSWENO.f90)

#if(ENABLE_HYDRA)
#  # Add Tangram kernels
#  SET(PlasCom2_KernelSources ${PlasCom2_KernelSources} hydra/tangram.cpp hydra/tangram.cu)
#endif()

if(ENABLE_HYDRA)
  SET(PlasCom2_HydraSources src/Hydra.C)
endif()

#MESSAGE(STATUS "${Green}${PROJECT_NAME}: Library Sources: ${PlasCom2_LibSources}${ColorReset}")
MESSAGE(STATUS "${Green}${PROJECT_NAME}: Kernel Sources: ${PlasCom2_KernelSources}${ColorReset}")


#
# Preprocess the sources with ICE
#
#if(USE_ICE)
#
#  ADD_DEFINITIONS(-DHAVE_ICE)
#
#  # list of includes that are candidates for ICE-Transformation
#  FILE(GLOB PlasCom2_IncludeSources include/*.H)
#
#  list(APPEND tuningFileIn ${CMAKE_CURRENT_SOURCE_DIR}/ice/tuneHand.opt)
#  if(USE_MOYA)
#    list(APPEND tuningFileIn ${CMAKE_CURRENT_SOURCE_DIR}/ice/tuneHandMoya.yaml)
#  endif()
#
#  # comment this out for now ... timothy is working on it
#  if(HAVE_ICE_ROSE)
#    list(APPEND tuningFileIn ${CMAKE_CURRENT_SOURCE_DIR}/ice/tuneWenoRose.opt)
#  endif()

#  # we should be able to move this tuning file work
#  #   into a function provided by CSTUF
#  # get any user defined tuning files
#  if(TUNING_FILE)
#    foreach(file ${TUNING_FILE})
#      list(APPEND tuningFileIn ${file})
#    endforeach()
#  endif()
#
#  set(tuningFile ${CMAKE_CURRENT_BINARY_DIR}/ice/tuningFile.opt)
#  set(masterTuningFile ${CMAKE_CURRENT_BINARY_DIR}/ice/masterTuningFile.opt)
#
#  FILE(REMOVE ${masterTuningFile})
#  foreach(file ${tuningFileIn})
#    FILE(READ ${file} tmpFile)
#    FILE(APPEND ${masterTuningFile} "${tmpFile}")
#  endforeach()
#
#  configure_file(${masterTuningFile} ${tuningFile} @ONLY)
#
#  # LibSources
#  foreach(fileIn ${PlasCom2_LibSources})
#    PreProcessICE(${CMAKE_CURRENT_SOURCE_DIR}/${fileIn} ${tuningFile} 
#      ${CMAKE_CURRENT_BINARY_DIR}/ice/src fileOut "")
#    list(APPEND PlasCom2_LibSources-ICE ${fileOut})
#  endforeach()
#  set(PlasCom2_LibSources ${PlasCom2_LibSources-ICE})
#
#  # generate full path for kernel sources for dependency search
#  foreach(fileIn ${PlasCom2_KernelSources})
#    list(APPEND moduleCandidates "${CMAKE_CURRENT_SOURCE_DIR}/${fileIn}")
#  endforeach()
#
#  # KernelSources
#  foreach(fileIn ${moduleCandidates})
#    FortranFileDependence("${fileIn}" "${moduleCandidates}" dependents)
#
#    PreProcessICE(${fileIn} ${tuningFile} 
#      ${CMAKE_CURRENT_BINARY_DIR}/ice/kernels fileOut "${dependents}")
#    list(APPEND PlasCom2_KernelSources-ICE ${fileOut})
#  endforeach()
#  set(PlasCom2_KernelSources ${PlasCom2_KernelSources-ICE})
#
#  # IncludeSources
#  foreach(fileIn ${PlasCom2_IncludeSources})
#    PreProcessICE(${fileIn} ${tuningFile} 
#      ${CMAKE_CURRENT_BINARY_DIR}/ice/include fileOut "")
#    list(APPEND PlasCom2_IncludeSources-ICE ${fileOut})
#  endforeach()
#  add_custom_target(iceProcessIncludeFiles ALL DEPENDS ${PlasCom2_IncludeSources-ICE})
#endif()
#
#
#IF(ENABLE_HYDRA AND USE_ICE)
#  message("Using ICE for RKSUM offloading")
#  PreProcessICE(${CMAKE_CURRENT_SOURCE_DIR}/kernels/RungeKutta.f90 
#    ${CMAKE_CURRENT_SOURCE_DIR}/ice/addHydraRK.opt
#    ${CMAKE_CURRENT_BINARY_DIR}/ice/kernels fileOut  "")
#  list(APPEND PlasCom2_KernelSources ${fileOut})
#ELSE()
#  list(APPEND PlasCom2_KernelSources kernels/RungeKutta.f90)
#ENDIF()


# Project-specific constructs
#IF(BUILD_STATIC)
#  ADD_LIBRARY(PlasCom2 STATIC ${PlasCom2_LibSources})
#  set_target_properties(PlasCom2 PROPERTIES LINK_SEARCH_START_STATIC 1)
#  set_target_properties(PlasCom2 PROPERTIES LINK_SEARCH_END_STATIC 1)
#ELSE()
#  ADD_LIBRARY(PlasCom2 SHARED ${PlasCom2_LibSources})
#ENDIF()
#target_compile_features(PlasCom2 PUBLIC cxx_std_11) # Require C++ 11

# the ice part needs to come first so I can get any modified includes
#if(USE_ICE)
#  TARGET_INCLUDE_DIRECTORIES(PlasCom2 PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/ice/include)
#  ADD_DEPENDENCIES(PlasCom2 iceProcessIncludeFiles)
#ENDIF()

# Enable pre-built prometheus mechanisms
#TARGET_INCLUDE_DIRECTORIES(PlasCom2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Prometheus/mechanisms/include)
#
#TARGET_INCLUDE_DIRECTORIES(PlasCom2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
#ADD_DEPENDENCIES(PlasCom2 initializedomainheader initgridsolnheader gengridheader processgeometryheader)
#
#IF(BUILD_PROVENANCE)
#  ADD_DEPENDENCIES(PlasCom2 repoheader)
#ELSE()
#  FILE(WRITE ${CMAKE_BINARY_DIR}/include/repositoryinfo.h "")
#ENDIF()
#TARGET_LINK_LIBRARIES(PlasCom2 ${PCPP_LIB} ${IX_LIB} ${PAPI_LIB} ${Overkit_LIBRARIES}
#  ${OverkitFortran_LIBRARIES} ${PETSC_LIBRARIES} ${STELLA_LIB} ${CVODE_LIB} ${NVEC_SER_LIB} ${NVEC_PAR_LIB} ${SUNLINSOL_LIB} ${SUNMATRIX_LIB} ${PLATO_LIB} ${CUDA_LIBRARIES})
#
#
#IF(USE_HDF)
#  TARGET_LINK_LIBRARIES(PlasCom2 ${HDF5_LIBRARIES})
#ENDIF()

#IF(ENABLE_CANTERA)
#  TARGET_LINK_LIBRARIES(PlasCom2 ${CANTERA_LIBRARIES})
#ENDIF()

#INSTALL(TARGETS PlasCom2 LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

IF(STATIC_KERNELS) 
  if(ENABLE_HYDRA)
    CUDA_ADD_LIBRARY(OK STATIC ${PlasCom2_KernelSources})
  else()
    ADD_LIBRARY(OK STATIC ${PlasCom2_KernelSources})
  endif()
ELSE()
  ADD_LIBRARY(OK SHARED ${PlasCom2_KernelSources})
ENDIF()
SET_TARGET_PROPERTIES(OK PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include")
INSTALL(TARGETS OK LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

if(ENABLE_HYDRA)
  ADD_LIBRARY(Hydra STATIC ${PlasCom2_HydraSources})
  TARGET_INCLUDE_DIRECTORIES(Hydra
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/include
  )
  INSTALL(TARGETS Hydra LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
endif()

#IF(USE_AMPI) 
#  # Disable renaming std::exit for AMPI
#  set_source_files_properties(src/PC2Util.C PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DAMPI_RENAME_EXIT=0")
#  IF(ENV{AMPI_HOME})
#    SET(AMPI_HOME $ENV{AMPI_HOME})
#  ELSE()
#    GET_FILENAME_COMPONENT(AMPI_BIN ${CMAKE_CXX_COMPILER} PATH)
#    SET(AMPI_HOME "${AMPI_BIN}/..")
#  ENDIF()
#  #FILE(COPY ${AMPI_HOME}/include/mpif.h DESTINATION ${CMAKE_BINARY_DIR}/include)
#  IF(NOT BUILD_STATIC)
#	  #set(CMAKE_CXX_STANDARD_LIBRARIES "${CMAKE_CXX_STANDARD_LIBRARIES} -charm-shared")
#	  target_link_libraries(PlasCom2 INTERFACE -charm-shared)
#	  target_link_libraries(IX INTERFACE -charm-shared)
#  endif()
#ENDIF()
#
#IF(ENABLE_MPIP)
#  SET(MPI_LINK_FLAGS "${MPI_LINK_FLAGS} -lmpiP")
#ENDIF()

#set_source_files_properties(src/PC2Main.C
#    PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
#ADD_EXECUTABLE(plascom2x src/PC2Main.C)
#TARGET_LINK_LIBRARIES(plascom2x PlasCom2 Kernels ${PCPP_LIB} ${IX_LIB}
#  ${MPI_CXX_LIBRARIES} ${MPI_LIBRARIES} ${LIBEFENCE})
#SET_TARGET_PROPERTIES(plascom2x PROPERTIES LINK_FLAGS "${LINK_FLAGS} ${MPI_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")
#
#INSTALL(TARGETS plascom2x RUNTIME DESTINATION bin)
#
#ADD_EXECUTABLE(pc2fileinfo src/pc2fileinfo.C)
#TARGET_LINK_LIBRARIES(pc2fileinfo PlasCom2 Kernels ${PCPP_LIB} ${IX_LIB} ${MPI_LIBRARIES} ${LIBEFENCE})
#SET_TARGET_PROPERTIES(pc2fileinfo PROPERTIES LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
#INSTALL(TARGETS pc2fileinfo RUNTIME DESTINATION bin)
#
#ADD_EXECUTABLE(pc2compare src/pc2compare.C)
#TARGET_LINK_LIBRARIES(pc2compare PlasCom2 Kernels ${PCPP_LIB} ${IX_LIB} ${MPI_LIBRARIES} ${LIBEFENCE})
#SET_TARGET_PROPERTIES(pc2compare PROPERTIES LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
#INSTALL(TARGETS pc2compare RUNTIME DESTINATION bin)

#if(ENABLE_HYDRA)
#  # -lxlcufdevice is required for allocate() support in offloaded Fortran regions
#  TARGET_LINK_LIBRARIES(PlasCom2 -lxlcufdevice)
#  TARGET_LINK_LIBRARIES(plascom2x   Hydra)
#  TARGET_LINK_LIBRARIES(pc2fileinfo Hydra)
#  TARGET_LINK_LIBRARIES(pc2compare  Hydra)
#  # SET_TARGET_PROPERTIES(plascom2x   PROPERTIES LINK_FLAGS "${LINK_FLAGS} -qcuda")
#  # SET_TARGET_PROPERTIES(pc2fileinfo PROPERTIES LINK_FLAGS "${LINK_FLAGS} -qcuda")
#  # SET_TARGET_PROPERTIES(pc2compare  PROPERTIES LINK_FLAGS "${LINK_FLAGS} -qcuda")
#endif()

# Add target to build the documentation
IF(BUILD_DOCUMENTATION)
  add_subdirectory(doc) 
endif()

# ========= USE IN EXISTING PROJECT ==============

IF(BUILD_TESTS) 
  ADD_SUBDIRECTORY(Testing)
ENDIF()


# Print configuration
get_cmake_property(_variableNames VARIABLES)
list (SORT _variableNames)

cmake_policy(SET CMP0054 NEW)
foreach (v ${_variableNames})
    if("${v}" MATCHES "^ENABLE_" OR "${v}" MATCHES "^USE_")
        if("${${v}}" STREQUAL "on" OR "${${v}}" STREQUAL "ON" OR "${${v}}" STREQUAL "1" OR "${${v}}" STREQUAL "TRUE" OR "${${v}}" STREQUAL "true")
            string(REPLACE "ENABLE_" "E_" _name ${v})
            string(REPLACE "USE_" "U_" _name ${_name})
            set(opts_enabled "${opts_enabled}${_name} ")
        endif()
    endif()
endforeach()

message("${Green}")
message("==============================")
message("${PROJECT_NAME} configuration: ")
message("  C compiler:      ${CMAKE_C_COMPILER} [${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}]")
message("    C flags:         ${CMAKE_C_FLAGS}")
message("  C++ compiler:    ${CMAKE_CXX_COMPILER} [${CMAKE_CXX_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}]")
message("    C++ flags:       ${CMAKE_CXX_FLAGS}")
message("  F90 compiler:    ${CMAKE_Fortran_COMPILER} [${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}]")
message("    F90 flags:       ${CMAKE_Fortran_FLAGS}")
message("  Enabled options: ${opts_enabled}")
message("==============================" )
message("${ColorReset}")
