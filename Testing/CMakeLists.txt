
FILE(GLOB LIB_CPP_SRCS "src/Test*.C")
FILE(GLOB LIB_F90_SRCS "src/Test*.f90")
FILE(WRITE ${CMAKE_SOURCE_DIR}/utils/DisableSuites.txt "")
IF(NOT USE_ICE)
  FILE(GLOB ICE_TEST_FILES "src/*ICE*")
  LIST(REMOVE_ITEM LIB_CPP_SRCS ${ICE_TEST_FILES})
  FILE(APPEND ${CMAKE_SOURCE_DIR}/utils/DisableSuites.txt "ICE\n")
  MESSAGE(STATUS "${Yellow}Testing framework excluding disabled ICE tests in:${ColorReset} ${ICE_TEST_FILES}")
ENDIF()
SET(Testing_KernelSourcesCPP src/kernels/CSTestKernelsCXX.C)
SET(Testing_KernelSourcesFortran src/kernels/CSTestKernelsFortran.f90 src/kernels/CSTestKernelsFortranAltDesc.f90)

if(USE_ICE)
  # build a list of all the optimization files to be used
  set(tuningFileIn ${CMAKE_CURRENT_SOURCE_DIR}/ice/tuneTesting.opt)
  if(HAVE_ICE_ROSE)
    list(APPEND tuningFileIn ${CMAKE_CURRENT_SOURCE_DIR}/ice/tuneTestingRose.opt)
  endif()

  set(masterTuningFile ${CMAKE_CURRENT_BINARY_DIR}/ice/masterTuningFile.opt)
  file(WRITE ${masterTuningFile} "")
  foreach(tuningFile ${tuningFileIn})
    cat(${tuningFile} ${masterTuningFile})
  endforeach()
  #if(USE_MOYA)
    #set(tuningFileIn ${CMAKE_CURRENT_SOURCE_DIR}/ice/tuneHandMoya.yaml)
  #endif()
  set(tuningFile ${CMAKE_CURRENT_BINARY_DIR}/ice/tuningFile.opt)
  #configure_file(${tuningFileIn} ${tuningFile} @ONLY)
  configure_file(${masterTuningFile} ${tuningFile} @ONLY)

  # C++ LibSources
  foreach(fileIn ${LIB_CPP_SRCS})
    PreProcessICE(${fileIn} ${tuningFile} 
      ${CMAKE_CURRENT_BINARY_DIR}/ice/src fileOut "")
    list(APPEND LIB_CPP_SRCS-ICE ${fileOut})
  endforeach()
  set(LIB_CPP_SRCS ${LIB_CPP_SRCS-ICE})

  # KernelSources
  foreach(fileIn ${Testing_KernelSourcesFortran})
    FortranFileDependence("${CMAKE_CURRENT_SOURCE_DIR}/${fileIn}" "${moduleCandidates}" dependents)

    PreProcessICE(${CMAKE_CURRENT_SOURCE_DIR}/${fileIn} ${tuningFile} 
      ${CMAKE_CURRENT_BINARY_DIR}/ice/kernels fileOut "${dependents}")
    list(APPEND Testing_KernelSourcesFortran-ICE ${fileOut})
  endforeach()
  set(Testing_KernelSourcesFortran ${Testing_KernelSourcesFortran-ICE})

  # KernelSources
  foreach(fileIn ${Testing_KernelSourcesCPP})
    PreProcessICE(${CMAKE_CURRENT_SOURCE_DIR}/${fileIn} ${tuningFile} 
      ${CMAKE_CURRENT_BINARY_DIR}/ice/kernels fileOut "")
      # doing something here to pass in library dependencies
      #${CMAKE_CURRENT_BINARY_DIR}/ice/kernels fileOut "ice.C" "TestingKernels")
    list(APPEND Testing_KernelSourcesCPP-ICE ${fileOut})
  endforeach()
  set(Testing_KernelSourcesCPP ${Testing_KernelSourcesCPP-ICE})

  ## Fortran LibSources
  #foreach(fileIn ${LIB_F90_SRCS})
    #PreProcessICE(${fileIn} ${tuningFile} 
      #${CMAKE_CURRENT_BINARY_DIR}/ice/src fileOut )
    #list(APPEND LIB_F90_SRCS-ICE ${fileOut})
  #endforeach()
  #set(LIB_F90_SRCS ${LIB_F90_SRCS-ICE})
endif()

set (LIB_SRCS ${LIB_CPP_SRCS} ${LIB_F90_SRCS})
set (ALL_SRCS ${LIB_SRCS} src/OKTestDriver.C src/OKParallelTestDriver.C)

IF(NOT BUILD_STATIC)
  set_source_files_properties(${ALL_SRCS} PROPERTIES COMPILE_FLAGS "-fPIC" )
ENDIF()
MESSAGE(STATUS "ALL TEST SOURCES: ${ALL_SRCS}")
IF(USE_ICE)
  MESSAGE(STATUS "ALL ICE C++ KERNEL SOURCES: ${Testing_KernelSourcesCPP}")
  MESSAGE(STATUS "ALL ICE Fortran KERNEL SOURCES: ${Testing_KernelSourcesFortran}")
ENDIF()
# rpath settings
IF(NOT BUILD_STATIC)
  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
ENDIF()

#ADD_CUSTOM_TARGET(
# testingheader
# COMMAND python ${CMAKE_SOURCE_DIR}/utils/DetectTests.py WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/utils
# )

IF(BUILD_STATIC)
  ADD_LIBRARY(OKUnitTest STATIC ${LIB_SRCS})
  set_target_properties(OKUnitTest PROPERTIES LINK_SEARCH_START_STATIC 1)
  set_target_properties(OKUnitTest PROPERTIES LINK_SEARCH_END_STATIC 1)
ELSE()
  ADD_LIBRARY(OKUnitTest SHARED ${LIB_SRCS})
ENDIF()
TARGET_INCLUDE_DIRECTORIES(OKUnitTest 
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/utils/include
)

#ADD_DEPENDENCIES(PC2UnitTest testingheader)

TARGET_LINK_LIBRARIES(OKUnitTest OK )
IF(USE_ICE)
  IF(BUILD_STATIC)
    ADD_LIBRARY(TestingKernels STATIC ${Testing_KernelSourcesFortran} ${Testing_KernelSourcesCPP})
  ELSE()
    ADD_LIBRARY(TestingKernels SHARED ${Testing_KernelSourcesFortran} ${Testing_KernelSourcesCPP})
  ENDIF()
  SET_TARGET_PROPERTIES(TestingKernels PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include")
  TARGET_INCLUDE_DIRECTORIES(TestingKernels 
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/include
  )
  TARGET_LINK_LIBRARIES(PC2UnitTest TestingKernels)
ENDIF()

ADD_EXECUTABLE(plascom2_test src/PC2TestDriver.C)
TARGET_INCLUDE_DIRECTORIES(plascom2_test 
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
TARGET_LINK_LIBRARIES(ok_test OKUnitTest OK ${NVEC_LIB})
#TARGET_LINK_LIBRARIES(plascom2_test PC2UnitTest PlasCom2 ${PCPP_LIB} ${IX_LIB} Kernels)
if(ENABLE_HYDRA)
  TARGET_LINK_LIBRARIES(ok_test Hydra)
endif()

ADD_DEFINITIONS(-D_PC2_PARALLEL_)

ADD_EXECUTABLE(ok_parallel_test src/OKParallelTestDriver.C)
#TARGET_LINK_LIBRARIES(plascom2_parallel_test PC2UnitTest PlasCom2 ${PCPP_LIB} ${IX_LIB} Kernels)
TARGET_LINK_LIBRARIES(ok_parallel_test OKUnitTest OK ${NVEC_LIB})
IF(ENABLE_MPIP)
  SET_TARGET_PROPERTIES(ok_parallel_test PROPERTIES LINK_FLAGS "${LINK_FLAGS} -lmpiP")
ENDIF()
if(ENABLE_HYDRA)
  TARGET_LINK_LIBRARIES(ok_parallel_test Hydra)
endif()
TARGET_INCLUDE_DIRECTORIES(ok_parallel_test 
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)


#ADD_EXECUTABLE(pc2twenke src/pc2twenke.C ../src/Stencil.C)
#TARGET_LINK_LIBRARIES(pc2twenke PlasCom2 ${PCPP_LIB} ${IX_LIB} Kernels)
#if(ENABLE_HYDRA)
#  TARGET_LINK_LIBRARIES(pc2twenke Hydra)
#endif()
#TARGET_INCLUDE_DIRECTORIES(pc2twenke 
#  PRIVATE
#    ${CMAKE_CURRENT_SOURCE_DIR}/include
#    ${CMAKE_CURRENT_SOURCE_DIR}/../include
#)

#ADD_EXECUTABLE(pc2twenke2 src/pc2twenke2.C ../src/Stencil.C)
#TARGET_LINK_LIBRARIES(pc2twenke2 PlasCom2 ${PCPP_LIB} ${IX_LIB} Kernels)
#if(ENABLE_HYDRA)
#  TARGET_LINK_LIBRARIES(pc2twenke2 Hydra)
#endif()
#TARGET_INCLUDE_DIRECTORIES(pc2twenke2
#  PRIVATE
#    ${CMAKE_CURRENT_SOURCE_DIR}/include
#    ${CMAKE_CURRENT_SOURCE_DIR}/../include
#)

#IF(MPI_LINK_FLAGS)
if(ENABLE_HYDRA)
 SET_TARGET_PROPERTIES(ok_parallel_test PROPERTIES LINK_FLAGS "-qcuda")
 SET_TARGET_PROPERTIES(ok_test PROPERTIES LINK_FLAGS "-qcuda")
# SET_TARGET_PROPERTIES(pc2twenke PROPERTIES LINK_FLAGS "-qcuda")
endif()
#ENDIF()

#
# Build the suites file that specifies which tests to run
# First read the defaults from the master list
#
SET(numSerialSuites 0)
message(STATUS "${Blue}Discovering serial unit suites ...${ColorReset}")
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/Scripts/serialUnitSuites.txt serialTestSuites)
string(REGEX REPLACE "\n" " " serialTestSuiteList ${serialTestSuites})
#set(testSuiteList ${testSuites})
separate_arguments(serialTestSuiteList)

#
# Add optional suites that are enabled by various options
#
IF(NOT DISABLE_FAILING_TESTS)
  list(APPEND serialTestSuiteList FAILS)
ENDIF()

IF(ENABLE_CUDA)
  list(APPEND serialTestSuiteList CUDA)
ENDIF()


# ICE tests
if(USE_ICE)
  list(APPEND serialTestSuiteList ICE)
  # tests specific to modules enabled by ICE
  # Rose
  if(HAVE_ICE_ROSE)
    list(APPEND serialTestSuiteList ICE-Rose)
  endif()
  # Pips
  #if(ICE_HAVE_PIPS)
    #list(APPEND serialTestSuiteList ICE-Pips)
  #endif()
endif()

#
# write the list out to the file that is read by the test framework
#
SET(serialSuiteFile ${CMAKE_CURRENT_BINARY_DIR}/serialUnitSuites)
FILE(WRITE ${serialSuiteFile} "")
FOREACH(suiteName ${serialTestSuiteList})
  message(STATUS "${Cyan}\t Running tests from suite ${suiteName}${ColorReset}")
  file(APPEND ${serialSuiteFile} "${suiteName} \n")
  MATH(EXPR numSerialSuites "${numSerialSuites} + 1") 
ENDFOREACH()

#
# parallel test suites
#
SET(numParallelSuites 0)
message(STATUS "${Blue}Discovering parallel unit testing suites ...${ColorReset}")
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/Scripts/parallelUnitSuites.txt parallelTestSuites)
string(REGEX REPLACE "\n" " " parallelTestSuiteList ${parallelTestSuites})
#set(testSuiteList ${testSuites})
separate_arguments(parallelTestSuiteList)

#
# Add optional suites that are enabled by various options
#
IF(NOT DISABLE_FAILING_TESTS)
  list(APPEND parallelTestSuiteList FAILS)
ENDIF()

IF(ENABLE_PAPI)
  LIST(APPEND parallelTestSuiteList PAPIX)
ENDIF()

#
# write the list out to the file that is read by the test framework
#
SET(parallelSuiteFile ${CMAKE_CURRENT_BINARY_DIR}/parallelUnitSuites)
FILE(WRITE ${parallelSuiteFile} "")
FOREACH(suiteName ${parallelTestSuiteList})
  MESSAGE(STATUS "${Cyan}\t Running tests from suite ${suiteName}${ColorReset}")
  FILE(APPEND ${parallelSuiteFile} "${suiteName} \n")
  MATH(EXPR numParallelSuites "${numParallelSuites} + 1")
ENDFOREACH()


#ADD_TEST(TestSetup ${RUNTEST} -f ${CMAKE_CURRENT_SOURCE_DIR}/Scripts/test_setup.sh -s ${CMAKE_CURRENT_SOURCE_DIR} -b ${CMAKE_BINARY_DIR}/bin -o setup.txt)

#ADD_TEST(RunTestScripts ${RUNTEST} -l ${CMAKE_CURRENT_SOURCE_DIR}/Scripts/tests.list -s ${CMAKE_CURRENT_SOURCE_DIR} -b ${CMAKE_BINARY_DIR}/bin -v 2 -o plascom2_testresults.txt)
#SET_PROPERTY(TEST RunTestScripts APPEND PROPERTY DEPENDS TestSetup)
#SET_PROPERTY(TEST PlasCom2:RunParallelPlatformTests APPEND PROPERTY TIMEOUT 7200)

#ADD_TEST(RunSpecificTestScript    ${RUNTEST} -f ${CMAKE_CURRENT_SOURCE_DIR}/Platforms/test_stub.sh -s ${CMAKE_CURRENT_SOURCE_DIR} -b ${CMAKE_BINARY_DIR}/bin -o plascom2_testresults.txt)
#SET_PROPERTY(TEST RunSpecificTestScript APPEND PROPERTY DEPENDS TestSetup)

#ADD_TEST(PlasCom2:RunSerialPlatformTests ${RUNTEST} -p ${CMAKE_CURRENT_SOURCE_DIR}/Platforms/platforms -s ${CMAKE_CURRENT_SOURCE_DIR} -b ${CMAKE_BINARY_DIR}/bin -v 3 -o ok_testresults.txt)
#SET_PROPERTY(TEST PlasCom2:RunSerialPlatformTests APPEND PROPERTY DEPENDS TestSetup)

# Run these one suite at a time
FOREACH(suiteName ${parallelTestSuiteList})
  #MESSAGE(STATUS "${Cyan}\t Running tests from suite ${suiteName}${ColorReset}")
  ADD_TEST(PlasCom2:RunParallelPlatformTests:${suiteName} ${RUNTEST} -p ${CMAKE_CURRENT_SOURCE_DIR}/Platforms/parallel_platforms -s ${CMAKE_CURRENT_SOURCE_DIR} -b ${CMAKE_BINARY_DIR}/bin -a ${suiteName} -v 3 -o plascom2_testresults.txt)
  SET_PROPERTY(TEST PlasCom2:RunParallelPlatformTests:${suiteName} APPEND PROPERTY DEPENDS TestSetup)
ENDFOREACH()

#ADD_TEST(PlasCom2:RunsInSerial ${TESTRESULT} PlasCom2:RunsInSerial plascom2_testresults.txt)
#ADD_TEST(PlasCom2:TestStubWorks ${TESTRESULT} PlascTestStubWorks plascom2_testresults.txt)
#ADD_TEST(PlasCom2:RunsInParallel ${TESTRESULT} PlasCom2:RunsInParallel plascom2_testresults.txt)


#ADD_TEST(Compare:Legacy:Same ${EXECUTABLE_OUTPUT_PATH}/pc2compare PlasComCM/hdf5/legacyCMFile.h5 PlasComCM/hdf5/legacyCMFile.h5)
#ADD_TEST(Compare:Legacy:Different ${EXECUTABLE_OUTPUT_PATH}/pc2compare PlasComCM/hdf5/flocm2d_005.h5 InviscidTests/AcousticPulsePeriodic2D/CM/flocm_acousticpulse_1000.h5)
#ADD_TEST(Compare:Legacy:Incompatible ${EXECUTABLE_OUTPUT_PATH}/pc2compare HDF5Examples/MultipleGrids2D.h5 HDF5Examples/MultipleGrids3D.h5)
#ADD_TEST(Compare:PC2:Same ${EXECUTABLE_OUTPUT_PATH}/pc2compare HDF5Examples/pc2_2d_0.h5 HDF5Examples/pc2_2d_0.h5)
#ADD_TEST(Compare:PC2:Different ${EXECUTABLE_OUTPUT_PATH}/pc2compare HDF5Examples/pc2_2d_1000.h5 HDF5Examples/pc2_2d_0.h5)
#ADD_TEST(Compare:PC2:Incompatible ${EXECUTABLE_OUTPUT_PATH}/pc2compare HDF5Examples/MultipleGrids2D.h5 HDF5Examples/pc2_2d_1000.h5)
#ADD_TEST(Compare:LegacyPC2:Same ${EXECUTABLE_OUTPUT_PATH}/pc2compare HDF5Examples/pc2_2d_1000.h5 InviscidTests/AcousticPulsePeriodic2D/CM/flocm_acousticpulse_1000.h5 1e-4)
#ADD_TEST(Compare:LegacyPC2:Different ${EXECUTABLE_OUTPUT_PATH}/pc2compare HDF5Examples/pc2_2d_1000.h5 PlasComCM/hdf5/flocm_2d_005.h5)
#ADD_TEST(Compare:LegacyPC2:Incompatible ${EXECUTABLE_OUTPUT_PATH}/pc2compare HDF5Examples/pc2_2d_1000.h5 InviscidTests/AcousticPulsePeriodic3D/CM/flocm_acousticpulse_5.h5)

#SET_TESTS_PROPERTIES(Compare:Legacy:Different Compare:Legacy:Incompatible 
#                     Compare:PC2:Different Compare:PC2:Incompatible 
#                     Compare:LegacyPC2:Different Compare:LegacyPC2:Incompatible PROPERTIES WILL_FAIL true)

#ADD_TEST(ComparisonUtil:Legacy:Same ${EXECUTABLE_OUTPUT_PATH}/pc2compare legacyCMFile.h5 legacyCMFile.h5)
#ADD_TEST(ComparisonUtil:Legacy:Different ${EXECUTABLE_OUTPUT_PATH}/pc2compare legacyCMFile.h5 legacyCMFile.h5)
#ADD_TEST(ComparisonUtil:Different pc2compare legacyCMFile.h5 )



#
# discover and add serial unit tests
#
SET(numSerialTests 0)
MESSAGE(STATUS "${Blue}Discovering serial unit tests ...${ColorReset}")
EXECUTE_PROCESS(COMMAND grep -I -r --exclude=CMakeLists.txt --exclude=*~ --exclude=\#* "serialUnitResults.UpdateResult" "${CMAKE_CURRENT_SOURCE_DIR}/src"
                COMMAND grep -v "//"
                COMMAND cut -d "(" -f 2
                COMMAND cut -d "," -f 1
                OUTPUT_VARIABLE allSerialUnitTests
                )
#message(${allSerialUnitTests})
STRING(REPLACE "\"" "" allSerialUnitTests ${allSerialUnitTests})
STRING(REGEX REPLACE "\n" " " allSerialUnitTests ${allSerialUnitTests})
SET(serialUnitTestList ${allSerialUnitTests})
SEPARATE_ARGUMENTS(serialUnitTestList)
LIST(SORT serialUnitTestList)
LIST(REMOVE_DUPLICATES serialUnitTestList)

FILE(WRITE ${CMAKE_SOURCE_DIR}/utils/EnabledSerialTests.txt "")

FOREACH(indTest ${serialUnitTestList})
  # only add tests from the active suites
  FOREACH(suiteName ${serialTestSuiteList})
    IF("${indTest}" MATCHES "^${suiteName}:")
      IF(((DISABLE_FAILING_TESTS) AND ("${indTest}" MATCHES "FAILS")))
        IF(SHOWTESTS)
          MESSAGE(STATUS "${Yellow}\t excluding failing test ${indTest}${ColorResest}")
        ENDIF()
      ELSEIF((NOT USE_CANTERA) AND ("${indTest}" MATCHES "Cantera"))
        IF(SHOWTESTS)
          MESSAGE(STATUS "${Yellow}\t excluding Cantera test ${indTest}${ColorResest}")
        ENDIF()
      ELSEIF((NOT ENABLE_PAPI) AND ("${indTest}" MATCHES "PAPIX"))
        IF(SHOWTESTS)
          MESSAGE(STATUS "${Yellow}\t excluding PAPI interface test ${indTest}${ColorResest}")
        ENDIF()
      ELSEIF((NOT ENABLE_Y6_TESTS) AND ("${indTest}" MATCHES "Y6"))
        IF(SHOWTESTS)
          MESSAGE(STATUS "${Yellow}\t excluding Y6 test ${indTest}${ColorResest}")
        ENDIF()
#      ELSEIF((NOT USE_CHEMGEN) AND ("${indTest}" MATCHES "Prometheus"))
#        IF(SHOWTESTS)
#          MESSAGE(STATUS "${Yellow}\t excluding Prometheus test ${indTest}${ColorResest}")
#        ENDIF()
      ELSE()
        IF(SHOWTESTS) 
	  message(STATUS "${Cyan}\t adding test ${indTest}${ColorReset}")
        ENDIF()
        ADD_TEST(${indTest} ${TESTRESULT} ${indTest} plascom2_testresults.txt)
        FILE(APPEND ${CMAKE_SOURCE_DIR}/utils/EnabledSerialTests.txt "${indTest}\n")
        MATH(EXPR numSerialTests "${numSerialTests} + 1")
      ENDIF()
    ENDIF()
  ENDFOREACH()
ENDFOREACH()
MESSAGE(STATUS "${Blue}Found ${numSerialTests} tests in ${numSerialSuites} serial suites.${ColorReset}")


#
# discover and add parallel unit tests
#
SET(numParallelTests 0)
MESSAGE(STATUS "${Blue}Discovering parallel unit tests ...${ColorReset}")
EXECUTE_PROCESS(COMMAND grep -I -r --exclude=CMakeLists.txt --exclude=*~ --exclude=\#* "parallelUnitResults.UpdateResult" "${CMAKE_CURRENT_SOURCE_DIR}/src"
                COMMAND grep -v "//"
                COMMAND cut -d "(" -f 2
                COMMAND cut -d "," -f 1
                OUTPUT_VARIABLE allParallelUnitTests
                )
string(REPLACE "\"" "" allParallelUnitTests ${allParallelUnitTests})
string(REGEX REPLACE "\n" " " allParallelUnitTests ${allParallelUnitTests})
set(parallelUnitTestList ${allParallelUnitTests})
separate_arguments(parallelUnitTestList)
list(SORT parallelUnitTestList)
list(REMOVE_DUPLICATES parallelUnitTestList)

FILE(WRITE ${CMAKE_SOURCE_DIR}/utils/EnabledParallelTests.txt "")
FOREACH(indTest ${parallelUnitTestList})
  # only add tests from the active suites
  FOREACH(suiteName ${parallelTestSuiteList})
    IF("${indTest}" MATCHES "^${suiteName}:")
      IF((DISABLE_FAILING_TESTS) AND ("${indTest}" MATCHES "FAILS"))
        IF(SHOWTESTS)
          MESSAGE(STATUS "${Yellow}\t excluding failing test ${indTest}${ColorResest}")
        ENDIF()
      ELSEIF((NOT USE_CANTERA) AND ("${indTest}" MATCHES "Cantera"))
        IF(SHOWTESTS)
          MESSAGE(STATUS "${Yellow}\t excluding Cantera test ${indTest}${ColorResest}")
        ENDIF()
      ELSEIF((NOT ENABLE_Y6_TESTS) AND ("${indTest}" MATCHES "Y6"))
        IF(SHOWTESTS)
          MESSAGE(STATUS "${Yellow}\t excluding Y6 test ${indTest}${ColorResest}")
        ENDIF()
#      ELSEIF((NOT USE_CHEMGEN) AND ("${indTest}" MATCHES "Prometheus"))
#        IF(SHOWTESTS)
#          MESSAGE(STATUS "${Yellow}\t excluding Prometheus test ${indTest}${ColorResest}")
#        ENDIF()
      ELSE()
        IF(SHOWTESTS) 
          message(STATUS "${Cyan}\t adding test ${indTest}${ColorReset}")
        ENDIF()
        add_test(${indTest} ${TESTRESULT} ${indTest} plascom2_testresults.txt)
        FILE(APPEND ${CMAKE_SOURCE_DIR}/utils/EnabledParallelTests.txt "${indTest}\n")
        MATH(EXPR numParallelTests "${numParallelTests} + 1")
      ENDIF()
    ENDIF()
  ENDFOREACH()
ENDFOREACH()
MESSAGE(STATUS "${Blue}Found ${numParallelTests} tests in ${numParallelSuites} parallel suites.${ColorReset}")

#INCLUDE(parallel_unit_tests.cmake)
