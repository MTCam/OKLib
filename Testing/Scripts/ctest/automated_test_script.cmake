# ====== PLATFORM SPECIFIC SETTINGS =======
# THESE VARIABLES MUST BE SET FOR EACH PLATFORM ON WHICH TESTING IS ENABLED
# The following variables are inherited from the environment
# --------
# SVNCOMMAND = "/path/to/svn" 
# GITCOMMAND = "/path/to/git" 
# PROJECT_SOURCE = "/path/to/project/source" (e.g. /home/mtcampbe/AutomatedTesting/Codelet3)
# PROJECT_ROOT = "full/repository/path/to/project" (e.g. git@bitbucket.org:/xpacc-dev/Codelet3)
# These are picked up from the user's environment
# CMAKECOMMAND = "/home/mtcampbe/VTK/Install/bin/cmake"
# CTESTCOMMAND = "/home/mtcampbe/VTK/Install/bin/ctest"
# This setting is optional
# PROJECT_CONFIGURATION_OPTIONS = "-DBOOSTROOT=/path/to/boost"
# --------
#SET (_XPACC_REPO_ "svn://xpaccsvn/SourceRepository")
#SET (_XPACC_REPO_PUBLIC_ "https://bitbucket.org/xpacc-dev")
SET (_CONTINUOUS_INTERVAL_ $ENV{CONTINUOUS_INTERVAL})
SET (_CONTINUOUS_DURATION_ $ENV{CONTINUOUS_DURATION})
SET (_PATH_TO_SVN_ $ENV{SVNCOMMAND})
SET (_PATH_TO_GIT_ $ENV{GITCOMMAND})
SET (_PATH_TO_PROJECT_SOURCE_ $ENV{PROJECT_SOURCE})
SET (_PROJECT_ROOT_ $ENV{PROJECT_ROOT})
SET (_CTEST_COMMAND_ $ENV{CTESTCOMMAND})
SET (_CMAKE_COMMAND_ $ENV{CMAKECOMMAND})
SET (_PATH_TO_PROJECT_BINARY_ $ENV{PROJECT_BUILD_DIRECTORY})
SET (_REPO_TYPE_ $ENV{REPO_TYPE})
SET (_BRANCH_ $ENV{BRANCH})
SET (_PROJECT_CONFIG_TYPE_ $ENV{PROJECT_CONFIGURATION_TYPE})
#SET (_PROJECT_CONFIG_OPTIONS_ $ENV{PROJECT_CONFIGURATION_OPTIONS})
SET (_PROJECT_BUILD_OPTIONS_ $ENV{PROJECT_BUILD_OPTIONS})
# OPTIONAL
SET (_CMAKE_CONFIGURATION_OPTIONS_ $ENV{PROJECT_CONFIGURATION_OPTIONS})
# specific options related to building documentation
SET (_BUILD_DOCUMENTATION_ $ENV{BUILD_DOCS})
IF(${_BUILD_DOCUMENTATION_} MATCHES TRUE)
  SET (_DOCUMENTATION_TARGET_ $ENV{DOC_TARGET})
ENDIF(${_BUILD_DOCUMENTATION_} MATCHES TRUE)
# =========================================
find_program(MAKE NAMES make)
SET(CTEST_TEST_TIMEOUT 7200)

# -----------------------------------------------------------  
# -- Get environment
# -----------------------------------------------------------  
SET (CTEST_ENVIRONMENT
  "PATH=/home/mtcampbe/Software/Install/bin:${PATH}"
  "LD_LIBRARY_PATH=/home/mtcampbe/Software/Install/lib:${LD_LIBRARY_PATH}"
  )

## -- Set hostname
## --------------------------
find_program(HOSTNAME_CMD NAMES hostname)
macro(gethostname name flag)
  exec_program("${HOSTNAME_CMD}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
endmacro(gethostname)

gethostname(shorthost -s)
exec_program(${HOSTNAME_CMD} ARGS OUTPUT_VARIABLE HOSTNAME)
#exec_program(${HOSTNAME_CMD} -s OUTPUT_VARIABLE shorthost)
set(CTEST_SITE                          "${HOSTNAME}")
## -- Set site / build name
## --------------------------

find_program(UNAME NAMES uname)
macro(getuname name flag)
  exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
endmacro(getuname)

getuname(osname -s)
getuname(osrel  -r)
getuname(cpu    -m)

set(CTEST_BUILD_NAME                    "${osname}-${cpu}")

SET (CTEST_SOURCE_DIRECTORY "${_PATH_TO_PROJECT_SOURCE_}")
SET (CTEST_BINARY_DIRECTORY "${_PATH_TO_PROJECT_BINARY_}")
MESSAGE("Setting CTEST_SOURCE_DIRECTORY = ${CTEST_SOURCE_DIRECTORY}")

# should ctest wipe the binary tree before running
SET (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)

IF(${_REPO_TYPE_} MATCHES "svn")
  SET (CTEST_REPO_COMMAND "${_PATH_TO_SVN_}")
  SET (CTEST_REPO_CHECKOUT  "${CTEST_REPO_COMMAND} co ${_PROJECT_ROOT_} ${CTEST_SOURCE_DIRECTORY}")
ELSEIF(${_REPO_TYPE_} MATCHES "git")
  SET (CTEST_REPO_COMMAND "${_PATH_TO_GIT_}")
  SET (CTEST_REPO_CHECKOUT  "${CTEST_REPO_COMMAND} clone --recursive -b ${_BRANCH_} ${_PROJECT_ROOT_} ${CTEST_SOURCE_DIRECTORY}")
ELSE()
  MESSAGE(FATAL_ERROR "Unknown repository type, ${_REPO_TYPE_}, specified in projects file, exiting.")
ENDIF()

###
###
### FOREIGN CODE
###
##
## -- Checkout command
## -- Build Command
set(CTEST_BUILD_COMMAND                "${MAKE} ${_PROJECT_BUILD_OPTIONS_}")
###
### FOREIGN CODE
###
###



# which ctest command to use for running the dashboard
SET(MODEL Experimental)
IF(${CTEST_SCRIPT_ARG} MATCHES Nightly)
  SET(MODEL Nightly)
ENDIF(${CTEST_SCRIPT_ARG} MATCHES Nightly)
IF(${CTEST_SCRIPT_ARG} MATCHES Continuous)
  SET(MODEL Continuous)
#  SET(CTEST_CONTINUOUS_DURATION 60)
#  SET(CTEST_CONTINUOUS_MINIMUM_INTERVAL 10)
ENDIF(${CTEST_SCRIPT_ARG} MATCHES Continuous)
SET (CTEST_COMMAND "${_CTEST_COMMAND_} -D ${MODEL}")


find_program(CTEST_COVERAGE_COMMAND NAMES gcov)
find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
#SET(XPACC_SOURCE_REPOSITORY "file:///Projects/XPACC/SourceRepository")

if(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}")
    #set(CTEST_CHECKOUT_COMMAND "${CTEST_SVN_COMMAND} co ${_PROJECT_ROOT_} ${CTEST_SOURCE_DIRECTORY}")
    set(CTEST_CHECKOUT_COMMAND "${CTEST_REPO_CHECKOUT}")
endif()
SET(CTEST_UPDATE_COMMAND "${CTEST_REPO_COMMAND}")

# what cmake command to use for configuring this dashboard
SET (CTEST_CMAKE_COMMAND "${_CMAKE_COMMAND_}")

####################################################################
# The values in this section are optional you can either
# have them or leave them commented out
####################################################################

# this is the initial cache to use for the binary tree, be careful to escape
# any quotes inside of this string if you use it
SET (CTEST_INITIAL_CACHE "
CMAKE_GENERATOR:INTERNAL=Unix Makefiles
BUILDNAME:STRING=${CTEST_BUILD_NAME}
SITE:STRING=${CTEST_SITE}
SVN_UPDATE_OPTIONS:STRING=update
")

set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
IF(${_PROJECT_CONFIG_TYPE_} MATCHES "autotools")
  set(CTEST_CONFIGURE_COMMAND            "${CTEST_SOURCE_DIRECTORY}/configure ${_CMAKE_CONFIGURATION_OPTIONS_}")
#    configure_file(${CTEST_SOURCE_DIRECTORY}/CMake/CTestCustom.cmake   ${CTEST_BINARY_DIRECTORY}/CTestCustom.cmake)
#    ctest_read_custom_files("${CTEST_BINARY_DIRECTORY}")
ELSE()
  set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} ${_CMAKE_CONFIGURATION_OPTIONS_} ${CTEST_SOURCE_DIRECTORY}")
ENDIF()

#MESSAGE("\n*** printing some variables for debugging ***\n")
#MESSAGE(" CC ${CC} CXX ${CXX} FC ${FC} ")
#MESSAGE(" CC $ENV{CC} CXX $ENV{CXX} FC $ENV{FC} ")

IF(${MODEL} MATCHES Continuous)
  while (${CTEST_ELAPSED_TIME} LESS ${_CONTINUOUS_DURATION_})
    set (START_TIME ${CTEST_ELAPSED_TIME})
    ctest_start (Continuous)
    ctest_update (RETURN_VALUE count)
    IF(${_PROJECT_CONFIG_TYPE_} MATCHES "autotools")
      configure_file(${CTEST_SOURCE_DIRECTORY}/CMake/CTestTestfile.cmake ${CTEST_BINARY_DIRECTORY}/CTestTestfile.cmake)
#      configure_file(${CTEST_SOURCE_DIRECTORY}/CMake/CTestCustom.cmake   ${CTEST_BINARY_DIRECTORY}/CTestCustom.cmake)
      ctest_read_custom_files("${CTEST_BINARY_DIRECTORY}")
    ENDIF()
    ## -- Update git submodules
    IF (EXISTS "${CTEST_SOURCE_DIRECTORY}/.gitmodules")
      message (" -- Updating submodules --")
      execute_process (COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
	WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY})
      execute_process (COMMAND ${GIT_EXECUTABLE} submodule foreach git pull origin master
	WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY})
    ENDIF ()
    if (count GREATER 0)
      MESSAGE("\n****** Project Change Detected.... ******\n")
      MESSAGE("\n***** Configuring Project  *****\n")
      IF(${_PROJECT_CONFIG_TYPE_} MATCHES "autotools")
	
	message(" -- AutoGen ${MODEL} - ${CTEST_BUILD_NAME} --")
	message("CTEST_SOURCE_DIRECTORY = ${CTEST_SOURCE_DIRECTORY}")
	execute_process(COMMAND /bin/sh ./autogen.sh
	  WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY} 
	  RESULT_VARIABLE autogenResult 
	  OUTPUT_VARIABLE autogenLog 
	  ERROR_VARIABLE autogenLog)
	file(WRITE ${CTEST_BINARY_DIRECTORY}/Testing/autogen.log "${autogenLog}")
	
	message(" -- Autoreconf ${MODEL} - ${CTEST_BUILD_NAME} --")
	execute_process(COMMAND autoreconf -f -i
	  WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY} 
	  RESULT_VARIABLE autoreconfResult 
	  OUTPUT_VARIABLE autoreconfLog 
	  ERROR_VARIABLE autoreconfLog)
	file(WRITE ${CTEST_BINARY_DIRECTORY}/Testing/autoreconf.log "${autoreconfLog}")
	
	if( NOT ${autoreconfResult} )
	  
	  ## -- Configure
	  message(" -- Configure ${MODEL} - ${CTEST_BUILD_NAME} --")
	  ctest_configure(BUILD  "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
	  
	  ## -- BUILD
	  message(" -- Build ${MODEL} - ${CTEST_BUILD_NAME} --")
	  ctest_build(    BUILD  "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
	  
	  ## -- INSTALL
	  message(" -- Install ${MODEL} - ${CTEST_BUILD_NAME} --")
	  execute_process(COMMAND "${MAKE} install ${_PROJECT_BUILD_OPTIONS_}" WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY} 
	    RESULT_VARIABLE makeInstallResult OUTPUT_VARIABLE makeInstallLog ERROR_VARIABLE makeInstallLog)
	  file(WRITE ${CTEST_BINARY_DIRECTORY}/Testing/makeinstall.log "${makeInstallLog}")
	  
	  ## -- TEST
	  message(" -- Test ${MODEL} - ${CTEST_BUILD_NAME} --")
	  ctest_test(     BUILD  "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
	  
	endif( NOT ${autoreconfResult} )
      ELSE()
	ctest_configure()
	MESSAGE("\n***** Building Project  *****\n")
	ctest_build()
	MESSAGE("\n***** Running Project Tests *****\n")
	ctest_test()
      ENDIF()
      if (WITH_MEMCHECK AND CTEST_COVERAGE_COMMAND)
	MESSAGE("\n***** Checking Project Test Coverage *****\n")
	ctest_coverage()
      endif (WITH_MEMCHECK AND CTEST_COVERAGE_COMMAND)
      if (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
	ctest_memcheck()
      endif (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
      MESSAGE("\n***** Submiting Results to CDash *****\n")
      ctest_submit()
    endif ()
    ctest_sleep( ${START_TIME} ${_CONTINUOUS_INTERVAL_} ${CTEST_ELAPSED_TIME})
  endwhile()
ELSE() 
  ctest_start("${MODEL}")
  MESSAGE("\n***** Getting latest source from the repository *****\n")
  ctest_update()
  IF(${_PROJECT_CONFIG_TYPE_} MATCHES "autotools")
    configure_file(${CTEST_SOURCE_DIRECTORY}/CMake/CTestTestfile.cmake ${CTEST_BINARY_DIRECTORY}/CTestTestfile.cmake)
#    configure_file(${CTEST_SOURCE_DIRECTORY}/CMake/CTestCustom.cmake   ${CTEST_BINARY_DIRECTORY}/CTestCustom.cmake)
    ctest_read_custom_files("${CTEST_BINARY_DIRECTORY}")
  ENDIF()
  ## -- Update git submodules
  IF (EXISTS "${CTEST_SOURCE_DIRECTORY}/.gitmodules")
    message (" -- Updating submodules --")
    execute_process (COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
      WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY})
    execute_process (COMMAND ${GIT_EXECUTABLE} submodule foreach git pull origin master
      WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY})
  ENDIF ()
  IF(${_PROJECT_CONFIG_TYPE_} MATCHES "autotools")

    message(" -- AutoGen ${MODEL} - ${CTEST_BUILD_NAME} --")
    message("CTEST_SOURCE_DIRECTORY = ${CTEST_SOURCE_DIRECTORY}")
    execute_process(COMMAND /bin/sh ./autogen.sh
      WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY} 
      RESULT_VARIABLE autogenResult 
      OUTPUT_VARIABLE autogenLog 
      ERROR_VARIABLE autogenLog)
    file(WRITE ${CTEST_BINARY_DIRECTORY}/Testing/autogen.log "${autogenLog}")

    message(" -- Autoreconf ${MODEL} - ${CTEST_BUILD_NAME} --")
    execute_process(COMMAND autoreconf -f -i
      WORKING_DIRECTORY ${CTEST_SOURCE_DIRECTORY} 
      RESULT_VARIABLE autoreconfResult 
      OUTPUT_VARIABLE autoreconfLog 
      ERROR_VARIABLE autoreconfLog)
    file(WRITE ${CTEST_BINARY_DIRECTORY}/Testing/autoreconf.log "${autoreconfLog}")

    if( NOT ${autoreconfResult} )

      ## -- Configure
      message(" -- Configure ${MODEL} - ${CTEST_BUILD_NAME} --")
      ctest_configure(BUILD  "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)

      ## -- BUILD
      message(" -- Build ${MODEL} - ${CTEST_BUILD_NAME} --")
      ctest_build(    BUILD  "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)

      ## -- INSTALL
      message(" -- Install ${MODEL} - ${CTEST_BUILD_NAME} --")
      execute_process(COMMAND "${MAKE} install ${_PROJECT_BUILD_OPTIONS_}" WORKING_DIRECTORY ${CTEST_BINARY_DIRECTORY} 
	RESULT_VARIABLE makeInstallResult OUTPUT_VARIABLE makeInstallLog ERROR_VARIABLE makeInstallLog)
      file(WRITE ${CTEST_BINARY_DIRECTORY}/Testing/makeinstall.log "${makeInstallLog}")

      ## -- TEST
      message(" -- Test ${MODEL} - ${CTEST_BUILD_NAME} --")
      ctest_test(     BUILD  "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
    endif( NOT ${autoreconfResult} )
  ELSE()
    MESSAGE("\n***** Configuring *****\n")
    ctest_configure()
    MESSAGE("\n***** Building Project *****\n")
    ctest_build()
    IF(${_BUILD_DOCUMENTATION_} MATCHES TRUE)
      MESSAGE("\n***** Building documentation *****\n")
      ctest_build(TARGET "${_DOCUMENTATION_TARGET_}")
    ENDIF(${_BUILD_DOCUMENTATION_} MATCHES TRUE)
    MESSAGE("\n***** Running Project Tests *****\n")
    ctest_test()
  ENDIF(${_PROJECT_CONFIG_TYPE_} MATCHES "autotools")
  if (WITH_MEMCHECK AND CTEST_COVERAGE_COMMAND)
    MESSAGE("\n***** Checking Project Test Coverage *****\n")
    ctest_coverage()
  endif (WITH_MEMCHECK AND CTEST_COVERAGE_COMMAND)
  if (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
    ctest_memcheck()
  endif (WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
  MESSAGE("\n***** Submiting Results to CDash *****\n")
  ctest_submit()
ENDIF(${MODEL} MATCHES Continuous)