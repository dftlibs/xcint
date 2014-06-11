#
#  This file is based on FindOpenCL code
#
# - Try to find Accelerate Framework on OSX
# This module tries to find an Accelerate Framework implementation on your Mac system. 
#
# Once done this will define
#  ACCELERATE_FOUND        - system has Accelerate framework
#  ACCELERATE_INCLUDE_DIRS - the Accelerate framework include directory
#  ACCELERATE_LIBRARIES    - link these to use Accelerate framework
#

FIND_PACKAGE( PackageHandleStandardArgs )

IF (APPLE)

  FIND_LIBRARY(ACCELERATE_LIBRARIES Accelerate DOC "Accelerate Framework for OSX")
  FIND_PATH(ACCELERATE_INCLUDE_DIRS Accelerate/Accelerate.h DOC "Include for Accelerate Framework on OSX")

ENDIF (APPLE)

FIND_PACKAGE_HANDLE_STANDARD_ARGS( Accelerate DEFAULT_MSG ACCELERATE_LIBRARIES ACCELERATE_INCLUDE_DIRS )

MARK_AS_ADVANCED(
  ACCELERATE_INCLUDE_DIRS
)
