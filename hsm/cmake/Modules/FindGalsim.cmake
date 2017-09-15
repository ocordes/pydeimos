# - Try to find libgalsim
# Once done this will define
#  GALSIM_FOUND - System has libgalsim
#  GALSIM_INCLUDE_DIRS - The libgalsim include directories
#  GALSIM_LIBRARIES - The libraries needed to use LibXml2



#find_path(LIBXML2_INCLUDE_DIR libxml/xpath.h
#          HINTS ${PC_LIBXML_INCLUDEDIR} ${PC_LIBXML_INCLUDE_DIRS}
#          PATH_SUFFIXES libxml2 )


find_library(GALSIM_LIBRARY NAMES  galsim )


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GALSIM_FOUND to TRUE
# if all listed variables are TRUE
#find_package_handle_standard_args(GalSim  DEFAULT_MSG
                                #  GALSIM_LIBRARY GALSIM_INCLUDE_DIR)
find_package_handle_standard_args(GALSIM  DEFAULT_MSG GALSIM_LIBRARY)

mark_as_advanced(GALSIM_INCLUDE_DIR GALSIM_LIBRARY )

set(GALSIM_LIBRARIES ${GALSIM_LIBRARY} )
#set(LIBXML2_INCLUDE_DIRS ${LIBXML2_INCLUDE_DIR} )
