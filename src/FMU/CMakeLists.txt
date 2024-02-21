# SPDX-FileCopyrightText: 2023 SAP SE
#
# SPDX-License-Identifier: Apache-2.0
#
# This file is part of FEDEM - https://openfedem.org

# Library setup

set ( LIB_ID fedem_fmu )
set ( UNIT_ID ${DOMAIN_ID}_${PACKAGE_ID}_${LIB_ID} )

message ( STATUS "INFORMATION : Processing unit ${UNIT_ID}" )

string ( APPEND CMAKE_CXX_FLAGS_DEBUG " -DFMU_DEBUG" )
if ( FTENV_WARNINGS AND LINUX AND NOT USE_INTEL_FORTRAN )
  string ( APPEND CMAKE_CXX_FLAGS " -Wno-unused-parameter" )
endif ( FTENV_WARNINGS AND LINUX AND NOT USE_INTEL_FORTRAN )

# Build and install

include_directories ( "${CMAKE_CURRENT_SOURCE_DIR}/FMI2/headers" )
add_library ( ${LIB_ID}2 SHARED fmu2.C )
include_directories ( "${CMAKE_CURRENT_SOURCE_DIR}/FMI3/headers" )
add_library ( ${LIB_ID}3 SHARED fmu3.C )
add_dependencies ( all_solvers ${LIB_ID}2 ${LIB_ID}3 )

set ( CLOUD_DIR "${CMAKE_INSTALL_PREFIX}/bin/Templates/cloudsim" )
install ( TARGETS ${LIB_ID}2 ${LIB_ID}3
          RUNTIME DESTINATION "${CLOUD_DIR}"
          LIBRARY DESTINATION "${CLOUD_DIR}" )
