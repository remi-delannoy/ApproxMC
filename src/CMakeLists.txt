# Copyright (c) 2017, Mate Soos
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

include_directories(
    ${PROJECT_SOURCE_DIR}
)

if (NOT WIN32)
    add_cxx_flag_if_supported("-Wno-bitfield-constant-conversion")
    #add_cxx_flag_if_supported("-Wduplicated-cond")
    #add_cxx_flag_if_supported("-Wduplicated-branches")
    add_cxx_flag_if_supported("-Wlogical-op")
    add_cxx_flag_if_supported("-Wrestrict")
    add_cxx_flag_if_supported("-Wnull-dereference")
    add_cxx_flag_if_supported("-Wjump-misses-init")
    add_cxx_flag_if_supported("-Wdouble-promotion")
    add_cxx_flag_if_supported("-Wshadow")
    add_cxx_flag_if_supported("-Wformat=2")
    add_cxx_flag_if_supported("-Wextra-semi")
    add_cxx_flag_if_supported("-pedantic")
    #add_cxx_flag_if_supported("-Wdeprecated")
endif()

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CRYPTOMINISAT5_INCLUDE_DIRS})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" @ONLY)

add_library(approxmc-lib
          approxmc.cpp)

add_executable(approxmc
    approxmc.cpp
    main.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp
)

set(approxmc_exec_link_libs
    ${Boost_LIBRARIES}
)

IF (ZLIB_FOUND)
    SET(approxmc_exec_link_libs ${approxmc_exec_link_libs} ${ZLIB_LIBRARY})
ENDIF()

target_link_libraries(approxmc-lib
    ${approxmc_exec_link_libs}
)
GENERATE_EXPORT_HEADER(approxmc-lib
         BASE_NAME approxmc
)
set_target_properties(approxmc-lib PROPERTIES
            PUBLIC_HEADER "${approxmc_public_headers}")

target_link_libraries(approxmc
    ${approxmc_exec_link_libs}
)

if (STATICCOMPILE)
    SET_TARGET_PROPERTIES(approxmc PROPERTIES LINK_SEARCH_START_STATIC 1)
    SET(cms_libs ${CRYPTOMINISAT5_STATIC_LIBRARIES})
    message(STATUS "bimm: ${CRYPTOMINISAT5_STATIC_LIBRARIES}")
else()
    SET(cms_libs ${CRYPTOMINISAT5_LIBRARIES})
endif()
set_target_properties(approxmc PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    INSTALL_RPATH_USE_LINK_PATH TRUE)
if (STATICCOMPILE)
    SET_TARGET_PROPERTIES(approxmc PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

target_link_libraries(approxmc-lib
    ${cms_libs}
)
target_link_libraries(approxmc
    ${cms_libs}
)

install(TARGETS approxmc
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

cmsat_add_public_header(approxmc-lib ${CMAKE_CURRENT_SOURCE_DIR}/approxmc.h)
cmsat_add_public_header(approxmc-lib ${CMAKE_CURRENT_SOURCE_DIR}/approxmcconfig.h)

install(TARGETS approxmc-lib
    EXPORT ${APPROXMC_EXPORT_NAME}
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/approxmc"
)

SET(CPACK_PACKAGE_EXECUTABLES "approxmc")
