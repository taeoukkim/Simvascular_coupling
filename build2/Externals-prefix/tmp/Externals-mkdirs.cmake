# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/home/runner/work/3D-1D-coupling-dev/3D-1D-coupling-dev/Externals")
  file(MAKE_DIRECTORY "/home/runner/work/3D-1D-coupling-dev/3D-1D-coupling-dev/Externals")
endif()
file(MAKE_DIRECTORY
  "/home/runner/work/3D-1D-coupling-dev/3D-1D-coupling-dev/build2/Externals-build"
  "/home/runner/work/3D-1D-coupling-dev/3D-1D-coupling-dev/build2/Externals-prefix"
  "/home/runner/work/3D-1D-coupling-dev/3D-1D-coupling-dev/build2/Externals-prefix/tmp"
  "/home/runner/work/3D-1D-coupling-dev/3D-1D-coupling-dev/build2/Externals-prefix/src/Externals-stamp"
  "/home/runner/work/3D-1D-coupling-dev/3D-1D-coupling-dev/build2/Externals-prefix/src"
  "/home/runner/work/3D-1D-coupling-dev/3D-1D-coupling-dev/build2/Externals-prefix/src/Externals-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/runner/work/3D-1D-coupling-dev/3D-1D-coupling-dev/build2/Externals-prefix/src/Externals-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/runner/work/3D-1D-coupling-dev/3D-1D-coupling-dev/build2/Externals-prefix/src/Externals-stamp${cfgdir}") # cfgdir has leading slash
endif()
