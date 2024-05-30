cmake_minimum_required(VERSION 3.16.0)
project(BeBOT_cpp VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

# Add include directories
include_directories(${CMAKE_SOURCE_DIR}/include
                    ${CMAKE_SOURCE_DIR}/../ipopt/include
                    ${CMAKE_SOURCE_DIR}/../ipopt/src/Interfaces
                    /opt/intel/oneapi/mkl/latest/include
)

# List all source files for the BeBOT library
set(BEBOT_SOURCE_FILES
    ${CMAKE_SOURCE_DIR}/bebot/bebot.cpp
    ${CMAKE_SOURCE_DIR}/bebot/bernsteinpoly.cpp
    ${CMAKE_SOURCE_DIR}/bebot/degelevmatrix.cpp
    ${CMAKE_SOURCE_DIR}/bebot/bernsteinmatrix_a2b.cpp
    ${CMAKE_SOURCE_DIR}/bebot/bernsteindifferentialmatrix.cpp
    ${CMAKE_SOURCE_DIR}/bebot/nchoosek_mod.cpp
    ${CMAKE_SOURCE_DIR}/examples/bebot/ma_57/example_2/libbebot.cpp)

# Add the library target
add_library(bebot SHARED ${BEBOT_SOURCE_FILES})

# Set the output directory for the shared library
set_target_properties(bebot PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/bebot/ma_57/example_2
)

# Add the executable target
add_executable(BeBOT_cpp ${CMAKE_SOURCE_DIR}/examples/bebot/ma_57/example_2/libbebot_caller.cpp)

# Set C++ standard
set_property(TARGET BeBOT_cpp PROPERTY CXX_STANDARD 17)

# Specify the MKL libraries
set(MKL_LIBRARIES
    /opt/intel/oneapi/mkl/latest/lib/libmkl_intel_lp64.so
    /opt/intel/oneapi/mkl/latest/lib/libmkl_sequential.so
    /opt/intel/oneapi/mkl/latest/lib/libmkl_core.so
    /opt/intel/oneapi/mkl/latest/lib/libmkl_intel_thread.so
    /opt/intel/oneapi/mkl/latest/lib/libmkl_rt.so
    /opt/intel/oneapi/mkl/latest/lib/libmkl_avx2.so.2
    /opt/intel/oneapi/mkl/latest/lib/libmkl_def.so.2
    /opt/intel/oneapi/2024.1/lib/libiomp5.so)

# Specify the Ipopt libraries
set(IPOPT_LIBRARIES /home/camilla/dev_cpp/ipopt/lib/libipopt.so
                    /home/camilla/dev_cpp/ipopt/lib/libsipopt.so
                    /home/camilla/dev_cpp/ipopt/lib/libipoptamplinterface.so)

# Specify system libraries
set(SYSTEM_LIBRARIES pthread)

# Find the HSL library
find_library(HSL_LIBRARY NAMES hsl PATHS /usr/local/lib)

if(NOT HSL_LIBRARY)
    message(FATAL_ERROR "HSL library not found")
endif()

# Set RPATH to include the directories of libraries
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_SOURCE_DIR}/examples/bebot/ma_57/example_2;/usr/local/lib")

# Link the MKL, Ipopt, HSL, and system libraries to the library
target_link_libraries(bebot PRIVATE ${MKL_LIBRARIES} ${IPOPT_LIBRARIES} ${HSL_LIBRARY} ${SYSTEM_LIBRARIES})

target_link_libraries(BeBOT_cpp PRIVATE dl)

# Specify the output binary for the executable
set_target_properties(BeBOT_cpp PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_NAME test_bebot
)
