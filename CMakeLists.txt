cmake_minimum_required(VERSION 3.16.0)
project(BeBOT_cpp VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

# Add include directories
include_directories(${CMAKE_SOURCE_DIR}/include
                    ${CMAKE_SOURCE_DIR}/../ipopt/include
                    ${CMAKE_SOURCE_DIR}/../ipopt/src/Interfaces
                    /home/camilla/anaconda3/include)

# List all source files
set(SOURCE_FILES
    ${CMAKE_SOURCE_DIR}/examples/bebot/ma_57/example_2/libbebot_caller.cpp
    ${CMAKE_SOURCE_DIR}/bebot/bebot.cpp
    ${CMAKE_SOURCE_DIR}/bebot/bernsteinpoly.cpp
    ${CMAKE_SOURCE_DIR}/bebot/degelevmatrix.cpp
    ${CMAKE_SOURCE_DIR}/bebot/bernsteinmatrix_a2b.cpp
    ${CMAKE_SOURCE_DIR}/bebot/bernsteindifferentialmatrix.cpp
    ${CMAKE_SOURCE_DIR}/bebot/nchoosek_mod.cpp)

# Add the executable target
add_executable(BeBOT_cpp ${SOURCE_FILES})

# Set C++ standard
set_property(TARGET BeBOT_cpp PROPERTY CXX_STANDARD 17)

# Specify the MKL libraries
set(MKL_LIBRARIES /home/camilla/anaconda3/lib/libmkl_intel_lp64.so
                  /home/camilla/anaconda3/lib/libmkl_intel_lp64.so.2
                  /home/camilla/anaconda3/lib/libmkl_sequential.so
                  /home/camilla/anaconda3/lib/libmkl_core.so)

# Specify the Ipopt libraries
set(IPOPT_LIBRARIES /home/camilla/dev_cpp/ipopt/lib/libipopt.so
                    /home/camilla/dev_cpp/ipopt/lib/libsipopt.so
                    /home/camilla/dev_cpp/ipopt/lib/libipoptamplinterface.so)

# Specify the path to the BeBOT library
set(BEBOT_LIBRARIES ${CMAKE_SOURCE_DIR}/examples/bebot/ma_57/example_2/libbebot.so)

# Specify system libraries
set(SYSTEM_LIBRARIES dl pthread)

# Find the HSL library
find_library(HSL_LIBRARY NAMES hsl PATHS /usr/local/lib)

if(NOT HSL_LIBRARY)
    message(FATAL_ERROR "HSL library not found")
endif()

# Set RPATH to include the directories of libraries
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_SOURCE_DIR}/examples/bebot/ma_57/example_2;/usr/local/lib;/home/camilla/anaconda3/lib")

# Link the MKL, Ipopt, HSL, and system libraries to your target
target_link_libraries(BeBOT_cpp ${MKL_LIBRARIES} ${IPOPT_LIBRARIES} ${HSL_LIBRARY} ${SYSTEM_LIBRARIES} ${BEBOT_LIBRARIES})

# Specify the output binary
set_target_properties(BeBOT_cpp PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_NAME test_bebot
)
