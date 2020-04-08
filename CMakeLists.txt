
CMAKE_MINIMUM_REQUIRED(VERSION 3.13)

set(PROJECT_NAME "LocallyInjectiveParametrization")
project( ${PROJECT_NAME} )

# User-specific assignments to variables
include(./settings_local.cmake OPTIONAL)

# cgal
find_package(CGAL COMPONENTS Core REQUIRED)
include( ${CGAL_USE_FILE} )

add_subdirectory(code)

set_property(TARGET ${PROJECT_NAME} PROPERTY FOLDER ${PROJECT_NAME})

# include dirs
target_include_directories( ${PROJECT_NAME} PRIVATE ${GMM_INCLUDE_DIR} )

# matlab
target_include_directories( ${PROJECT_NAME} PRIVATE ${MATLAB_DIR}/extern/include ${MATLAB_DIR}/simulink/include ${MATLAB_DIR}/extern/include/cpp )
set(MATLAB_STATIC_DIR ${MATLAB_DIR}/extern/lib/${MATLAB_PLATFORM}/microsoft)
target_link_libraries( ${PROJECT_NAME} # Matlab engine
    ${MATLAB_STATIC_DIR}/libeng.lib
    ${MATLAB_STATIC_DIR}/libmat.lib
    ${MATLAB_STATIC_DIR}/libmx.lib )
target_link_libraries( ${PROJECT_NAME} # Matlab BLAS/LAPACK
    ${MATLAB_STATIC_DIR}/libmwblas.lib
    ${MATLAB_STATIC_DIR}/libmwlapack.lib )

set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME})
