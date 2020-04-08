
# Source Files
file(GLOB_RECURSE sources *.c *.cpp *.h *.hpp)

add_executable(${PROJECT_NAME} ${sources})

target_include_directories(${PROJECT_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/code")


