cmake_minimum_required(VERSION 3.16)

set(LIBRARY_NAME
    math_lib
)

set(GLOBAL_NAMESPACE
    xitren
)

project(${LIBRARY_NAME}
	DESCRIPTION "Math library"
    VERSION 1.0.1
	LANGUAGES CXX
)

add_library(${LIBRARY_NAME} INTERFACE)
add_library(${GLOBAL_NAMESPACE}::${LIBRARY_NAME} ALIAS ${LIBRARY_NAME})

target_include_directories(
        ${LIBRARY_NAME} INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>)

message(STATUS "Adding library project \"${LIBRARY_NAME}\"")

include(FetchContent)
FetchContent_Declare(
		circular_buffer
		GIT_REPOSITORY https://github.com/xitren/circular_buffer.git
		GIT_TAG        v1.0.1
)
FetchContent_MakeAvailable(circular_buffer)
target_link_libraries(${LIBRARY_NAME} INTERFACE ${GLOBAL_NAMESPACE}::circular_buffer)

enable_testing()
add_subdirectory(tests)