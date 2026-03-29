cmake_minimum_required(VERSION 3.10)
project(mpc_controller)

# 1. Define your Source Files
# We explicitly list your C files. 
# We DO NOT list header (.h) files here; CMake finds them automatically.
set(SOURCES
    mpc_pipeline.c
    my_timing.c
)

set(CND
    condensedAB.c
    my_timing.c
)

set(HELI
    helicoptor.c
    my_timing.c
)
set(HELIN
    helicopter_lin.c
    my_timing.c
)
set(PEND
    invert_pendulum.c
    my_timing.c
)
set(CRANE
    crane_2D.c
    my_timing.c
)
set(LTV
    crane_LTV.c
    my_timing.c
)
# 2. Create the App
add_executable(mpc_test ${SOURCES})
add_executable(condensed ${CND})
add_executable(helicoptor ${HELI})
add_executable(helin ${HELIN})
add_executable(pendulum ${PEND})
add_executable(crane ${CRANE})
add_executable(ltv ${LTV})

# 3. Setup Header Paths (The "Spotlight")
# This replaces the need for any CMakeLists.txt inside the include folder
target_include_directories(mpc_test PRIVATE
    "${CMAKE_SOURCE_DIR}/include"
    "${CMAKE_SOURCE_DIR}/include/public"
)

target_include_directories(condensed PRIVATE
    "${CMAKE_SOURCE_DIR}/include"
    "${CMAKE_SOURCE_DIR}/include/public"
)

target_include_directories(helicoptor PRIVATE
    "${CMAKE_SOURCE_DIR}/include"
    "${CMAKE_SOURCE_DIR}/include/public"
)
target_include_directories(helin PRIVATE
    "${CMAKE_SOURCE_DIR}/include"
    "${CMAKE_SOURCE_DIR}/include/public"
)
target_include_directories(pendulum PRIVATE
    "${CMAKE_SOURCE_DIR}/include"
    "${CMAKE_SOURCE_DIR}/include/public"
)
target_include_directories(crane PRIVATE
    "${CMAKE_SOURCE_DIR}/include"
    "${CMAKE_SOURCE_DIR}/include/public"
)
target_include_directories(ltv PRIVATE
    "${CMAKE_SOURCE_DIR}/include"
    "${CMAKE_SOURCE_DIR}/include/public"
)

# 4. Find and Link the Library
# This looks for 'libosqp.a' inside your 'lib' folder
find_library(OSQP_LIB osqp 
    PATHS "${CMAKE_SOURCE_DIR}/lib"
    NO_DEFAULT_PATH
)

# Safety Check: Stop immediately if the library is missing
if(NOT OSQP_LIB)
    message(FATAL_ERROR "CRITICAL ERROR: 'libosqp' not found in mpc_example/lib/! \n Did you copy the compiled library there?")
endif()

# Link it!
target_link_libraries(mpc_test PRIVATE ${OSQP_LIB})
target_link_libraries(condensed PRIVATE ${OSQP_LIB})
target_link_libraries(helicoptor PRIVATE ${OSQP_LIB})
target_link_libraries(helin PRIVATE ${OSQP_LIB})
target_link_libraries(pendulum PRIVATE ${OSQP_LIB})
target_link_libraries(crane PRIVATE ${OSQP_LIB})
target_link_libraries(ltv PRIVATE ${OSQP_LIB})