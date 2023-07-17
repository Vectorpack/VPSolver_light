cmake_minimum_required(VERSION 3.2)
project("VPSolver_light"
    VERSION 3.1.3 # Follows version from original VPSolver
    DESCRIPTION "Light C++ VPSolver library")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Where to find the Find*.cmake files
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

set(HEADER_FILES
    src/arcflow.hpp
    src/arcflowsol.hpp
    src/common.hpp
    src/graph.hpp
    src/instance.hpp
    src/vps_utils.hpp
)
set(SOURCE_FILES
    src/arcflow.cpp
    src/arcflowsol.cpp
    src/common.cpp
    src/graph.cpp
    src/instance.cpp
    src/vps_utils.cpp
)

# Gurobi is required
find_package(GUROBI REQUIRED)


# Creation of the library
add_library(${PROJECT_NAME} STATIC ${HEADER_FILES} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME}
    PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    PUBLIC ${GUROBI_INCLUDE_DIRS}

    PUBLIC $<INSTALL_INTERFACE:lib>
    PUBLIC $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
)

target_link_libraries(${PROJECT_NAME}
        PUBLIC ${GUROBI_LIBRARIES}
)

set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${HEADER_FILES}")

### Instal rules
install(TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}Config
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include/${PROJECT_NAME})

install(EXPORT ${PROJECT_NAME}Config DESTINATION lib/cmake/${PROJECT_NAME})
