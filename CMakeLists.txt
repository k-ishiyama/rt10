cmake_minimum_required(VERSION 3.12)
project (app)
file (GLOB MAIN_SOURCES "*.h" "*.cpp")
add_executable(app ${MAIN_SOURCES})
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT app)

if(MSVC)
	# Properties > C/C++ > Output Files > Assembler Output
	target_compile_options(app PUBLIC "/Fa${PROJECT_BINARY_DIR}/${ASSEMBLY_LISTING_FILE}" /FAs)
	target_compile_features(app PUBLIC cxx_std_20)

	# OpenMP
	find_package(OpenMP REQUIRED)
  target_compile_options(app PUBLIC /openmp)
  target_link_libraries(app PUBLIC OpenMP::OpenMP_CXX)

	# Properties > Debugging > Working Directory
	set_target_properties(app PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
endif()

# target_link_libraries(${project_name} )
include_directories("${PROJECT_SOURCE_DIR}/glm")
