cmake_minimum_required(VERSION 3.5)

project(engine)

############################################################
# Set compiler flags
############################################################
# If supported by your compiler, you can add the -Wall, -Wextra, –fstack-protector-all and -g3 flags here.
set(OWN_GXX_FLAGS "-std=c++17")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${OWN_GXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OWN_GXX_FLAGS}")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OWN_GXX_FLAGS}")

############################################################
# List all sources
############################################################
set(engine_sources
		easy_image.cc
		engine.cc
		ini_configuration.cc
		intro_exercises.cpp
		l_systems.cpp
		l_parser.cc
		Line2D.h
		lineDrawer.cpp
		Figure.cpp
		Transformations.cpp
		Projection.cpp
		lijntekeningen3D.cpp
		vector3d.cc
		draw_zbuffered_wireframes.cpp
		ZBuffer.cpp
		draw_triangulated_wireframes.cpp
		fractals3D.cpp
		Light.h
)

############################################################
# Create an executable
############################################################
set(exe_name "engine")
add_executable( ${exe_name} ${engine_sources} )
install( TARGETS ${exe_name} DESTINATION ${PROJECT_SOURCE_DIR}/ )