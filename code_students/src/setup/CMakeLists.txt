add_library(sim_setup
    grid1D.cpp
    grid.cpp
    fluid.cpp
    physics.cpp
)

target_include_directories(sim_setup
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    PUBLIC ${PROJECT_SOURCE_DIR}/include
)