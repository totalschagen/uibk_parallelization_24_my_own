add_library(IO 
    data_storage.cpp
)

target_include_directories(IO
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    PUBLIC ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(IO
    PUBLIC SerialHdf5
)