set(DS_SOURCES FL_table.hpp LF_table.hpp r_index.hpp)

add_library(ds OBJECT ${DS_SOURCES})
target_link_libraries(ds common sdsl ri)