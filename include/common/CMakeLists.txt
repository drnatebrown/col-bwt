set(COMMON_SOURCES common.hpp)

add_library(common OBJECT ${COMMON_SOURCES})
target_link_libraries(common sdsl)

set(IO_SOURCES io.hpp)

add_library(io OBJECT ${IO_SOURCES})
target_link_libraries(io klib)