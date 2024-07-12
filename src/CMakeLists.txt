add_executable(build_FL build_FL.cpp)
target_link_libraries(build_FL sdsl malloc_count divsufsort divsufsort64 ri)
target_include_directories(build_FL PUBLIC    "../include/ds" 
                                        "../include/common"
                                        )
target_compile_options(build_FL PUBLIC "-std=c++17")

add_executable(col_split col_split.cpp)
target_link_libraries(col_split sdsl malloc_count divsufsort divsufsort64 ri)
target_include_directories(col_split PUBLIC    "../include/ds" 
                                        "../include/common"
                                        )
target_compile_options(col_split PUBLIC "-std=c++17")