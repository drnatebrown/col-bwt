add_executable(build_FL build_FL.cpp)
target_link_libraries(build_FL sdsl malloc_count divsufsort divsufsort64 ri)
target_include_directories(build_FL PUBLIC    "../include/ds" 
                                        "../include/common"
                                        )
target_compile_options(build_FL PUBLIC "-std=c++17")

add_executable(col_split col_split.cpp)
target_link_libraries(col_split sdsl malloc_count divsufsort divsufsort64 ri)
target_include_directories(col_split PUBLIC    "../include/ds"
                                        "../include"
                                        "../include/common"
                                        )
target_compile_options(col_split PUBLIC "-std=c++17")

add_executable(build_col_bwt build_col_bwt.cpp)
target_link_libraries(build_col_bwt sdsl malloc_count divsufsort divsufsort64 ri)
target_include_directories(build_col_bwt PUBLIC    "../include/ds"
                                        "../include"
                                        "../include/common"
                                        )
target_compile_options(build_col_bwt PUBLIC "-std=c++17")

add_executable(pml_query pml_query.cpp)
target_link_libraries(pml_query sdsl malloc_count divsufsort divsufsort64 ri pthread klib)
target_include_directories(pml_query PUBLIC    "../include/ds"
                                        "../include"
                                        "../include/common"
                                        )
target_compile_options(pml_query PUBLIC "-std=c++17")

add_executable(rlbwt_to_bwt rlbwt_to_bwt.cpp)
target_link_libraries(rlbwt_to_bwt sdsl malloc_count divsufsort divsufsort64 ri)
target_include_directories(rlbwt_to_bwt PUBLIC    "../include/ds"
                                        "../include"
                                        "../include/common"
                                        )
target_compile_options(rlbwt_to_bwt PUBLIC "-std=c++17")