include(FetchContent)

## Add sdsl
FetchContent_Declare(
    sdsl
    # GIT_REPOSITORY https://github.com/elarielcl/sdsl-lite.git
    GIT_REPOSITORY https://github.com/simongog/sdsl-lite
)

FetchContent_GetProperties(sdsl)
if(NOT sdsl_POPULATED)
    FetchContent_Populate(sdsl)

    set(GENERATE_DOC OFF CACHE BOOL "Do not generate doxygen for sdsl-lite")
  
    add_subdirectory(${sdsl_SOURCE_DIR} ${sdsl_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

## Add klib
FetchContent_Declare(
  klib
  GIT_REPOSITORY https://github.com/drnatebrown/klib
)

FetchContent_GetProperties(klib)
if(NOT klib_POPULATED)
  FetchContent_Populate(klib)
  add_library(klib INTERFACE)

  target_include_directories(klib INTERFACE ${klib_SOURCE_DIR})
endif()

## Add r-index
FetchContent_Declare(
  r-index
  GIT_REPOSITORY https://github.com/maxrossi91/r-index.git
)

FetchContent_GetProperties(r-index)
if(NOT r-index_POPULATED)
  FetchContent_Populate(r-index)

  add_subdirectory(${r-index_SOURCE_DIR} ${r-index_BINARY_DIR} )#EXCLUDE_FROM_ALL)
  add_library(ri INTERFACE)
  target_link_libraries(ri INTERFACE klib z)
  target_include_directories(ri INTERFACE ${r-index_SOURCE_DIR}/internal)
endif()

## Add divsuffsort
FetchContent_Declare(
    divsufsort
    GIT_REPOSITORY https://github.com/simongog/libdivsufsort.git
    GIT_TAG        2.0.1
)

FetchContent_GetProperties(divsufsort)
if(NOT divsufsort_POPULATED)
    FetchContent_Populate(divsufsort)

    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Do not build a shared library for libdivsufsort")
    set(BUILD_EXAMPLES OFF CACHE BOOL "Do not build libdivsufsort example")
    set(BUILD_DIVSUFSORT64 ON CACHE BOOL "Build libdivsufsort in 64-bits mode")

    add_subdirectory(${divsufsort_SOURCE_DIR} ${divsufsort_BINARY_DIR} EXCLUDE_FROM_ALL)

    target_include_directories(divsufsort PUBLIC "${divsufsort_BINARY_DIR}/include")
    target_include_directories(divsufsort64 PUBLIC "${divsufsort_BINARY_DIR}/include")
endif()

## Add malloc_count
FetchContent_Declare(
    malloc_count
    GIT_REPOSITORY https://github.com/bingmann/malloc_count
)

FetchContent_GetProperties(malloc_count)
if(NOT malloc_count_POPULATED)
    FetchContent_Populate(malloc_count)

    add_library(malloc_count OBJECT ${malloc_count_SOURCE_DIR}/malloc_count.c ${malloc_count_SOURCE_DIR}/malloc_count.h)
    target_link_libraries(malloc_count dl)
    target_include_directories(malloc_count PUBLIC "${malloc_count_SOURCE_DIR}")

    add_library(memprofile OBJECT ${malloc_count_SOURCE_DIR}/memprofile.h)
    target_include_directories(memprofile PUBLIC "${malloc_count_SOURCE_DIR}")
endif()

# ## Add mumemto
# FetchContent_Declare(
#     mumemto
#     GIT_REPOSITORY https://github.com/drnatebrown/mumemto.git
#     GIT_TAG col_mums
# )
  
# FetchContent_GetProperties(mumemto)
# if(NOT mumemto_POPULATED)
#     FetchContent_Populate(mumemto)
#     # Ensure the mumemto project is configured and built correctly
#     set(MUMEMTO_BUILD_DIR ${mumemto_BINARY_DIR})
#     set(MUMEMTO_SOURCE_DIR ${mumemto_SOURCE_DIR})
#     add_subdirectory(${MUMEMTO_SOURCE_DIR} ${MUMEMTO_BUILD_DIR})
# endif()
# install(TARGETS mumemto
#     RUNTIME DESTINATION bin
#     LIBRARY DESTINATION lib
#     ARCHIVE DESTINATION lib
# )

# install(DIRECTORY ${mumemto_SOURCE_DIR}/include/ DESTINATION include)


# ## Add Movi
# FetchContent_Declare(
#     movi
#     GIT_REPOSITORY https://github.com/drnatebrown/Movi.git
#     GIT_TAG alt_reuse_byte
# )
  
# FetchContent_GetProperties(movi)
# if(NOT movi_POPULATED)
#     FetchContent_Populate(movi)
#     add_subdirectory(${movi_SOURCE_DIR} ${movi_BINARY_DIR})
# endif()