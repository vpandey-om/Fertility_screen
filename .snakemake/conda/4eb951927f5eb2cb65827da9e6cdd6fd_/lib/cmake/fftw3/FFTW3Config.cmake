# defined since 2.8.3
if (CMAKE_VERSION VERSION_LESS 2.8.3)
  get_filename_component (CMAKE_CURRENT_LIST_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)
endif ()

# Allows loading FFTW3 settings from another project
set (FFTW3_CONFIG_FILE "${CMAKE_CURRENT_LIST_FILE}")

set (FFTW3_LIBRARIES fftw3)
set (FFTW3_LIBRARY_DIRS /Users/vpandey/projects/githubs/Fertility_screen_final/Fertility_screen/.snakemake/conda/4eb951927f5eb2cb65827da9e6cdd6fd_/lib)
set (FFTW3_INCLUDE_DIRS /Users/vpandey/projects/githubs/Fertility_screen_final/Fertility_screen/.snakemake/conda/4eb951927f5eb2cb65827da9e6cdd6fd_/include)

include ("${CMAKE_CURRENT_LIST_DIR}/FFTW3LibraryDepends.cmake")

if (CMAKE_VERSION VERSION_LESS 2.8.3)
  set (CMAKE_CURRENT_LIST_DIR)
endif ()
