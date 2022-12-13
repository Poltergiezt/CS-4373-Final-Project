cmake_minimum_required(VERSION 3.23)
project(untitled C)

set(CMAKE_C_STANDARD 99)

include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories(SYSTEM ${MSMPI_INC})
include_directories(SYSTEM ${MSMPI_INC}/x64)
include_directories(SYSTEM ${MSMPI_BIN})
include_directories(SYSTEM ${MSMPI_LIB64})
find_package(MPI REQUIRED)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executable(Part1MPI DeterminantMPI.c)

add_executable(Part1OMP DeterminantOMP.c)


add_executable(Part1 Determinant.c)

add_executable(Part2SA TSP_SA.c)

add_executable(Part2SAMPI TSP_SA_MPI.c)



target_link_libraries(Part1MPI PUBLIC MPI::MPI_C)
target_link_libraries(Part2SAMPI PUBLIC MPI::MPI_C)