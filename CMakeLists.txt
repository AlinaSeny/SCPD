cmake_minimum_required(VERSION 3.26)
project(SKPD C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DN=1026m ,.1 -DREPETITIONS=1 -DDEBUG")

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
set(OMP_EXE_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

add_executable(SKPD source.c)
target_link_libraries(SKPD PRIVATE OpenMP::OpenMP_C)
set_target_properties(SKPD PROPERTIES COMPILE_FLAGS ${OMP_EXE_FLAGS})

add_executable(var36 var36.c)
target_link_libraries(var36 PRIVATE OpenMP::OpenMP_C)
set_target_properties(var36 PROPERTIES COMPILE_FLAGS ${OMP_EXE_FLAGS})

add_executable(var36_omp_for var36_omp_for.c)
target_link_libraries(var36_omp_for PRIVATE OpenMP::OpenMP_C)
set_target_properties(var36_omp_for PROPERTIES COMPILE_FLAGS ${OMP_EXE_FLAGS})

add_executable(var36_omp_task var36_omp_task.c)
target_link_libraries(var36_omp_task PRIVATE OpenMP::OpenMP_C)
set_target_properties(var36_omp_task PROPERTIES COMPILE_FLAGS ${OMP_EXE_FLAGS})

add_executable(var36_omp_taskloop var36_omp_taskloop.c)
target_link_libraries(var36_omp_taskloop PRIVATE OpenMP::OpenMP_C)
set_target_properties(var36_omp_taskloop PROPERTIES COMPILE_FLAGS ${OMP_EXE_FLAGS})

add_executable(var36_mpi var36_mpi.c)
target_link_libraries(var36_mpi ${MPI_LIBRARIES})

include_directories(${MPI_INCLUDE_PATH})
