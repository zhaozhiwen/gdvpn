cmake_minimum_required(VERSION 2.6)

project(gdvpn)

set(LIBRARY_OUTPUT_PATH            "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH         "${CMAKE_BINARY_DIR}")
set(CMAKE_RUNTINE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

set(SOURCE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

# Compiler optimizations
set(CMAKE_CXX_FLAGS "-Ofast")

# Cuba library is used for integration
set(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -lcuba")

# Libraries in this package
set(LIBS "vmeson_params" "nnscatter" "dwf" "gdvpn")

# Compile libraries
foreach(LIB ${LIBS})
  add_library(${LIB} SHARED ${SOURCE_PATH}/${LIB}.cc)
endforeach()

# gdvpn depends on the other libraries
foreach(LIB ${LIBS})
  if(NOT ${LIB} MATCHES "gdvpn")
    target_link_libraries(gdvpn PUBLIC ${LIB})
  endif()
endforeach()

# Copy input data into build directory
set(TABLES "B_ampl_pn.tbl" "B_ampl_pp.tbl" "B_ampl_pp_nocoul.tbl")
foreach(TABLE ${TABLES})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${TABLE}
    ${CMAKE_BINARY_DIR}/${TABLE} COPYONLY)
endforeach()


# Compile and link the driver program
add_executable(driver.elf ${SOURCE_PATH}/driver.cc)
foreach(LIB ${LIBS})
  target_link_libraries(driver.elf ${LIB})
endforeach()
