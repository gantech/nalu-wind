target_sources(nalu PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/AeroContainer.C
)

add_subdirectory(actuator)

add_subdirectory(aero_utils)

if(ENABLE_OPENFAST_FSI)
add_subdirectory(fsi)
endif()

