cmake_minimum_required(VERSION 3.1)
project(SoundTouch VERSION 2.3.0 LANGUAGES CXX)

include(GNUInstallDirs)

if(MSVC)
  set(COMPILE_DEFINITIONS CMAKE /O2 /fp:fast)
  set(COMPILE_OPTIONS )
else()
  set(COMPILE_OPTIONS -Ofast)
  set(COMPILE_DEFINITIONS CMAKE)
endif()

#####################
# SoundTouch library

add_library(SoundTouch
  source/SoundTouch/AAFilter.cpp
  source/SoundTouch/BPMDetect.cpp
  source/SoundTouch/cpu_detect_x86.cpp
  source/SoundTouch/FIFOSampleBuffer.cpp
  source/SoundTouch/FIRFilter.cpp
  source/SoundTouch/InterpolateCubic.cpp
  source/SoundTouch/InterpolateLinear.cpp
  source/SoundTouch/InterpolateShannon.cpp
  source/SoundTouch/mmx_optimized.cpp
  source/SoundTouch/PeakFinder.cpp
  source/SoundTouch/RateTransposer.cpp
  source/SoundTouch/SoundTouch.cpp
  source/SoundTouch/sse_optimized.cpp
  source/SoundTouch/TDStretch.cpp
)
target_include_directories(SoundTouch PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
   $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_definitions(SoundTouch PRIVATE ${COMPILE_DEFINITIONS})
target_compile_options(SoundTouch PRIVATE ${COMPILE_OPTIONS})

if(WIN32 AND BUILD_SHARED_LIBS)
  set_target_properties(SoundTouch PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  if(SOUNDTOUCH_DLL)
    target_compile_definitions(SoundTouchDLL PRIVATE DLL_EXPORTS)
  endif()
endif()

option(INTEGER_SAMPLES "Use integers instead of floats for samples" OFF)
if(INTEGER_SAMPLES)
  target_compile_definitions(SoundTouch PRIVATE SOUNDTOUCH_INTEGER_SAMPLES)
else()
  target_compile_definitions(SoundTouch PRIVATE SOUNDTOUCH_FLOAT_SAMPLES)
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(armv7.*|armv8.*)$")
  set(NEON_CPU ON)
else()
  set(NEON_CPU OFF)
endif()

option(NEON "Use ARM Neon SIMD instructions if in ARM CPU" ON)
if(${NEON} AND ${NEON_CPU})
  target_compile_definitions(SoundTouch PRIVATE SOUNDTOUCH_USE_NEON)
endif()

install(
  FILES
    include/BPMDetect.h
    include/FIFOSampleBuffer.h
    include/FIFOSamplePipe.h
    include/STTypes.h
    include/SoundTouch.h
  DESTINATION
    "${CMAKE_INSTALL_INCLUDEDIR}/soundtouch"
)

install(TARGETS SoundTouch
  EXPORT SoundTouchTargets
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

#######################
# soundstretch utility

add_executable(soundstretch
  source/SoundStretch/main.cpp
  source/SoundStretch/RunParameters.cpp
  source/SoundStretch/WavFile.cpp
)
target_include_directories(soundstretch PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_compile_definitions(soundstretch PRIVATE ${COMPILE_DEFINITIONS})
target_compile_options(soundstretch PRIVATE ${COMPILE_OPTIONS})
target_link_libraries(soundstretch PRIVATE SoundTouch)

install(TARGETS soundstretch
  DESTINATION bin
)

########################
# SoundTouchDll library

option(SOUNDTOUCH_DLL "Build SoundTouchDLL C wrapper library" OFF)
if(SOUNDTOUCH_DLL)
  add_library(SoundTouchDLL SHARED
    source/SoundTouchDLL/SoundTouchDLL.cpp
    source/SoundTouchDLL/SoundTouchDLL.rc
  )
  target_compile_options(SoundTouch PRIVATE -fPIC)
  target_compile_options(SoundTouchDLL PRIVATE -fPIC ${COMPILE_OPTIONS})
  target_compile_definitions(SoundTouchDLL PRIVATE ${COMPILE_DEFINITIONS})
  set_target_properties(SoundTouchDLL PROPERTIES CXX_VISIBILITY_PRESET hidden)
  target_include_directories(SoundTouchDLL PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
  target_link_libraries(SoundTouchDLL PRIVATE SoundTouch)
  install(FILES source/SoundTouchDLL/SoundTouchDLL.h DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/soundtouch")
  install(TARGETS SoundTouchDLL EXPORT SoundTouchTargets)
endif()

########################

# pkgconfig
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(execprefix "\${prefix}")
set(libdir "\${prefix}/${CMAKE_INSTALL_LIBDIR}")
set(includedir "\${prefix}/${CMAKE_INSTALL_INCLUDEDIR}")
set(VERSION "${CMAKE_PROJECT_VERSION}")
configure_file(soundtouch.pc.in "${CMAKE_CURRENT_BINARY_DIR}/soundtouch.pc" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/soundtouch.pc" DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

# CMake config
include(CMakePackageConfigHelpers)
set(SOUNDTOUCH_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/SoundTouch")
install(
  EXPORT SoundTouchTargets
  FILE SoundTouchTargets.cmake
  NAMESPACE SoundTouch::
  DESTINATION "${SOUNDTOUCH_INSTALL_CMAKEDIR}"
)
configure_package_config_file(SoundTouchConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/SoundTouchConfig.cmake"
  INSTALL_DESTINATION "${SOUNDTOUCH_INSTALL_CMAKEDIR}"
)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/SoundTouchConfigVersion.cmake"
  VERSION "${CMAKE_PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion
)
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/SoundTouchConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/SoundTouchConfigVersion.cmake"
  DESTINATION "${SOUNDTOUCH_INSTALL_CMAKEDIR}"
)
