cmake_minimum_required(VERSION 3.5.0)
project(LaikaBxDFs VERSION 1.0)

# install directly to the project folder
SET(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/release")

# Please update variables for your OS
if (MSVC)
    # WINDOWS
    SET(HFS "C:/Program Files/Side Effects Software/Houdini 19.5.605" CACHE PATH "Houdini Installation Folder")
    SET(RMANTREE "C:/Program Files/Pixar/RenderManProServer-25.2" CACHE PATH "RendermanProServer Installation Folder")
    SET(RFHTREE "C:/Program Files/Pixar/RenderManForHoudini-25.2" CACHE PATH "RendermanForHoudini Installation Folder")
    SET(RMAN_EXAMPLES "C:/Pixar/PixarRenderMan-Examples-25.2" CACHE PATH "Renderman Examples need to be installed separately")
else()
    # LINUX
    SET(HFS "/opt/hfs19.5.640" CACHE PATH "Houdini Installation Folder")
    SET(RMANTREE "/opt/pixar/RenderManProServer-25.2" CACHE PATH "RendermanProServer Installation Folder")
    SET(RFHTREE "/opt/pixar/RenderManForHoudini-25.2" CACHE PATH "RendermanForHoudini Installation Folder")
    SET(RMAN_EXAMPLES "/opt/pixar/PixarRenderMan-Examples-25.2" CACHE PATH "Renderman Examples need to be installed separately")

endif()

# prevent attempt to build right in the project
if ("${PROJECT_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message (FATAL_ERROR "Not allowed to run in-source build!")
endif ()

# required to build HDAs
find_package(Houdini
    PATHS ${HFS}/toolkit/cmake)

find_program(ARGS2HDA args2hda.py
    HINTS ${RFHTREE}/${_houdini_python_version}/${Houdini_VERSION}/python${_houdini_python_version}libs)
if(ARGS2HDA_NOTFOUND)
    message(FATAL_ERROR "args2hda.py not found!")
endif()


include_directories(SYSTEM ${RMANTREE}/include ${RMAN_EXAMPLES}/plugins ${PROJECT_SOURCE_DIR}/cpp)
add_definitions(-DNDEBUG)

# Linking only on Windows
# On Unix-based platforms, plugins are linked such that symbols that resolve
# to entrypoints inlibprman.so or libprman.dylib are left unresolved.
if (MSVC)
    link_directories( ${RMANTREE}/lib )
    link_libraries(libprman libpxrcore libstats)
endif()


#################
# Build Plugins #
#################

set (bxdfs
    HGScatter
    KPVelvet
    PraterFuzz
    PraterScatter)

foreach (bxdf ${bxdfs})
    add_library( ${bxdf} SHARED cpp/bxdf/${bxdf}.cpp)
    set_target_properties(${bxdf} PROPERTIES PREFIX "")
endforeach ()

if (MSVC)
    #Prevent installation .lib files
    install(TARGETS ${bxdfs}
            RUNTIME COMPONENT Runtime
            DESTINATION bxdf)
else()
    install(TARGETS ${bxdfs}
            DESTINATION bxdf)
endif()

install(DIRECTORY
    ${PROJECT_SOURCE_DIR}/cpp/bxdf/Args
    DESTINATION bxdf)


#############
# Build HDA #
#############

file (GLOB args_files "cpp/bxdf/Args/*.args")

file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/otls)
if (MSVC)
    add_custom_command(OUTPUT LaikaBxDFs.hda
        COMMAND
        ${CMAKE_COMMAND} -E env
        "PYTHONPATH=${RFHTREE}/${_houdini_python_version}/${Houdini_VERSION}/python${_houdini_python_version}libs;${RMANTREE}/bin;${RMANTREE}/lib/python${_houdini_python_version}/Lib/site-packages"
        ${HFS}/bin/hython.exe ${ARGS2HDA} -O -m "Laika" -L ${CMAKE_CURRENT_BINARY_DIR}/LaikaBxDFs.hda ${args_files}
        DEPENDS ${args_files} )
else()
    add_custom_command(OUTPUT LaikaBxDFs.hda
        COMMAND ${CMAKE_COMMAND} -E env
        PATH="${HFS}/bin:$ENV{PATH}"
        PYTHONPATH="${RFHTREE}/${_houdini_python_version}/${Houdini_VERSION}/python${_houdini_python_version}libs:${RMANTREE}/bin:${RMANTREE}/lib/python${_houdini_python_version}/site-packages"
        hython ${ARGS2HDA} -O -m "Laika" -L ${CMAKE_CURRENT_BINARY_DIR}/LaikaBxDFs.hda ${args_files}
        DEPENDS ${args_files} )
endif()

add_custom_target(HDA ALL
                DEPENDS LaikaBxDFs.hda
                SOURCES ${args_files})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/LaikaBxDFs.hda
        DESTINATION otls)
