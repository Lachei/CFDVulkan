##################################
##### Compiling the Shaders  #####
##################################

if (${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Windows")
    set(GLSL_VALIDATOR "$ENV{VULKAN_SDK}/Bin/glslangValidator.exe")
elseif (${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Linux")
    set(GLSL_VALIDATOR "$ENV{VULKAN_SDK}/bin/glslangValidator")
else ()
    message("Unsupported build platform")
endif ()

file(GLOB_RECURSE GLSL_SOURCE_FILES
        "*.frag"
        "*.vert"
        "*.comp"
)

foreach (GLSL_FILE ${GLSL_SOURCE_FILES})
    get_filename_component(FILE_NAME ${GLSL_FILE} NAME)
    set(SPIRV_FILE "${PROJECT_BINARY_DIR}/shaders/${FILE_NAME}.spv")
    message("${GLSL_VALIDATOR} -V --target-env spirv1.3 ${GLSL_FILE} -o ${SPIRV_FILE}")
    add_custom_command(
            OUTPUT ${SPIRV_FILE}
            COMMAND ${CMAKE_COMMAND} -E make_directory "${PROJECT_BINARY_DIR}/shaders/"
            COMMAND ${GLSL_VALIDATOR} -V --target-env spirv1.3 ${GLSL_FILE} -o ${SPIRV_FILE}
            DEPENDS ${GLSL_FILE})
    list(APPEND SPIRV_BINARY_FILES ${SPIRV_FILE})
endforeach (GLSL_FILE)

file(GLOB SHADER_HEADERS
    "*.h"
)

list(APPEND SPIRV_BINARY_FILES ${SHADER_HEADERS})

add_custom_target(
        build_shader ALL
        DEPENDS ${SPIRV_BINARY_FILES}
)