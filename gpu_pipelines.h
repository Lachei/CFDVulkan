#pragma once
#include <vulkan/vulkan.h>
#include "vk_utils.h"
#include "grid.hpp"

class GpuPipelines {
public:
    GpuPipelines(VkUtil::VkContext& context);
    ~GpuPipelines();

    struct Pipeline {
        VkPipeline pipeline;
        VkPipelineLayout pipelineLayout;
        VkDescriptorSetLayout descriptorSetLayout;
        VkDescriptorSet descriptorSet[2];   //every pipeline has two, as some pipelines have front and backbuffering
    };
    VkCommandBuffer simulate_commands[2];   // 0: update from back buffers to front buffer, 1: update form front to back buffer
    Pipeline calc_dt, calc_alpha, boundary_vals, calc_temp, calc_species, calc_fg, calc_rs, calc_sor, reset_sor, calc_uv;
    void init_descriptor_sets_and_commands(GpuGrid& grid, GpuSimulationInfos& simulation_info, SimulationContext& sim_ctx);
private:
    VkUtil::VkContext vk_context;
    VkSampler sampler;
    const char string_dt[30]= "shaders/calculate_dt.comp.spv";
    const char string_alpha[33]= "shaders/calculate_alpha.comp.spv";
    const char string_boundary[31]= "shaders/boundary_vals.comp.spv";
    const char string_temp[32]= "shaders/calculate_temp.comp.spv";
    const char string_species[35]= "shaders/calculate_species.comp.spv";
    const char string_fg[30]= "shaders/calculate_fg.comp.spv";
    const char string_rs[30]= "shaders/calculate_rs.comp.spv";
    const char string_sor[31]= "shaders/calculate_sor.comp.spv";
    const char string_sor_reset[31] = "shaders/reset_sor.comp.spv";
    const char string_uv[30]= "shaders/calculate_uv.comp.spv";
};