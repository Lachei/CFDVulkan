#include "gpu_pipelines.h"

GpuPipelines::GpuPipelines(VkUtil::VkContext& context):
    vk_context(context)
{
    VkUtil::create_image_sampler(vk_context, VK_SAMPLER_ADDRESS_MODE_BEGIN_RANGE, VK_FILTER_LINEAR, 0, 1, &sampler);
    // -----------------------------------------------------------------------------------------------------------------------
    // calc_dt pipeline
    // -----------------------------------------------------------------------------------------------------------------------
    VkShaderModule computeModule = VkUtil::create_shader_module(context, VkUtil::read_byte_file(string_dt));

    std::vector<VkDescriptorSetLayoutBinding> bindings;
    VkDescriptorSetLayoutBinding binding = {};
    binding.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;

    binding.binding = 0;								//velocities
    binding.descriptorCount = 1;
    binding.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
    bindings.push_back(binding);

    binding.binding = 1;								//sim_context
    binding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    bindings.push_back(binding);

    VkUtil::create_descriptor_set_layout(context, bindings, &calc_dt.descriptorSetLayout);
    std::vector<VkDescriptorSetLayout> layouts{ calc_dt.descriptorSetLayout };

    VkUtil::create_compute_pipeline(context, computeModule, layouts, &calc_dt.pipelineLayout, &calc_dt.pipeline);

    VkUtil::create_descriptor_sets(context, layouts, calc_dt.descriptorSet);
    // -----------------------------------------------------------------------------------------------------------------------
    // calc_alpha pipeline
    // -----------------------------------------------------------------------------------------------------------------------
    computeModule = VkUtil::create_shader_module(context, VkUtil::read_byte_file(string_alpha));

    VkUtil::create_descriptor_set_layout(context, bindings, &calc_alpha.descriptorSetLayout);
    layouts = { calc_alpha.descriptorSetLayout };

    VkUtil::create_compute_pipeline(context, computeModule, layouts, &calc_alpha.pipelineLayout, &calc_alpha.pipeline);

    VkUtil::create_descriptor_sets(context, layouts, calc_alpha.descriptorSet);
    // -----------------------------------------------------------------------------------------------------------------------
    // boundary_val pipeline
    // -----------------------------------------------------------------------------------------------------------------------
    computeModule = VkUtil::create_shader_module(context, VkUtil::read_byte_file(string_boundary));

    bindings.clear();

    binding.binding = 0;								//velocities
    binding.descriptorCount = 1;
    binding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
    bindings.push_back(binding);

    binding.binding = 1;                                //pressure
    bindings.push_back(binding);

    binding.binding = 2;                                //tempetature
    bindings.push_back(binding);

    binding.binding = 3;                                //species_concentrations
    bindings.push_back(binding);

    binding.binding = 4;                                //cell_flags
    bindings.push_back(binding);

    binding.binding = 5;                                //neighbours
    bindings.push_back(binding);

    binding.binding = 6;								//sim_context
    binding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    bindings.push_back(binding);
        
    binding.binding = 7;                                //inserter
    bindings.push_back(binding);

    binding.binding = 8;                                //converter
    bindings.push_back(binding);

    VkUtil::create_descriptor_set_layout(context, bindings, &boundary_vals.descriptorSetLayout);
    layouts = { boundary_vals.descriptorSetLayout, boundary_vals.descriptorSetLayout };

    VkUtil::create_compute_pipeline(context, computeModule, layouts, &boundary_vals.pipelineLayout, &boundary_vals.pipeline);

    VkUtil::create_descriptor_sets(context, layouts, boundary_vals.descriptorSet);
    // -----------------------------------------------------------------------------------------------------------------------
    // calc_temp pipeline
    // -----------------------------------------------------------------------------------------------------------------------
    computeModule = VkUtil::create_shader_module(context, VkUtil::read_byte_file(string_temp));

    bindings.clear();

    binding.binding = 0;								//velocities
    binding.descriptorCount = 1;
    binding.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
    bindings.push_back(binding);

    binding.binding = 1;                                //old_temperature
    bindings.push_back(binding);

    binding.binding = 2;                                //new_tempetature
    binding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
    bindings.push_back(binding);

    binding.binding = 3;								//sim_context
    binding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    bindings.push_back(binding);

    VkUtil::create_descriptor_set_layout(context, bindings, &calc_temp.descriptorSetLayout);
    layouts = { calc_temp.descriptorSetLayout, calc_temp.descriptorSetLayout };

    VkUtil::create_compute_pipeline(context, computeModule, layouts, &calc_temp.pipelineLayout, &calc_temp.pipeline);

    VkUtil::create_descriptor_sets(context, layouts, calc_temp.descriptorSet);
    // -----------------------------------------------------------------------------------------------------------------------
    // calc_species pipeline (includes normalization and density calculation)
    // -----------------------------------------------------------------------------------------------------------------------
    computeModule = VkUtil::create_shader_module(context, VkUtil::read_byte_file(string_species));

    bindings.clear();

    binding.binding = 0;								//velocities
    binding.descriptorCount = 1;
    binding.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
    bindings.push_back(binding);

    binding.binding = 1;                                //old_species
    bindings.push_back(binding);

    binding.binding = 2;                                //new_species
    binding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
    bindings.push_back(binding);

    binding.binding = 3;                                //temperature
    binding.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
    bindings.push_back(binding);

    binding.binding = 4;                                //density
    binding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
    bindings.push_back(binding);

    binding.binding = 5;                                //cell_flags
    binding.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
    bindings.push_back(binding);

    binding.binding = 6;								//sim_context
    binding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    bindings.push_back(binding);

    binding.binding = 7;                                //multispecies_converter
    bindings.push_back(binding);

    VkUtil::create_descriptor_set_layout(context, bindings, &calc_species.descriptorSetLayout);
    layouts = { calc_species.descriptorSetLayout, calc_species.descriptorSetLayout };

    VkUtil::create_compute_pipeline(context, computeModule, layouts, &calc_species.pipelineLayout, &calc_species.pipeline);

    VkUtil::create_descriptor_sets(context, layouts, calc_species.descriptorSet);
    // -----------------------------------------------------------------------------------------------------------------------
    // calc_fg pipeline
    // -----------------------------------------------------------------------------------------------------------------------
    computeModule = VkUtil::create_shader_module(context, VkUtil::read_byte_file(string_fg));

    bindings.clear();

    binding.binding = 0;								//velocities
    binding.descriptorCount = 1;
    //binding.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
    binding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
    bindings.push_back(binding);

    binding.binding = 1;                                //fg
    bindings.push_back(binding);

    binding.binding = 2;                                //density
    binding.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
    bindings.push_back(binding);

    binding.binding = 3;                                //cell_flags
    bindings.push_back(binding);

    binding.binding = 4;                                //neighbourhood
    bindings.push_back(binding);

    binding.binding = 5;								//sim_context
    binding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    bindings.push_back(binding);

    VkUtil::create_descriptor_set_layout(context, bindings, &calc_fg.descriptorSetLayout);
    layouts = { calc_fg.descriptorSetLayout };

    VkUtil::create_compute_pipeline(context, computeModule, layouts, &calc_fg.pipelineLayout, &calc_fg.pipeline);

    VkUtil::create_descriptor_sets(context, layouts, calc_fg.descriptorSet);
    // -----------------------------------------------------------------------------------------------------------------------
    // calc_rs pipeline
    // -----------------------------------------------------------------------------------------------------------------------
    computeModule = VkUtil::create_shader_module(context, VkUtil::read_byte_file(string_rs));

    bindings.clear();

    binding.binding = 0;								//rs
    binding.descriptorCount = 1;
    binding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
    bindings.push_back(binding);

    binding.binding = 1;                                //fg
    binding.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
    bindings.push_back(binding);

    binding.binding = 2;								//sim_context
    binding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    bindings.push_back(binding);

    VkUtil::create_descriptor_set_layout(context, bindings, &calc_rs.descriptorSetLayout);
    layouts = { calc_rs.descriptorSetLayout };

    VkUtil::create_compute_pipeline(context, computeModule, layouts, &calc_rs.pipelineLayout, &calc_rs.pipeline);

    VkUtil::create_descriptor_sets(context, layouts, calc_rs.descriptorSet);
    // -----------------------------------------------------------------------------------------------------------------------
    // calc_sor pipeline
    // -----------------------------------------------------------------------------------------------------------------------
    computeModule = VkUtil::create_shader_module(context, VkUtil::read_byte_file(string_sor));

    bindings.clear();

    binding.binding = 0;								//rs
    binding.descriptorCount = 1;
    binding.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
    bindings.push_back(binding);

    binding.binding = 1;                                //presssure
    binding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
    bindings.push_back(binding);

    binding.binding = 2;                                //cell_flags
    binding.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
    bindings.push_back(binding);

    binding.binding = 3;                                //neighbourhood
    bindings.push_back(binding);

    binding.binding = 4;								//sim_context
    binding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    bindings.push_back(binding);

    VkUtil::create_descriptor_set_layout(context, bindings, &calc_sor.descriptorSetLayout);
    layouts = { calc_sor.descriptorSetLayout };

    VkPushConstantRange push_constant;
    push_constant.offset = 0;
    push_constant.size = sizeof(uint32_t);
    push_constant.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
    std::vector<VkPushConstantRange> push_constants{ push_constant };

    VkUtil::create_compute_pipeline(context, computeModule, layouts, push_constants, &calc_sor.pipelineLayout, &calc_sor.pipeline);

    VkUtil::create_descriptor_sets(context, layouts, calc_sor.descriptorSet);
    // -----------------------------------------------------------------------------------------------------------------------
   // calc_sor pipeline
   // -----------------------------------------------------------------------------------------------------------------------
    computeModule = VkUtil::create_shader_module(context, VkUtil::read_byte_file(string_sor_reset));

    bindings.clear();

    binding.binding = 0;								//sim_context
    binding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    bindings.push_back(binding);

    VkUtil::create_descriptor_set_layout(context, bindings, &reset_sor.descriptorSetLayout);
    layouts = { reset_sor.descriptorSetLayout };

    VkUtil::create_compute_pipeline(context, computeModule, layouts, &reset_sor.pipelineLayout, &reset_sor.pipeline);

    VkUtil::create_descriptor_sets(context, layouts, reset_sor.descriptorSet);
    // -----------------------------------------------------------------------------------------------------------------------
    // calc_uv pipeline
    // -----------------------------------------------------------------------------------------------------------------------
    computeModule = VkUtil::create_shader_module(context, VkUtil::read_byte_file(string_uv));

    bindings.clear();

    binding.binding = 0;								//velocities
    binding.descriptorCount = 1;
    binding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
    bindings.push_back(binding);

    binding.binding = 1;                                //fg
    binding.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
    bindings.push_back(binding);

    binding.binding = 2;                                //pressure
    bindings.push_back(binding);

    binding.binding = 3;                                //cell_flags
    bindings.push_back(binding);

    binding.binding = 4;								//sim_context
    binding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    bindings.push_back(binding);

    VkUtil::create_descriptor_set_layout(context, bindings, &calc_uv.descriptorSetLayout);
    layouts = { calc_uv.descriptorSetLayout };

    VkUtil::create_compute_pipeline(context, computeModule, layouts, &calc_uv.pipelineLayout, &calc_uv.pipeline);

    VkUtil::create_descriptor_sets(context, layouts, calc_uv.descriptorSet);
}

GpuPipelines::~GpuPipelines()
{
    if (sampler) vkDestroySampler(vk_context.device, sampler, vk_context.allocator);

    if (calc_dt.pipeline) vkDestroyPipeline(vk_context.device, calc_dt.pipeline, vk_context.allocator);
    if (calc_alpha.pipeline) vkDestroyPipeline(vk_context.device, calc_alpha.pipeline, vk_context.allocator);
    if (boundary_vals.pipeline) vkDestroyPipeline(vk_context.device, boundary_vals.pipeline, vk_context.allocator);
    if (calc_temp.pipeline) vkDestroyPipeline(vk_context.device, calc_temp.pipeline, vk_context.allocator);
    if (calc_species.pipeline) vkDestroyPipeline(vk_context.device, calc_species.pipeline, vk_context.allocator);
    if (calc_fg.pipeline) vkDestroyPipeline(vk_context.device, calc_fg.pipeline, vk_context.allocator);
    if (calc_rs.pipeline) vkDestroyPipeline(vk_context.device, calc_rs.pipeline, vk_context.allocator);
    if (calc_sor.pipeline) vkDestroyPipeline(vk_context.device, calc_sor.pipeline, vk_context.allocator);
    if (reset_sor.pipeline) vkDestroyPipeline(vk_context.device, reset_sor.pipeline, vk_context.allocator);
    if (calc_uv.pipeline) vkDestroyPipeline(vk_context.device, calc_uv.pipeline, vk_context.allocator);

    if (calc_dt.pipelineLayout) vkDestroyPipelineLayout(vk_context.device, calc_dt.pipelineLayout, vk_context.allocator);
    if (calc_alpha.pipelineLayout) vkDestroyPipelineLayout(vk_context.device, calc_alpha.pipelineLayout, vk_context.allocator);
    if (boundary_vals.pipelineLayout) vkDestroyPipelineLayout(vk_context.device, boundary_vals.pipelineLayout, vk_context.allocator);
    if (calc_temp.pipelineLayout) vkDestroyPipelineLayout(vk_context.device, calc_temp.pipelineLayout, vk_context.allocator);
    if (calc_species.pipelineLayout) vkDestroyPipelineLayout(vk_context.device, calc_species.pipelineLayout, vk_context.allocator);
    if (calc_fg.pipelineLayout) vkDestroyPipelineLayout(vk_context.device, calc_fg.pipelineLayout, vk_context.allocator);
    if (calc_rs.pipelineLayout) vkDestroyPipelineLayout(vk_context.device, calc_rs.pipelineLayout, vk_context.allocator);
    if (calc_sor.pipelineLayout) vkDestroyPipelineLayout(vk_context.device, calc_sor.pipelineLayout, vk_context.allocator);
    if (reset_sor.pipelineLayout) vkDestroyPipelineLayout(vk_context.device, reset_sor.pipelineLayout, vk_context.allocator);
    if (calc_uv.pipelineLayout) vkDestroyPipelineLayout(vk_context.device, calc_uv.pipelineLayout, vk_context.allocator);

    if (calc_dt.descriptorSetLayout) vkDestroyDescriptorSetLayout(vk_context.device, calc_dt.descriptorSetLayout, vk_context.allocator);
    if (calc_alpha.descriptorSetLayout) vkDestroyDescriptorSetLayout(vk_context.device, calc_alpha.descriptorSetLayout, vk_context.allocator);
    if (boundary_vals.descriptorSetLayout) vkDestroyDescriptorSetLayout(vk_context.device, boundary_vals.descriptorSetLayout, vk_context.allocator);
    if (calc_temp.descriptorSetLayout) vkDestroyDescriptorSetLayout(vk_context.device, calc_temp.descriptorSetLayout, vk_context.allocator);
    if (calc_species.descriptorSetLayout) vkDestroyDescriptorSetLayout(vk_context.device, calc_species.descriptorSetLayout, vk_context.allocator);
    if (calc_fg.descriptorSetLayout) vkDestroyDescriptorSetLayout(vk_context.device, calc_fg.descriptorSetLayout, vk_context.allocator);
    if (calc_rs.descriptorSetLayout) vkDestroyDescriptorSetLayout(vk_context.device, calc_rs.descriptorSetLayout, vk_context.allocator);
    if (calc_sor.descriptorSetLayout) vkDestroyDescriptorSetLayout(vk_context.device, calc_sor.descriptorSetLayout, vk_context.allocator);
    if (reset_sor.descriptorSetLayout) vkDestroyDescriptorSetLayout(vk_context.device, reset_sor.descriptorSetLayout, vk_context.allocator);
    if (calc_uv.descriptorSetLayout) vkDestroyDescriptorSetLayout(vk_context.device, calc_uv.descriptorSetLayout, vk_context.allocator);
}

void GpuPipelines::init_descriptor_sets_and_commands(GpuGrid& grid, GpuSimulationInfos& simulation_info, SimulationContext& sim_ctx)
{
    VkImage velocity, temperature_front, temperature_back, species_front, species_back, pressure, density, cell_flags, neighbourhood, fg, rs;
    VkImageView velocity_view, temperature_front_view, temperature_back_view, species_front_view, species_back_view, pressure_view, density_view, cell_flags_view, neighbourhood_view, fg_view, rs_view;
    grid.velocity.get_image(velocity, velocity_view);
    grid.pressure.get_image(pressure, pressure_view);
    grid.temperature.get_cur_image(temperature_front, temperature_front_view);
    grid.temperature.get_prev_image(temperature_back, temperature_back_view);
    grid.density.get_image(density, density_view);
    grid.cell_flags.get_image(cell_flags, cell_flags_view);
    grid.neighbourhood.get_image(neighbourhood, neighbourhood_view);
    grid.species_concentrations.get_cur_image(species_front, species_front_view);
    grid.species_concentrations.get_prev_image(species_back, species_back_view);
    grid.fg.get_image(fg, fg_view);
    grid.rs.get_image(rs, rs_view);

    uint32_t binding = 0;
    VkUtil::update_image_descriptor_set(vk_context, sampler, velocity_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_dt.descriptorSet[0]);
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_simulation_context_buffer(), simulation_info.get_simulation_context_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, calc_dt.descriptorSet[0]);

    binding = 0;
    VkUtil::update_image_descriptor_set(vk_context, sampler, velocity_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_alpha.descriptorSet[0]);
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_simulation_context_buffer(), simulation_info.get_simulation_context_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, calc_alpha.descriptorSet[0]);


    binding = 0;
    // set back images
    VkUtil::update_image_descriptor_set(vk_context, NULL, velocity_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, boundary_vals.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, NULL, pressure_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, boundary_vals.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, NULL, temperature_back_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, boundary_vals.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, NULL, species_back_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, boundary_vals.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, NULL, cell_flags_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, boundary_vals.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, NULL, neighbourhood_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, boundary_vals.descriptorSet[0]);
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_simulation_context_buffer(), simulation_info.get_simulation_context_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, boundary_vals.descriptorSet[0]);
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_inserter_buffer(), simulation_info.get_inserter_buffer_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, boundary_vals.descriptorSet[0]);
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_converter_buffer(), simulation_info.get_converter_buffer_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, boundary_vals.descriptorSet[0]);

    binding = 0;
    // set front images
    VkUtil::update_image_descriptor_set(vk_context, NULL, velocity_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, boundary_vals.descriptorSet[1]);
    VkUtil::update_image_descriptor_set(vk_context, NULL, pressure_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, boundary_vals.descriptorSet[1]);
    VkUtil::update_image_descriptor_set(vk_context, NULL, temperature_front_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, boundary_vals.descriptorSet[1]);
    VkUtil::update_image_descriptor_set(vk_context, NULL, species_front_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, boundary_vals.descriptorSet[1]);
    VkUtil::update_image_descriptor_set(vk_context, NULL, cell_flags_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, boundary_vals.descriptorSet[1]);
    VkUtil::update_image_descriptor_set(vk_context, NULL, neighbourhood_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, boundary_vals.descriptorSet[1]);
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_simulation_context_buffer(), simulation_info.get_simulation_context_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, boundary_vals.descriptorSet[1]);
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_inserter_buffer(), simulation_info.get_inserter_buffer_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, boundary_vals.descriptorSet[1]);
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_converter_buffer(), simulation_info.get_converter_buffer_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, boundary_vals.descriptorSet[1]);

    binding = 0;
    // update from back to front
    VkUtil::update_image_descriptor_set(vk_context, sampler, velocity_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_temp.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, sampler, temperature_back_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_temp.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, NULL, temperature_front_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, calc_temp.descriptorSet[0]);
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_simulation_context_buffer(), simulation_info.get_simulation_context_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, calc_temp.descriptorSet[0]);

    binding = 0;
    // update from front to back
    VkUtil::update_image_descriptor_set(vk_context, sampler, velocity_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_temp.descriptorSet[1]);
    VkUtil::update_image_descriptor_set(vk_context, sampler, temperature_front_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_temp.descriptorSet[1]);
    VkUtil::update_image_descriptor_set(vk_context, NULL, temperature_back_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, calc_temp.descriptorSet[1]);
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_simulation_context_buffer(), simulation_info.get_simulation_context_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, calc_temp.descriptorSet[1]);

    binding = 0;
    // update from back to front
    VkUtil::update_image_descriptor_set(vk_context, sampler, velocity_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_species.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, sampler, species_back_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_species.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, NULL, species_front_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, calc_species.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, sampler, temperature_front_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_species.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, NULL, density_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, calc_species.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, sampler, cell_flags_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_species.descriptorSet[0]);
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_simulation_context_buffer(), simulation_info.get_simulation_context_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, calc_species.descriptorSet[0]);
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_converter_buffer(), simulation_info.get_converter_buffer_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, calc_species.descriptorSet[0]);

    binding = 0;
    // update from front to back
    VkUtil::update_image_descriptor_set(vk_context, sampler, velocity_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_species.descriptorSet[1]);
    VkUtil::update_image_descriptor_set(vk_context, sampler, species_front_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_species.descriptorSet[1]);
    VkUtil::update_image_descriptor_set(vk_context, NULL, species_back_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, calc_species.descriptorSet[1]);
    VkUtil::update_image_descriptor_set(vk_context, sampler, temperature_back_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_species.descriptorSet[1]);
    VkUtil::update_image_descriptor_set(vk_context, NULL, density_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, calc_species.descriptorSet[1]);
    VkUtil::update_image_descriptor_set(vk_context, sampler, cell_flags_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_species.descriptorSet[1]);
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_simulation_context_buffer(), simulation_info.get_simulation_context_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, calc_species.descriptorSet[1]);
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_converter_buffer(), simulation_info.get_converter_buffer_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, calc_species.descriptorSet[1]);

    binding = 0;
    VkUtil::update_image_descriptor_set(vk_context, NULL, velocity_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, calc_fg.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, NULL, fg_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, calc_fg.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, sampler, density_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_fg.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, sampler, cell_flags_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_fg.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, sampler, neighbourhood_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_fg.descriptorSet[0]);
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_simulation_context_buffer(), simulation_info.get_simulation_context_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, calc_fg.descriptorSet[0]);

    binding = 0;
    VkUtil::update_image_descriptor_set(vk_context, NULL, rs_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, calc_rs.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, sampler, fg_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_rs.descriptorSet[0]);
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_simulation_context_buffer(), simulation_info.get_simulation_context_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, calc_rs.descriptorSet[0]);

    binding = 0;
    VkUtil::update_image_descriptor_set(vk_context, sampler, rs_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_sor.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, NULL, pressure_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, calc_sor.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, sampler, cell_flags_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_sor.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, sampler, neighbourhood_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_sor.descriptorSet[0]);
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_simulation_context_buffer(), simulation_info.get_simulation_context_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, calc_sor.descriptorSet[0]);

    binding = 0;
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_simulation_context_buffer(), simulation_info.get_simulation_context_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, reset_sor.descriptorSet[0]);

    binding = 0;
    VkUtil::update_image_descriptor_set(vk_context, NULL, velocity_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, binding++, calc_uv.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, sampler, fg_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_uv.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, sampler, pressure_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_uv.descriptorSet[0]);
    VkUtil::update_image_descriptor_set(vk_context, sampler, cell_flags_view, VK_IMAGE_LAYOUT_GENERAL, VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, binding++, calc_uv.descriptorSet[0]);
    VkUtil::update_descriptor_set(vk_context, simulation_info.get_simulation_context_buffer(), simulation_info.get_simulation_context_byte_size(), binding++, VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, calc_uv.descriptorSet[0]);

    //creating the command buffer. One for updating the front buffers and one for updating the back buffer
    uint32_t dispatch_x = (grid._imax + 31) / 32, dispatch_y = (grid._jmax + 31) / 32;
    uint32_t dispatch_x_half = (dispatch_x + 1) / 2;
    VkUtil::create_command_buffer(vk_context, simulate_commands);
    VkUtil::create_command_buffer(vk_context, simulate_commands + 1);

    vkCmdBindPipeline(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, calc_dt.pipeline);
    vkCmdBindPipeline(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, calc_dt.pipeline);
    vkCmdBindDescriptorSets(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, calc_dt.pipelineLayout, 0, 1, calc_dt.descriptorSet, 0, nullptr);
    vkCmdBindDescriptorSets(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, calc_dt.pipelineLayout, 0, 1, calc_dt.descriptorSet, 0, nullptr);
    vkCmdDispatch(simulate_commands[0], dispatch_x, dispatch_y, 1);
    vkCmdDispatch(simulate_commands[1], dispatch_x, dispatch_y, 1);

    vkCmdPipelineBarrier(simulate_commands[0], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 0, 0, 0, 0);
    vkCmdPipelineBarrier(simulate_commands[1], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 0, 0, 0, 0);

    vkCmdBindPipeline(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, calc_alpha.pipeline);
    vkCmdBindPipeline(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, calc_alpha.pipeline);
    vkCmdBindDescriptorSets(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, calc_alpha.pipelineLayout, 0, 1, calc_alpha.descriptorSet, 0, nullptr);
    vkCmdBindDescriptorSets(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, calc_alpha.pipelineLayout, 0, 1, calc_alpha.descriptorSet, 0, nullptr);
    vkCmdDispatch(simulate_commands[0], dispatch_x, dispatch_y, 1);
    vkCmdDispatch(simulate_commands[1], dispatch_x, dispatch_y, 1);

    vkCmdPipelineBarrier(simulate_commands[0], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 0, 0, 0, 0);
    vkCmdPipelineBarrier(simulate_commands[1], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 0, 0, 0, 0);

    vkCmdBindPipeline(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, boundary_vals.pipeline);
    vkCmdBindPipeline(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, boundary_vals.pipeline);
    vkCmdBindDescriptorSets(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, boundary_vals.pipelineLayout, 0, 1, boundary_vals.descriptorSet, 0, nullptr);
    vkCmdBindDescriptorSets(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, boundary_vals.pipelineLayout, 0, 1, boundary_vals.descriptorSet + 1, 0, nullptr);
    vkCmdDispatch(simulate_commands[0], dispatch_x, dispatch_y, 1);
    vkCmdDispatch(simulate_commands[1], dispatch_x, dispatch_y, 1);

    vkCmdPipelineBarrier(simulate_commands[0], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 0, 0, 0, 0);
    vkCmdPipelineBarrier(simulate_commands[1], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 0, 0, 0, 0);

    vkCmdBindPipeline(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, calc_temp.pipeline);
    vkCmdBindPipeline(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, calc_temp.pipeline);
    vkCmdBindDescriptorSets(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, calc_temp.pipelineLayout, 0, 1, calc_temp.descriptorSet, 0, nullptr);
    vkCmdBindDescriptorSets(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, calc_temp.pipelineLayout, 0, 1, calc_temp.descriptorSet + 1, 0, nullptr);
    vkCmdDispatch(simulate_commands[0], dispatch_x, dispatch_y, 1);
    vkCmdDispatch(simulate_commands[1], dispatch_x, dispatch_y, 1);

    vkCmdPipelineBarrier(simulate_commands[0], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 0, 0, 0, 0);
    vkCmdPipelineBarrier(simulate_commands[1], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 0, 0, 0, 0);

    vkCmdBindPipeline(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, calc_species.pipeline);
    vkCmdBindPipeline(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, calc_species.pipeline);
    vkCmdBindDescriptorSets(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, calc_species.pipelineLayout, 0, 1, calc_species.descriptorSet, 0, nullptr);
    vkCmdBindDescriptorSets(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, calc_species.pipelineLayout, 0, 1, &calc_species.descriptorSet[1], 0, nullptr);
    vkCmdDispatch(simulate_commands[0], dispatch_x, dispatch_y, 1);
    vkCmdDispatch(simulate_commands[1], dispatch_x, dispatch_y, 1);

    vkCmdPipelineBarrier(simulate_commands[0], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 0, 0, 0, 0);
    vkCmdPipelineBarrier(simulate_commands[1], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 0, 0, 0, 0);

    vkCmdBindPipeline(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, calc_fg.pipeline);
    vkCmdBindPipeline(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, calc_fg.pipeline);
    vkCmdBindDescriptorSets(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, calc_fg.pipelineLayout, 0, 1, calc_fg.descriptorSet, 0, nullptr);
    vkCmdBindDescriptorSets(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, calc_fg.pipelineLayout, 0, 1, calc_fg.descriptorSet, 0, nullptr);
    vkCmdDispatch(simulate_commands[0], dispatch_x, dispatch_y, 1);
    vkCmdDispatch(simulate_commands[1], dispatch_x, dispatch_y, 1);

    vkCmdPipelineBarrier(simulate_commands[0], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 0, 0, 0, 0);
    vkCmdPipelineBarrier(simulate_commands[1], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 0, 0, 0, 0);

    vkCmdBindPipeline(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, calc_rs.pipeline);
    vkCmdBindPipeline(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, calc_rs.pipeline);
    vkCmdBindDescriptorSets(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, calc_rs.pipelineLayout, 0, 1, calc_rs.descriptorSet, 0, nullptr);
    vkCmdBindDescriptorSets(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, calc_rs.pipelineLayout, 0, 1, calc_rs.descriptorSet, 0, nullptr);
    vkCmdDispatch(simulate_commands[0], dispatch_x, dispatch_y, 1);
    vkCmdDispatch(simulate_commands[1], dispatch_x, dispatch_y, 1);

    vkCmdPipelineBarrier(simulate_commands[0], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 0, 0, 0, 0);
    vkCmdPipelineBarrier(simulate_commands[1], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 0, 0, 0, 0);

    uint32_t push_constant0 = 0, push_constant1 = 1;
    for (int i = 0; i < sim_ctx.itermax; ++i) {
        vkCmdBindPipeline(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, calc_sor.pipeline);
        vkCmdBindPipeline(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, calc_sor.pipeline);
        vkCmdBindDescriptorSets(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, calc_sor.pipelineLayout, 0, 1, calc_sor.descriptorSet, 0, nullptr);
        vkCmdBindDescriptorSets(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, calc_sor.pipelineLayout, 0, 1, calc_sor.descriptorSet, 0, nullptr);
        vkCmdPushConstants(simulate_commands[0], calc_sor.pipelineLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(uint32_t), &push_constant0);
        vkCmdPushConstants(simulate_commands[1], calc_sor.pipelineLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(uint32_t), &push_constant0);
        vkCmdDispatchIndirect(simulate_commands[0], simulation_info.get_simulation_context_buffer(), offsetof(GpuSimulationContext, patch_x));
        vkCmdDispatchIndirect(simulate_commands[1], simulation_info.get_simulation_context_buffer(), offsetof(GpuSimulationContext, patch_x));
    
        vkCmdPipelineBarrier(simulate_commands[0], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 0, 0, 0, 0);
        vkCmdPipelineBarrier(simulate_commands[1], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 0, 0, 0, 0);
    
        vkCmdPushConstants(simulate_commands[0], calc_sor.pipelineLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(uint32_t), &push_constant1);
        vkCmdPushConstants(simulate_commands[1], calc_sor.pipelineLayout, VK_SHADER_STAGE_COMPUTE_BIT, 0, sizeof(uint32_t), &push_constant1);
        vkCmdDispatchIndirect(simulate_commands[0], simulation_info.get_simulation_context_buffer(), offsetof(GpuSimulationContext, patch_x));
        vkCmdDispatchIndirect(simulate_commands[1], simulation_info.get_simulation_context_buffer(), offsetof(GpuSimulationContext, patch_x));

        VkMemoryBarrier memBarrier{};
        memBarrier.sType = VK_STRUCTURE_TYPE_MEMORY_BARRIER;
        memBarrier.srcAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;
        memBarrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;

        vkCmdPipelineBarrier(simulate_commands[0], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memBarrier, 0, 0, 0, 0);
        vkCmdPipelineBarrier(simulate_commands[1], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memBarrier, 0, 0, 0, 0);
    
        vkCmdBindPipeline(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, reset_sor.pipeline);
        vkCmdBindPipeline(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, reset_sor.pipeline);
        vkCmdBindDescriptorSets(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, reset_sor.pipelineLayout, 0, 1, reset_sor.descriptorSet, 0, nullptr);
        vkCmdBindDescriptorSets(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, reset_sor.pipelineLayout, 0, 1, reset_sor.descriptorSet, 0, nullptr);
        vkCmdDispatch(simulate_commands[0], 1, 1, 1);
        vkCmdDispatch(simulate_commands[1], 1, 1, 1);

        vkCmdPipelineBarrier(simulate_commands[0], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memBarrier, 0, 0, 0, 0);
        vkCmdPipelineBarrier(simulate_commands[1], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 1, &memBarrier, 0, 0, 0, 0);
    }

    vkCmdBindPipeline(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, calc_uv.pipeline);
    vkCmdBindPipeline(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, calc_uv.pipeline);
    vkCmdBindDescriptorSets(simulate_commands[0], VK_PIPELINE_BIND_POINT_COMPUTE, calc_uv.pipelineLayout, 0, 1, calc_uv.descriptorSet, 0, nullptr);
    vkCmdBindDescriptorSets(simulate_commands[1], VK_PIPELINE_BIND_POINT_COMPUTE, calc_uv.pipelineLayout, 0, 1, calc_uv.descriptorSet, 0, nullptr);
    vkCmdDispatch(simulate_commands[0], dispatch_x, dispatch_y, 1);
    vkCmdDispatch(simulate_commands[1], dispatch_x, dispatch_y, 1);

    vkCmdPipelineBarrier(simulate_commands[0], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 0, 0, 0, 0);
    vkCmdPipelineBarrier(simulate_commands[1], VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT, 0, 0, 0, 0, 0, 0, 0);

    vkEndCommandBuffer(simulate_commands[0]);
    vkEndCommandBuffer(simulate_commands[1]);
}