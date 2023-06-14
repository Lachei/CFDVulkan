#pragma once
#include <vulkan/vulkan.h>
#include <vector>
#include <string>
#include <cassert>
#include <fstream>
#include <iostream>
#include <cstring>

#ifdef _DEBUG
#define VULKAN_DEBUG_REPORT
#endif

namespace VkUtil {
	struct VkContext {
		VkQueueFlagBits          queue_type = VK_QUEUE_FLAG_BITS_MAX_ENUM;
		uint32_t                 queue_family = (uint32_t)-1;
		VkInstance               instance;
		VkAllocationCallbacks*   allocator = nullptr;
		VkDebugReportCallbackEXT debug_report;
		VkPhysicalDevice         physical_device;
		VkDevice                 device;
		VkQueue                  queue;
		VkDescriptorPool         descriptor_pool;
        VkCommandPool            command_pool;
	};

    uint32_t find_memory_type(VkContext& context, uint32_t typeFilter, VkMemoryPropertyFlags properties);

	void check_vk_result(VkResult err);
	void create_vk_context(std::vector<const char*>& extensions, uint32_t gpu_index, VkContext& context);
    void destroy_vk_context(VkContext& context);
    void create_compute_pipeline(VkContext& context, VkShaderModule& shaderModule, std::vector<VkDescriptorSetLayout> descriptorLayouts, VkPipelineLayout* pipelineLayout, VkPipeline* pipeline);
    void create_compute_pipeline(VkContext& context, VkShaderModule& shaderModule, std::vector<VkDescriptorSetLayout> descriptorLayouts, std::vector<VkPushConstantRange>& push_constants, VkPipelineLayout* pipelineLayout, VkPipeline* pipeline);
    void create_descriptor_set_layout(VkContext& context, const std::vector<VkDescriptorSetLayoutBinding>& bindings, VkDescriptorSetLayout* descriptorSetLayout);
    void create_descriptor_setLayout_partially_bound(VkContext& context, const std::vector<VkDescriptorSetLayoutBinding>& bindings, const std::vector<bool>& enableValidation, VkDescriptorSetLayout* descriptorSetLayout);
    void create_descriptor_sets(VkContext& context, const std::vector<VkDescriptorSetLayout>& layouts, VkDescriptorSet* descriptorSetArray);
    void update_descriptor_set(VkContext& context, VkBuffer buffer, uint32_t size, uint32_t binding, VkDescriptorType descriptorType, VkDescriptorSet descriptorSet);
    void update_image_descriptor_set(VkContext& context, VkSampler sampler, VkImageView imageView, VkImageLayout imageLayout, VkDescriptorType imageType, uint32_t binding, VkDescriptorSet descriptorSet);
    void create_image(VkContext& context, uint32_t width, uint32_t height, VkFormat imageFormat, VkImageUsageFlags usageFlags, VkImage* image);
    void create_image_view(VkContext& context, VkImage image, VkFormat imageFormat, uint32_t mipLevelCount, VkImageAspectFlags aspectMask, VkImageView* imageView);
    void create_image_array(VkContext& context, uint32_t width, uint32_t height, uint32_t array_count, VkFormat imageFormat, VkImageUsageFlags usageFlags, VkImage* image);
    void create_image_view_array(VkContext& context, VkImage image, VkFormat imageFormat, uint32_t mipLevelCount, uint32_t array_count, VkImageAspectFlags aspectMask, VkImageView* imageView);
    void create_buffer(VkContext& context, uint32_t size, VkBufferUsageFlags usage, VkBuffer* buffer);
    void create_image_sampler(VkContext& context, VkSamplerAddressMode adressMode, VkFilter filter, uint16_t maxAnisotropy, uint16_t mipLevels, VkSampler* sampler);
    void create_command_buffer(VkContext& context, VkCommandBuffer* commandBuffer);
    void commit_command_buffer(VkContext& context, VkCommandBuffer commandBuffer);
    void execute_command_buffer(VkContext& context, VkCommandBuffer commandBuffer);

    void copy_buffer_to_image(VkCommandBuffer commandBuffer, VkBuffer buffer, VkImageLayout layout, VkImage image, uint32_t width, uint32_t height, uint32_t depth, uint32_t array_layer);
    void copy_image_to_buffer(VkCommandBuffer commandBuffer, VkBuffer buffer, VkImageLayout layout, VkImage image, uint32_t width, uint32_t height, uint32_t depth, uint32_t array_layer);
    void upload_data(VkContext& context, VkDeviceMemory memory, uint32_t offset, uint32_t byteSize, void* data);
    void download_data(VkContext& context, VkDeviceMemory memory, uint32_t offset, uint32_t byteSize, void* data);
    void upload_image_data(VkContext& context, VkImage image, VkImageLayout imageLayout, VkFormat imageFormat, uint32_t x, uint32_t y, uint32_t z, uint32_t array_layer, void* data, uint32_t byteSize);
    void download_image_data(VkContext& context, VkImage image, VkImageLayout layout, VkFormat imageFormat, uint32_t x, uint32_t y, uint32_t z, uint32_t array_layer, void* data, uint32_t byteSize);

    void transition_image_layout(VkCommandBuffer commandBuffer, VkImage image, VkFormat format, uint32_t layer_count, VkImageLayout oldLayout, VkImageLayout newLayout);
    void pipeline_barrier(VkCommandBuffer commandBuffer, VkPipelineStageFlags stage_flags);

    VkShaderModule create_shader_module(VkContext& context, const std::vector<char>& byteArr);
    std::vector<char> read_byte_file(const std::string& filename);
};