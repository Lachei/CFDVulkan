#pragma once
#include "vk_utils.h"
#include "datastructures.hpp"
#include <cassert>

class GpuDoubleImage {
public:
    GpuDoubleImage();
    GpuDoubleImage(VkUtil::VkContext& vk_context, int width, int height, VkFormat image_format, VkImageUsageFlags image_usage, VkMemoryPropertyFlags memory_properties);
    ~GpuDoubleImage();

    // swaps current and previous image
    void swap_images();
    void get_cur_image(VkImage& image, VkImageView& imageView);
    void get_prev_image(VkImage& image, VkImageView& imageView);

    void set_cur_image(void* data, uint32_t byte_size);
    void get_cur_image(void* data, uint32_t byte_size);
private:
    VkUtil::VkContext vk_context;
    VkFormat format_image;
    uint32_t offset_images[2];
    VkImage images[2];
    VkImageView imageViews[2];
    VkDeviceMemory memory_images;

    uint32_t current_image_index;
    uint32_t width, height;
};

class GpuDoubleArrayImage {
public:
    GpuDoubleArrayImage();
    GpuDoubleArrayImage(VkUtil::VkContext& vk_context, int width, int height, int amt_of_images, VkFormat image_format, VkImageUsageFlags image_usage, VkMemoryPropertyFlags memory_properties);
    ~GpuDoubleArrayImage();

    // swaps current and previous image
    void swap_images();
    void get_cur_image(VkImage& image, VkImageView& imageView);
    void get_prev_image(VkImage& image, VkImageView& imageView);

    void set_cur_image(uint32_t array_layer, void* data, uint32_t byte_size);
    void get_cur_image(uint32_t array_layer, void* data, uint32_t byte_size);
private:
    VkUtil::VkContext vk_context;
    VkFormat format_image;
    uint32_t offset_images[2];
    VkImage images[2];
    VkImageView imageViews[2];
    VkDeviceMemory memory_images;

    uint32_t current_image_index;
    uint32_t width, height, amt_of_images;
};

class GpuImage {
public:
    GpuImage();
    GpuImage(VkUtil::VkContext& vk_context, int width, int height, VkFormat image_format, VkImageUsageFlags image_usage, VkMemoryPropertyFlags memory_properties);
    ~GpuImage();

    void get_image(VkImage& image, VkImageView& imageView);

    void set_image(void* data, uint32_t byte_size);
    void get_image(void* data, uint32_t byte_size);
private:
    VkUtil::VkContext vk_context;
    VkFormat format_image;
    VkImage image;
    VkImageView imageView;
    VkDeviceMemory memory_image;

    uint32_t width, height;
};

class GpuBuffer {
public:
    GpuBuffer();
    GpuBuffer(VkUtil::VkContext& vk_context, uint32_t byte_size, VkBufferUsageFlags buffer_usage);
    ~GpuBuffer();

    VkBuffer buffer;

    void set_data(void* data, uint32_t byte_size);
    void get_data(void* data, uint32_t byte_size);
private:
    VkUtil::VkContext vk_context;
    VkDeviceMemory memory;
    uint32_t byte_size;
};

class GpuSimulationInfos {
public:
    GpuSimulationInfos();
    GpuSimulationInfos(VkUtil::VkContext& vk_context, SimulationContext& sim_context);
    ~GpuSimulationInfos();

    void get_simulation_context_data(void* data, uint32_t bytesize);
    void set_simulation_context_data(void* data, uint32_t bytesize);
    VkBuffer get_simulation_context_buffer();
    uint32_t get_simulation_context_byte_size();
    VkBuffer get_converter_buffer();
    uint32_t get_converter_buffer_byte_size();
    VkBuffer get_inserter_buffer();
    uint32_t get_inserter_buffer_byte_size();
private:
    uint32_t sim_ctx_byte_size, converter_byte_size, inserter_byte_size;
    GpuBuffer sim_ctx, converter, inserter;
};