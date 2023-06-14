#include "gpu_images.h"

GpuDoubleImage::GpuDoubleImage()
{
}

GpuDoubleImage::GpuDoubleImage(VkUtil::VkContext& vk_context, int width, int height, VkFormat image_format, VkImageUsageFlags image_usage, VkMemoryPropertyFlags memory_properties):
    vk_context(vk_context), format_image(image_format), width(width), height(height), current_image_index(0)
{
    /* creating the images*/
    VkUtil::create_image(vk_context, width, height, image_format, image_usage, images);
    VkUtil::create_image(vk_context, width, height, image_format, image_usage, images + 1);

    VkMemoryAllocateInfo memAlloc{ VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO };
    VkMemoryRequirements memReq;
    offset_images[0] = 0;
    vkGetImageMemoryRequirements(vk_context.device, images[0], &memReq);
    memAlloc.allocationSize = memReq.size;
    offset_images[1] = memAlloc.allocationSize;
    vkGetImageMemoryRequirements(vk_context.device, images[1], &memReq);
    memAlloc.allocationSize += memReq.size;
    memAlloc.memoryTypeIndex = VkUtil::find_memory_type(vk_context, memReq.memoryTypeBits, memory_properties);
    VkResult err = vkAllocateMemory(vk_context.device, &memAlloc, vk_context.allocator, &memory_images);    VkUtil::check_vk_result(err);

    err = vkBindImageMemory(vk_context.device, images[0], memory_images, offset_images[0]);                 VkUtil::check_vk_result(err);
    err = vkBindImageMemory(vk_context.device, images[1], memory_images, offset_images[1]);                 VkUtil::check_vk_result(err);

    VkUtil::create_image_view(vk_context, images[0], image_format, 1, VK_IMAGE_ASPECT_COLOR_BIT, imageViews);
    VkUtil::create_image_view(vk_context, images[1], image_format, 1, VK_IMAGE_ASPECT_COLOR_BIT, imageViews + 1);
    
    /* transform image layout to general layout*/
    VkCommandBuffer commands;
    VkUtil::create_command_buffer(vk_context, &commands);
    VkUtil::transition_image_layout(commands, images[0], image_format, 1, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_GENERAL);
    VkUtil::transition_image_layout(commands, images[1], image_format, 1, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_GENERAL);
    VkUtil::commit_command_buffer(vk_context, commands);
    vkQueueWaitIdle(vk_context.queue);
    vkFreeCommandBuffers(vk_context.device, vk_context.command_pool, 1, &commands);
}

GpuDoubleImage::~GpuDoubleImage()
{
    if (images[0]) vkDestroyImage(vk_context.device, images[0], vk_context.allocator);
    if (images[1]) vkDestroyImage(vk_context.device, images[1], vk_context.allocator);
    if (imageViews[0]) vkDestroyImageView(vk_context.device, imageViews[0], vk_context.allocator);
    if (imageViews[1]) vkDestroyImageView(vk_context.device, imageViews[1], vk_context.allocator);
    if (memory_images) vkFreeMemory(vk_context.device, memory_images, vk_context.allocator);
}

void GpuDoubleImage::swap_images()
{
    current_image_index ^= 1;
}

void GpuDoubleImage::get_cur_image(VkImage& image, VkImageView& imageView)
{
    image = images[current_image_index];
    imageView = imageViews[current_image_index];
}

void GpuDoubleImage::get_prev_image(VkImage& image, VkImageView& imageView)
{
    image = images[current_image_index ^ 1];
    imageView = imageViews[current_image_index ^ 1];
}

void GpuDoubleImage::set_cur_image(void* data, uint32_t byte_size)
{
    VkUtil::upload_image_data(vk_context, images[current_image_index], VK_IMAGE_LAYOUT_GENERAL, format_image, width, height, 1, 0, data, byte_size);
}

void GpuDoubleImage::get_cur_image(void* data, uint32_t byte_size)
{
    VkUtil::download_image_data(vk_context, images[current_image_index], VK_IMAGE_LAYOUT_GENERAL, format_image, width, height, 1, 0, data, byte_size);
}

GpuDoubleArrayImage::GpuDoubleArrayImage()
{
}

GpuDoubleArrayImage::GpuDoubleArrayImage(VkUtil::VkContext& vk_context, int width, int height, int amt_of_images, VkFormat image_format, VkImageUsageFlags image_usage, VkMemoryPropertyFlags memory_properties):
    vk_context(vk_context), format_image(image_format), width(width), height(height), amt_of_images(amt_of_images), current_image_index(0)
{
    /* creating the images*/
    VkUtil::create_image_array(vk_context, width, height, amt_of_images, image_format, image_usage, images);
    VkUtil::create_image_array(vk_context, width, height, amt_of_images, image_format, image_usage, images + 1);

    VkMemoryAllocateInfo memAlloc{ VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO };
    VkMemoryRequirements memReq;
    offset_images[0] = 0;
    vkGetImageMemoryRequirements(vk_context.device, images[0], &memReq);
    memAlloc.allocationSize = memReq.size;
    offset_images[1] = memAlloc.allocationSize;
    vkGetImageMemoryRequirements(vk_context.device, images[1], &memReq);
    memAlloc.allocationSize += memReq.size;
    memAlloc.memoryTypeIndex = VkUtil::find_memory_type(vk_context, memReq.memoryTypeBits, memory_properties);
    VkResult err = vkAllocateMemory(vk_context.device, &memAlloc, vk_context.allocator, &memory_images);    VkUtil::check_vk_result(err);

    err = vkBindImageMemory(vk_context.device, images[0], memory_images, offset_images[0]);                 VkUtil::check_vk_result(err);
    err = vkBindImageMemory(vk_context.device, images[1], memory_images, offset_images[1]);                 VkUtil::check_vk_result(err);

    VkUtil::create_image_view_array(vk_context, images[0], image_format, 1, amt_of_images, VK_IMAGE_ASPECT_COLOR_BIT, imageViews);
    VkUtil::create_image_view_array(vk_context, images[1], image_format, 1, amt_of_images, VK_IMAGE_ASPECT_COLOR_BIT, imageViews + 1);

    /* transform image layout to general layout*/
    VkCommandBuffer commands;
    VkUtil::create_command_buffer(vk_context, &commands);
    VkUtil::transition_image_layout(commands, images[0], image_format, amt_of_images, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_GENERAL);
    VkUtil::transition_image_layout(commands, images[1], image_format, amt_of_images, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_GENERAL);
    VkUtil::commit_command_buffer(vk_context, commands);
    vkQueueWaitIdle(vk_context.queue);
    vkFreeCommandBuffers(vk_context.device, vk_context.command_pool, 1, &commands);
}

GpuDoubleArrayImage::~GpuDoubleArrayImage()
{
    if (images[0]) vkDestroyImage(vk_context.device, images[0], vk_context.allocator);
    if (images[1]) vkDestroyImage(vk_context.device, images[1], vk_context.allocator);
    if (imageViews[0]) vkDestroyImageView(vk_context.device, imageViews[0], vk_context.allocator);
    if (imageViews[1]) vkDestroyImageView(vk_context.device, imageViews[1], vk_context.allocator);
    if (memory_images) vkFreeMemory(vk_context.device, memory_images, vk_context.allocator);
}

void GpuDoubleArrayImage::swap_images()
{
    current_image_index ^= 1;
}

void GpuDoubleArrayImage::get_cur_image(VkImage& image, VkImageView& imageView)
{
    image = images[current_image_index];
    imageView = imageViews[current_image_index];
}

void GpuDoubleArrayImage::get_prev_image(VkImage& image, VkImageView& imageView)
{
    image = images[current_image_index ^ 1];
    imageView = imageViews[current_image_index ^ 1];
}

void GpuDoubleArrayImage::set_cur_image(uint32_t array_layer, void* data, uint32_t byte_size)
{
    assert(array_layer < amt_of_images);
    VkUtil::upload_image_data(vk_context, images[current_image_index], VK_IMAGE_LAYOUT_GENERAL, format_image, width, height, 1, array_layer, data, byte_size);
}

void GpuDoubleArrayImage::get_cur_image(uint32_t array_layer, void* data, uint32_t byte_size)
{
    VkUtil::download_image_data(vk_context, images[current_image_index], VK_IMAGE_LAYOUT_GENERAL, format_image, width, height, 1, array_layer, data, byte_size);
}

GpuImage::GpuImage():
    image(0), imageView(0), memory_image(0)
{
}

GpuImage::GpuImage(VkUtil::VkContext& vk_context, int width, int height, VkFormat image_format, VkImageUsageFlags image_usage, VkMemoryPropertyFlags memory_properties) :
    vk_context(vk_context), width(width), height(height), format_image(image_format)
{
    /* creating the images*/
    VkUtil::create_image(vk_context, width, height, image_format, image_usage, &image);

    VkMemoryAllocateInfo memAlloc{ VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO };
    VkMemoryRequirements memReq;
    vkGetImageMemoryRequirements(vk_context.device, image, &memReq);
    memAlloc.allocationSize = memReq.size;
    memAlloc.memoryTypeIndex = VkUtil::find_memory_type(vk_context, memReq.memoryTypeBits, memory_properties);
    VkResult err = vkAllocateMemory(vk_context.device, &memAlloc, vk_context.allocator, &memory_image);     VkUtil::check_vk_result(err);

    err = vkBindImageMemory(vk_context.device, image, memory_image, 0);                                     VkUtil::check_vk_result(err);

    VkUtil::create_image_view(vk_context, image, image_format, 1, VK_IMAGE_ASPECT_COLOR_BIT, &imageView);

    /* transform image layout to general layout*/
    VkCommandBuffer commands;
    VkUtil::create_command_buffer(vk_context, &commands);
    VkUtil::transition_image_layout(commands, image, image_format, 1, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_GENERAL);
    VkUtil::commit_command_buffer(vk_context, commands);
    vkQueueWaitIdle(vk_context.queue);
    vkFreeCommandBuffers(vk_context.device, vk_context.command_pool, 1, &commands);
}

GpuImage::~GpuImage()
{
    if (image) vkDestroyImage(vk_context.device, image, vk_context.allocator);
    if (imageView) vkDestroyImageView(vk_context.device, imageView, vk_context.allocator);
    if (memory_image) vkFreeMemory(vk_context.device, memory_image, vk_context.allocator);
}

void GpuImage::get_image(VkImage& image, VkImageView& imageView)
{
    image = this->image;
    imageView = this->imageView;
}

void GpuImage::set_image(void* data, uint32_t byte_size)
{
    VkUtil::upload_image_data(vk_context, image, VK_IMAGE_LAYOUT_GENERAL, format_image, width, height, 1, 0, data, byte_size);
}

void GpuImage::get_image(void* data, uint32_t byte_size)
{
    VkUtil::download_image_data(vk_context, image, VK_IMAGE_LAYOUT_GENERAL, format_image, width, height, 1, 0, data, byte_size);
}

GpuBuffer::GpuBuffer()
{
}

GpuBuffer::GpuBuffer(VkUtil::VkContext& vk_context, uint32_t byte_size, VkBufferUsageFlags buffer_usage):
    byte_size(byte_size), vk_context(vk_context)
{
    VkUtil::create_buffer(vk_context, byte_size, buffer_usage, &buffer);
    VkMemoryAllocateInfo allocInfo{ VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO };
    VkMemoryRequirements memReq;
    vkGetBufferMemoryRequirements(vk_context.device, buffer, &memReq);
    allocInfo.allocationSize = memReq.size;
    allocInfo.memoryTypeIndex = VkUtil::find_memory_type(vk_context, memReq.memoryTypeBits, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);
    VkResult r = vkAllocateMemory(vk_context.device, &allocInfo, vk_context.allocator, &memory);    VkUtil::check_vk_result(r);
    r = vkBindBufferMemory(vk_context.device, buffer, memory, 0);                                   VkUtil::check_vk_result(r);
}

GpuBuffer::~GpuBuffer()
{
    if (buffer) vkDestroyBuffer(vk_context.device, buffer, vk_context.allocator);
    if (memory) vkFreeMemory(vk_context.device, memory, vk_context.allocator);
}

void GpuBuffer::set_data(void* data, uint32_t byte_size)
{
    assert(byte_size <= this->byte_size);
    VkUtil::upload_data(vk_context, memory, 0, byte_size, data);
}

void GpuBuffer::get_data(void* data, uint32_t byte_size)
{
    assert(byte_size <= this->byte_size);
    VkUtil::download_data(vk_context, memory, 0, byte_size, data);
}

GpuSimulationInfos::GpuSimulationInfos()
{
}

GpuSimulationInfos::GpuSimulationInfos(VkUtil::VkContext& vk_context, SimulationContext& sim_context):
    sim_ctx_byte_size(sizeof(GpuSimulationContext) + sim_context.species.size() * sizeof(float) * 2),
    converter_byte_size((sim_context.species_converters.size() * (sim_context.species.size() + sim_context.species.size() * sim_context.species.size()) + 1) * sizeof(float) + sizeof(int)),
    inserter_byte_size(sim_context.multi_inserters.size() * (sim_context.species.size() + 3) * sizeof(float) + sizeof(int)),
    sim_ctx(vk_context, sim_ctx_byte_size, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT | VK_BUFFER_USAGE_INDIRECT_BUFFER_BIT),
    converter(vk_context, converter_byte_size, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT),
    inserter(vk_context, converter_byte_size, VK_BUFFER_USAGE_STORAGE_BUFFER_BIT)
{
    float* tmp = new float[sim_ctx_byte_size / sizeof(float)];
    GpuSimulationContext* g = (GpuSimulationContext*)tmp;
    uint32_t dispatch_x = (sim_context.imax + 2 + 31) / 32, dispatch_y = (sim_context.jmax + 2 + 31) / 32;
    uint32_t dispatch_x_half = (dispatch_x + 1) / 2;
    g->Re = sim_context.Re;
    g->UI = sim_context.UI;
    g->VI = sim_context.VI;
    g->PI = sim_context.PI;
    g->TI = sim_context.TI;
    g->UIn = sim_context.UIn;
    g->VIn = sim_context.VIn;
    g->WTN = sim_context.WTI.N;
    g->WTE = sim_context.WTI.E;
    g->WTS = sim_context.WTI.S;
    g->WTW = sim_context.WTI.W;
    g->GX = sim_context.GX;
    g->GY = sim_context.GY;
    g->t_end = sim_context.t_end;
    g->dt = sim_context.dt;
    g->dx = sim_context.dx;
    g->dy = sim_context.dy;
    g->xlenght = sim_context.xlength;
    g->ylength = sim_context.ylength;
    g->imax = sim_context.imax;
    g->jmax = sim_context.jmax;
    g->alpha = 0;
    g->omg = sim_context.omg;
    g->tau = sim_context.tau;
    g->itermax = sim_context.itermax;
    g->eps = 1e-3;
    g->dt_value = sim_context.dt_value;
    g->Pr = sim_context.Pr;
    g->beta = sim_context.beta;
    g->sor_eps = sim_context.eps;
    g->cur_sor_eps = 0;
    g->prev_sor_eps = sim_context.eps + 1;
    g->divider = 0;
    g->lock = 0;
    g->sor_counter = 0;
    g->new_dt = *((uint32_t*)(&g->t_end));
    g->patch_x = dispatch_x_half;
    g->patch_x_origin = dispatch_x_half;
    g->patch_y = dispatch_y;
    g->patch_y_origin = dispatch_y;
    g->patch_z = 1;
    g->patch_z_origin = 1;
    g->model = TEST_CASE;
    g->o2index = 0;
    for (; g->o2index < sim_context.species.size() && sim_context.species[g->o2index].name != "O2"; g->o2index++) {}
    g->amt_species = sim_context.species.size();
    float* sp = (float*)(((char*)tmp) + sizeof(GpuSimulationContext));
    uint32_t counter = 0;
    for (auto& spe : sim_context.species) {
        sp[counter++] = spe.molar_mass;
        sp[counter++] = spe.initial_concentration;
    }
    sim_ctx.set_data(tmp, sim_ctx_byte_size);
    delete[] tmp;
    
    tmp = new float[converter_byte_size / sizeof(float)];
    int size = sim_context.species_converters.size();
    tmp[0] = *((float*)(&size));
    counter = 1;
    for (auto& conv : sim_context.species_converters) {
        for (auto& r : conv.reduction_per_qm_sec) {
            tmp[counter++] = r;
        }
        for (auto& x : conv.conversion_ratio) {
            for (auto& y : x) {
                tmp[counter++] = y;
            }
        }
        tmp[counter++] = conv.temp;
    }
    assert(counter == converter_byte_size / sizeof(float));
    converter.set_data(tmp, converter_byte_size);
    delete[] tmp;

    tmp = new float[inserter_byte_size / sizeof(float)];
    size = sim_context.multi_inserters.size();
    tmp[0] = *((float*)(&size));
    counter = 1;
    for (auto& inserter : sim_context.multi_inserters) {
        tmp[counter++] = inserter.temp;
        for (auto& c : inserter.concentrations) {
            tmp[counter++] = c;
        }
        tmp[counter++] = inserter.velocities[0];
        tmp[counter++] = inserter.velocities[1];
    }
    assert(counter == inserter_byte_size / sizeof(float));
    inserter.set_data(tmp, inserter_byte_size);
    delete[] tmp;
}

GpuSimulationInfos::~GpuSimulationInfos()
{
}

void GpuSimulationInfos::get_simulation_context_data(void* data, uint32_t bytesize)
{
    sim_ctx.get_data(data, bytesize);
}

void GpuSimulationInfos::set_simulation_context_data(void* data, uint32_t bytesize)
{
    sim_ctx.set_data(data, bytesize);
}

VkBuffer GpuSimulationInfos::get_simulation_context_buffer()
{
    return sim_ctx.buffer;
}

uint32_t GpuSimulationInfos::get_simulation_context_byte_size()
{
    return sim_ctx_byte_size;
}

VkBuffer GpuSimulationInfos::get_converter_buffer()
{
    return converter.buffer;
}

uint32_t GpuSimulationInfos::get_converter_buffer_byte_size()
{
    return converter_byte_size;
}

VkBuffer GpuSimulationInfos::get_inserter_buffer()
{
    return inserter.buffer;
}

uint32_t GpuSimulationInfos::get_inserter_buffer_byte_size()
{
    return inserter_byte_size;
}
