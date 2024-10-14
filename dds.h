#pragma once
#include <filesystem>
#include <cassert>

unsigned const DDS_FORMAT_R32G32B32A32_FLOAT = 2;
unsigned const DDS_FORMAT_R32G32B32_FLOAT = 6;
unsigned const DDS_FORMAT_R32G32_FLOAT = 16;
unsigned const DDS_FORMAT_R16G16_FLOAT = 34;
unsigned const DDS_FORMAT_R32_FLOAT = 41;
unsigned const DDS_FORMAT_R8_UNORM = 61;

// 2D
bool SaveDDS(
  const std::filesystem::path& inFilePath,
  unsigned inFormat, unsigned inTexelSizeInBytes,
  unsigned inWidth, unsigned inHeight,
  void const* outData
);

// 3D
bool SaveDDS(
  const std::filesystem::path& inFilePath,
  unsigned inFormat, unsigned inTexelSizeInBytes,
  unsigned inWidth, unsigned inHeight, unsigned depth,
  void const* outData
);

// 3D
bool LoadDDS(
  const std::filesystem::path& inFilePath,
  unsigned* outFormat, unsigned* outTexelSizeInBytes,
  unsigned* outWidth, unsigned* outHeight, unsigned* outDepth,
  std::vector<uint8_t>& outData
);

// 3D(template)
template<typename T>
bool LoadDDS(
  const std::filesystem::path& inFilePath,
  unsigned* outFormat,
  unsigned* outTexelSizeInBytes,
  unsigned* outWidth,
  unsigned* outHeight,
  unsigned* outDepth,
  std::vector<T>& outData)
{
  static_assert(std::is_trivially_copyable<T>::value, "T must be trivially copyable");

  std::vector<uint8_t> data;
  bool result = LoadDDS(inFilePath, outFormat, outTexelSizeInBytes, outWidth, outHeight, outDepth, data);
  if (!result)
  {
    return result;
  }

  T* ptr = reinterpret_cast<T*>(data.data());
  if (data.size() % sizeof(T) != 0)
  {
    assert(false && "Data size is not a multiple of the element size");
    return false;
  }
  if (reinterpret_cast<std::uintptr_t>(ptr) % alignof(T) != 0)
  {
    assert(false && "Data alignment is incorrect for the specified type");
    return false;
  }

  size_t count = data.size() / sizeof(T);
  outData.reserve(count);
  outData.assign(ptr, ptr + count);

  return true;
}