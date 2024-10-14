#include "dds.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <istream>
#include <cassert>

// Reference: https://github.com/Microsoft/DirectXTex

struct DDS_PIXELFORMAT
{
  uint32_t dwSize;
  uint32_t dwFlags;
  uint32_t dwFourCC;
  uint32_t dwRGBBitCount;
  uint32_t dwRBitMask;
  uint32_t dwGBitMask;
  uint32_t dwBBitMask;
  uint32_t dwABitMask;
};

struct DDS_HEADER
{
  uint32_t        dwSize;
  uint32_t        dwFlags;
  uint32_t        dwHeight;
  uint32_t        dwWidth;
  uint32_t        dwPitchOrLinearSize;
  uint32_t        dwDepth;
  uint32_t        dwMipMapCount;
  uint32_t        dwReserved1[11];
  DDS_PIXELFORMAT ddspf;
  uint32_t        dwCaps;
  uint32_t        dwCaps2;
  uint32_t        dwCaps3;
  uint32_t        dwCaps4;
  uint32_t        dwReserved2;
};
static_assert(sizeof(DDS_HEADER) == 124, "DDS Header size mismatch");

struct DDS_HEADER_DXT10
{
  uint32_t dxgiFormat;
  uint32_t resourceDimension;
  uint32_t miscFlag;
  uint32_t arraySize;
  uint32_t miscFlags2;
};
static_assert(sizeof(DDS_HEADER_DXT10) == 20, "DDS DX10 Extended Header size mismatch");

DDS_PIXELFORMAT const DDSPF_DX10 = { sizeof(DDS_PIXELFORMAT), 0x00000004, 0x30315844, 0, 0, 0, 0, 0 };

uint32_t const DDS_MAGIC = ' ' << 24 | 'S' << 16 | 'D' << 8 | 'D'; // "DDS "
uint32_t const DDS_HEADER_FLAGS_PITCH = 0x00000008;
uint32_t const DDS_SURFACE_FLAGS_TEXTURE = 0x00001000; // DDSCAPS_TEXTURE
uint32_t const DDS_RESOURCE_DIMENSION_TEXTURE2D = 3;
uint32_t const DDS_RESOURCE_DIMENSION_TEXTURE3D = 4;

bool SaveDDS(
  const std::filesystem::path& inFilePath,
  unsigned inFormat, unsigned inTexelSizeInBytes,
  unsigned inWidth, unsigned inHeight, void const* outData
)
{
  std::ofstream file;
  file.open(inFilePath, std::ios::out | std::ios::trunc | std::ios::binary);
  if (!file.is_open())
    return false;

  // Write DDS_MAGIC
  file.write(reinterpret_cast<const char*>(&DDS_MAGIC), sizeof(DDS_MAGIC));

  // DDSD_CAPS | DDSD_HEIGHT | DDSD_WIDTH | DDSD_PIXELFORMAT 
  uint32_t const DDS_HEADER_FLAGS_TEXTURE = 0x00001007;

  // Write DDS_HEADER
  DDS_HEADER header;
  memset(&header, 0, sizeof(DDS_HEADER));
  header.dwSize = sizeof(DDS_HEADER);
  header.dwFlags = DDS_HEADER_FLAGS_TEXTURE | DDS_HEADER_FLAGS_PITCH;
  header.dwHeight = inHeight;
  header.dwWidth = inWidth;
  header.dwDepth = 1;
  header.dwMipMapCount = 1;
  header.dwPitchOrLinearSize = inWidth * inTexelSizeInBytes;
  header.ddspf = DDSPF_DX10;
  header.dwCaps = DDS_SURFACE_FLAGS_TEXTURE;
  file.write(reinterpret_cast<const char*>(&header), sizeof(DDS_HEADER));

  // Write DDS_HEADER_DX10
  DDS_HEADER_DXT10 header_dxt10;
  memset(&header_dxt10, 0, sizeof(DDS_HEADER_DXT10));
  header_dxt10.dxgiFormat = inFormat;
  header_dxt10.resourceDimension = DDS_RESOURCE_DIMENSION_TEXTURE2D;
  header_dxt10.arraySize = 1;
  file.write(reinterpret_cast<const char*>(&header_dxt10), sizeof(DDS_HEADER_DXT10));

  // Write image data
  file.write(reinterpret_cast<const char*>(outData), inWidth * inHeight * inTexelSizeInBytes);

  file.close();

  return true;
}

bool SaveDDS(
  const std::filesystem::path& inFilePath,
  unsigned inFormat, unsigned inTexelSizeInBytes,
  unsigned inWidth, unsigned inHeight, unsigned inDepth,
  void const* outData
)
{
  std::ofstream file;
  file.open(inFilePath, std::ios::out | std::ios::trunc | std::ios::binary);
  if (!file.is_open())
    return false;

  // Write DDS_MAGIC
  file.write(reinterpret_cast<const char*>(&DDS_MAGIC), sizeof(DDS_MAGIC));

  // DDSD_CAPS | DDSD_HEIGHT | DDSD_WIDTH | DDSD_DEPTH | DDSD_PIXELFORMAT 
  uint32_t const DDS_HEADER_FLAGS_TEXTURE = 0x00801007;

  // Write DDS_HEADER
  DDS_HEADER dds_header;
  memset(&dds_header, 0, sizeof(DDS_HEADER));
  dds_header.dwSize = sizeof(DDS_HEADER);
  dds_header.dwFlags = DDS_HEADER_FLAGS_TEXTURE | DDS_HEADER_FLAGS_PITCH;
  dds_header.dwHeight = inHeight;
  dds_header.dwWidth = inWidth;
  dds_header.dwDepth = inDepth;
  dds_header.dwMipMapCount = 1;
  dds_header.dwPitchOrLinearSize = inWidth * inTexelSizeInBytes;
  dds_header.ddspf = DDSPF_DX10;
  dds_header.dwCaps = DDS_SURFACE_FLAGS_TEXTURE;
  file.write(reinterpret_cast<const char*>(&dds_header), sizeof(DDS_HEADER));

  // Write DDS_HEADER_DX10
  DDS_HEADER_DXT10 dds_header_dxt10;
  memset(&dds_header_dxt10, 0, sizeof(DDS_HEADER_DXT10));
  dds_header_dxt10.dxgiFormat = inFormat;
  dds_header_dxt10.resourceDimension = DDS_RESOURCE_DIMENSION_TEXTURE3D;
  dds_header_dxt10.arraySize = 1;
  file.write(reinterpret_cast<const char*>(&dds_header_dxt10), sizeof(DDS_HEADER_DXT10));

  // Write image data
  const size_t data_size = inWidth * inHeight * inDepth * inTexelSizeInBytes;
  file.write(reinterpret_cast<const char*>(outData), data_size);

  file.close();

  return true;
}

bool LoadDDS(
  const std::filesystem::path& inFilePath,
  unsigned* outFormat, unsigned* outTexelSizeInBytes, unsigned* outWidth, unsigned* outHeight, unsigned* outDepth,
  std::vector<uint8_t>& outData
  )
{
  std::ifstream file(inFilePath, std::ios::in | std::ios::binary);

  if (!file.is_open())
    return false;

  // Read DDS_MAGIC
  uint32_t dds_magic = 0;
  file.read(reinterpret_cast<char*>(&dds_magic), sizeof(DDS_MAGIC));
  assert(dds_magic == DDS_MAGIC);

  // Read DDS_HEADER
  DDS_HEADER dds_header;
  file.read(reinterpret_cast<char*>(&dds_header), sizeof(DDS_HEADER));

  // Read DDS_HEADER_DXT10
  DDS_HEADER_DXT10 dds_header_dxt10;
  file.read(reinterpret_cast<char*>(&dds_header_dxt10), sizeof(DDS_HEADER_DXT10));

  *outHeight = dds_header.dwHeight;
  *outWidth = dds_header.dwWidth;
  *outDepth = dds_header.dwDepth;
  *outTexelSizeInBytes = dds_header.dwPitchOrLinearSize / *outWidth;
  *outFormat = dds_header_dxt10.dxgiFormat;

  // The following options are restricted in this function
  assert(dds_header.dwMipMapCount == 1);
  assert(dds_header_dxt10.arraySize == 1);

  // Read image data
  size_t data_size = (*outWidth) * (*outHeight) * (*outDepth) * (*outTexelSizeInBytes);
  outData.resize(data_size);
  file.read(reinterpret_cast<char*>(outData.data()), data_size);

  file.close();

  return true;
}