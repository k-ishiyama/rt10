#pragma once
#include "dds.h"
#include "fpng.h"

// #define GLM_FORCE_MESSAGES
#define GLM_FORCE_AVX2
#define GLM_FORCE_INLINE
#define GLM_FORCE_SWIZZLE
#define GLM_FORCE_LEFT_HANDED
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <omp.h>
#include <filesystem>
#include <array>
#include <vector>
#include <limits>

namespace util {

template<typename T> constexpr T PI = T(3.14159265358979323846L);

constexpr float MIE_ABSORB_MUL = 0.11f;
const glm::vec3 RAYLEIGH_S = 1e-2f * glm::vec3(5.78f, 13.6f, 33.1f);

// https://chilliant.blogspot.com/2012/08/srgb-approximations-for-hlsl.html
glm::vec3 LinearToSRGBApprox(const glm::vec3& inLinearColor)
{
  const glm::vec3& c = inLinearColor;
  return glm::max(glm::vec3(1.055f) * glm::pow(c, glm::vec3(5.0f / 12.0f)) - glm::vec3(0.055f), glm::vec3(0));
}

// log depth
float LinearDepthToLogDepth(float linearDepth, float near, float far)
{
  float base = far / near;
  return log2(linearDepth / near) / log2(base);
}

// log depth
float LogDepthToLinearDepth(float logDepth, float near, float far)
{
  float base = far / near;
  return near * pow(base, logDepth);
}

// depth slice
unsigned int DepthToSlice(float linearDepth, float near, float far, unsigned int sliceCount)
{
  assert(sliceCount > 0);
  float logDepth = LinearDepthToLogDepth(linearDepth, near, far);
  unsigned int sliceIndex = static_cast<unsigned int>(logDepth * sliceCount);
  return std::min(sliceIndex, sliceCount - 1);
}

// depth slice
float SliceToDepth(unsigned int sliceIndex, float near, float far, unsigned int sliceCount)
{
  assert(sliceCount > 0);
  float logDepth = static_cast<float>(sliceIndex) / (sliceCount - 1);
  return LogDepthToLinearDepth(logDepth, near, far);
}

glm::vec3
HenyeyGreenstein(float mu, const glm::vec3& g)
{
  glm::vec3 base = glm::abs(glm::vec3(1.0) + g * g - glm::vec3(2.0) * g * mu);
  glm::vec3 denom = ((glm::vec3(4.0) * glm::pi<glm::vec3::value_type>()) * glm::pow(base, glm::vec3(1.5)));
  return (glm::vec3(1.0) - g * g) / denom;
}

float
Rayleigh(float mu)
{
  return 3.0f / 4.0f / (4.0f * PI<float>) * (1.0f + mu * mu);
}


struct Ray
{
public:
  Ray(const glm::vec3& origin, const glm::vec3& direction)
      : mOrigin(origin), mDirection(glm::normalize(direction)) {}

  const glm::vec3& GetOrigin() const { return mOrigin; }
  const glm::vec3& GetDirection() const { return mDirection; }

private:
  glm::vec3 mOrigin;
  glm::vec3 mDirection;
};

struct AABB
{
public:
  AABB(const glm::vec3& min, const glm::vec3& max)
    : mMin(min), mMax(max) {}

  const glm::vec3& GetMin() const { return mMin; }
  const glm::vec3& GetMax() const { return mMax; }

private:
  glm::vec3 mMin = glm::vec3(std::numeric_limits<float>::max());
  glm::vec3 mMax = glm::vec3(std::numeric_limits<float>::lowest());
};

std::tuple<bool, float, float>
IntersectRayAABB(const Ray& inRay, const AABB& inAABB)
{
  glm::vec3 t_min = (inAABB.GetMin() - inRay.GetOrigin()) / inRay.GetDirection();
  glm::vec3 t_max = (inAABB.GetMax() - inRay.GetOrigin()) / inRay.GetDirection();
  glm::vec3 t1 = glm::min(t_min, t_max);
  glm::vec3 t2 = glm::max(t_min, t_max);
  float t_near = std::max({ t1.x, t1.y, t1.z });
  float t_far = std::min({ t2.x, t2.y, t2.z });

  // (false: no intersection, near distance, far distance)
  return { (t_near <= t_far), t_near, t_far };
};

template<typename T> T WrapUV(const T& inUV)
{
  return glm::fract(inUV);
}

template<typename T> T ClampUV(const T& inUV)
{
  return glm::clamp(inUV, T(0), T(1));
}

glm::vec4 BilinearInterpWeights(const glm::vec2& inPosition)
{
  glm::vec2 f = glm::fract(inPosition);
  glm::vec4 w = f.xyxy() * glm::vec4(1,1,-1,-1) + glm::vec4(0,0,1,1);
  w = w.zxzx() * w.wwyy();
  return w;
}

template<typename T>
T BilinearInterpolate(
  const glm::vec2& inFracPosition,
  const T& inV00, const T& inV01, const T& inV10, const T& inV11)
{
  glm::vec4 w = BilinearInterpWeights(inFracPosition);  
  return w.x * inV00 + w.y * inV01 + w.z * inV10 + w.w * inV11;
}

//-----------------------------------------------------------------------------

template<typename Type, uint32_t NumDimensions> struct Image
{
  using type = Type;
  using iterator = typename std::vector<Type>::iterator;
  using const_iterator = typename std::vector<Type>::const_iterator;

  Image(Type* inBegin, Type* inEnd, const std::array<uint32_t, NumDimensions>& inDimension)
    : mDimensions(inDimension)
    , mElements(inBegin, inEnd)
  {
    assert(inEnd >= inBegin);
    
    uint32_t expected_size = 1;
    for (auto d : mDimensions)
    {
      expected_size *= d;
    }
    uint32_t actual_size = static_cast<uint32_t>(inEnd - inBegin);
    assert(expected_size == actual_size);
  }

  Image(Image<Type, NumDimensions>&& ioOther)
  {
    this->mElements = std::move(ioOther.mElements);
    this->mDimensions = std::move(ioOther.mDimensions);
  }
  Image(const std::array<uint32_t, NumDimensions>& inDimension)
  {
    for (uint32_t i = 0; i < NumDimensions; ++i)
    {
      SetDimension(i, inDimension[i]);
    }

    Resize();
  }

  uint32_t GetNumPixels() const { return static_cast<uint32_t>(mElements.size()); }
  Type* GetData() { return mElements.data(); }
  const Type* GetData() const { return mElements.data(); }
  std::vector<Type>& GetVector() { return mElements; }
  const std::vector<Type>& GetVector() const { return mElements; }
  uint32_t GetDimension(uint32_t i) const { return mDimensions[i]; }
  void SetDimension(uint32_t i, uint32_t length) { mDimensions[i] = length; }
  uint32_t GetSizeInBytes() const { return mElements.size() * sizeof(Type); }
  void Resize() {
    uint32_t size = 1;
    for (auto d : mDimensions)
    {
      size *= d;
    }
    mElements.resize(size);
  }
  void Fill(const Type& v)
  {
    std::fill(mElements.begin(), mElements.end(), v);
  }
  void Copy(const Image<Type, NumDimensions>& src)
  {
#ifdef _DEBUG
    assert(src.GetSize() == GetSizeInBytes());
#endif
    std::copy(mElements.begin(), mElements.end(), src.GetData().begin());
  }

  iterator begin() { return mElements.begin(); }
  iterator end() { return mElements.end(); }
  const_iterator begin() const { return mElements.begin(); }
  const_iterator end() const { return mElements.end(); }
  uint32_t GetIndex(iterator& it) { return static_cast<uint32_t>(std::distance(begin(), it)); }
  uint32_t GetIndex(const_iterator& it) const { return static_cast<uint32_t>(std::distance(begin(), it)); }

private:
  std::array<uint32_t, NumDimensions> mDimensions = {};
  std::vector<Type> mElements;
};

//-----------------------------------------------------------------------------

template<typename T> struct Image2D : Image<T, 2>
{
private:
  using Super = Image<T, 2>;

public:
  Image2D(const glm::uvec2& inSize)
    : Image<T, 2>({inSize.x, inSize.y})
  {}
  uint32_t GetWidth() const { return this->GetDimension(0); }
  uint32_t GetHeight() const { return this->GetDimension(1); }
  
  void Store(const glm::uvec2& pos, const T& v)
  {
#ifdef _DEBUG
    assert(pos.x < GetWidth());
    assert(pos.y < GetHeight());
#endif    
    this->GetVector()[GetWidth() * pos.y + pos.x] = v;
  }

  const T& Load(const glm::uvec2& pos) const
  {
#ifdef _DEBUG
    assert(pos.x < GetWidth());
    assert(pos.y < GetHeight());
#endif    
    return this->GetVector()[GetWidth() * pos.y + pos.x];
  }

  glm::uvec2 GetDimensions() const
  {
    return { this->GetDimension(0), this->GetDimension(1) };
  }

  glm::uvec2 GetPosition(int index) const
  {
    int x = index % GetWidth();
    int y = index / GetWidth();
#ifdef _DEBUG
    assert(x < (int)GetWidth());
    assert(y < (int)GetHeight());
#endif
    return { x, y };
  }
};

template<typename T, typename PositionPolicy> struct Texture2D : Image2D<T>
{
  Texture2D(const glm::uvec2& inSize)
    : Image2D<T>(inSize)
  {
  }


  T Sample(const glm::vec2& inUV) const
  {
    glm::vec2 half_pixel_offset = glm::vec2(0.5f);
    glm::vec2 size(this->GetWidth(), this->GetHeight());
    glm::vec2 pos = PositionPolicy::UV(inUV - half_pixel_offset / size) * size;
    glm::ivec2 grid_pos = glm::ivec2(pos);
    
    glm::ivec2 p00 = PositionPolicy::Position(grid_pos, size);
    glm::ivec2 p01 = PositionPolicy::Position({ grid_pos.x + 1, grid_pos.y }, size);
    glm::ivec2 p10 = PositionPolicy::Position({ grid_pos.x, grid_pos.y + 1 }, size);
    glm::ivec2 p11 = PositionPolicy::Position({ grid_pos.x + 1, grid_pos.y + 1 }, size);

    const T& c00 = this->Load(p00);
    const T& c01 = this->Load(p01);
    const T& c10 = this->Load(p10);
    const T& c11 = this->Load(p11);

    glm::vec2 frac_pos = glm::fract(pos);
    return BilinearInterpolate<T>(frac_pos, c00, c01, c10, c11);
  }

  T SamplePoint(const glm::vec2& inUV) const
  {
    glm::vec2 size(this->GetWidth(), this->GetHeight());
    glm::vec2 pos = inUV * size;
    glm::ivec2 grid_pos = glm::ivec2(pos);
    
    glm::ivec2 p = PositionPolicy::Position(grid_pos, size);
    const T& c = this->Load(p);
    return c;
  }

};

struct ClampPolicy2D {
  static glm::vec2 UV(const glm::vec2& inUV)
  {
    return glm::clamp(inUV, glm::vec2(0), glm::vec2(1));
  }

  static glm::ivec2 Position(const glm::ivec2& inPos, const glm::ivec2& size)
  {
    return glm::clamp(inPos, glm::ivec2(0), size - 1);
  }
};

struct WrapPolicy2D {
  static glm::vec2 UV(const glm::vec2& inUV)
  {
    return glm::mod(inUV, glm::vec2(1));
  }

  static glm::ivec2 Position(const glm::ivec2& inPos, const glm::ivec2& size)
  {
    return glm::mod(glm::vec2(inPos), glm::vec2(size));
  }
};

//-----------------------------------------------------------------------------

template<typename T> struct Image3D : Image<T, 3>
{
private:
  using Super = Image<T, 3>;

public:
  using ValueType = T;
  Image3D(T* inBegin, T* inEnd, const std::array<uint32_t, 3> inSize)
    : Image<T, 3>(inBegin, inEnd, inSize)
  {}
  Image3D(T* inBegin, T* inEnd, const glm::uvec3& inSize)
    : Image<T, 3>(inBegin, inEnd, {inSize.x, inSize.y, inSize.z})
  {}
  Image3D(const glm::uvec3& inSize)
    : Image<T, 3>({inSize.x, inSize.y, inSize.z})
  {}
  Image3D(Image3D<T>&& inOther)
    : Image<T, 3>(std::move(inOther))
  {}

  int32_t GetWidth() const { return this->GetDimension(0); }
  int32_t GetHeight() const { return this->GetDimension(1); }
  int32_t GetDepth() const { return this->GetDimension(2); }

  void Store(const glm::ivec3& pos, const T& v)
  {
#ifdef _DEBUG
    assert(pos.x < GetWidth());
    assert(pos.y < GetHeight());
    assert(pos.z < GetDepth());
#endif    
    const uint32_t index = (pos.z * GetHeight() + pos.y) * GetWidth() + pos.x;
    this->GetVector()[index] = v;
  }

  const T& Load(const glm::ivec3& pos) const
  {
#ifdef _DEBUG
    assert(pos.x < GetWidth());
    assert(pos.y < GetHeight());
    assert(pos.z < GetDepth());
#endif    
    const uint32_t index = (pos.z * GetHeight() + pos.y) * GetWidth() + pos.x;
    return this->GetVector()[index];
  }

  glm::uvec3 GetDimensions() const
  {
    return { this->GetDimension(0), this->GetDimension(1), this->GetDimension(2) };
  }

  glm::uvec3 GetPosition(int index) const
  {
    int z = index / (GetWidth() * GetHeight());
    int y = (index % (GetWidth() * GetHeight())) / GetWidth();
    int x = index % GetWidth();
#ifdef _DEBUG
    assert(x < GetWidth());
    assert(y < GetHeight());
    assert(z < GetDepth());
#endif
    return { x, y, z };
  }
};

template<typename T, typename PositionPolicy> struct Texture3D : Image3D<T>
{
  Texture3D(const glm::uvec3& inSize)
    : Image3D<T>(inSize)
  {}
  Texture3D(Image3D<T>&& inOther)
    : Image3D<T>(std::move(inOther))
  {}

  template<typename S>
  S Sample(const glm::vec3& inUVW) const
  {
    const glm::vec3 half_pixel_offset = glm::vec3(0.5f);
    const glm::vec3 size(this->GetWidth(), this->GetHeight(), this->GetDepth());
    glm::vec3 pos = PositionPolicy::UV(inUVW - half_pixel_offset / size) * size;
    glm::ivec3 ip = glm::ivec3(pos);

    const S c000 = static_cast<S>(this->Load(PositionPolicy::Position({ ip.x,   ip.y,   ip.z }, size)));
    const S c010 = static_cast<S>(this->Load(PositionPolicy::Position({ ip.x + 1, ip.y,   ip.z }, size)));
    const S c100 = static_cast<S>(this->Load(PositionPolicy::Position({ ip.x,   ip.y + 1, ip.z }, size)));
    const S c110 = static_cast<S>(this->Load(PositionPolicy::Position({ ip.x + 1, ip.y + 1, ip.z }, size)));
    const S c001 = static_cast<S>(this->Load(PositionPolicy::Position({ ip.x,   ip.y,   ip.z + 1 }, size)));
    const S c011 = static_cast<S>(this->Load(PositionPolicy::Position({ ip.x + 1, ip.y,   ip.z + 1 }, size)));
    const S c101 = static_cast<S>(this->Load(PositionPolicy::Position({ ip.x,   ip.y + 1, ip.z + 1 }, size)));
    const S c111 = static_cast<S>(this->Load(PositionPolicy::Position({ ip.x + 1, ip.y + 1, ip.z + 1 }, size)));

    glm::vec3 frac_pos = glm::fract(pos);
    S c0 = BilinearInterpolate<S>(frac_pos.xy(), c000, c010, c100, c110);
    S c1 = BilinearInterpolate<S>(frac_pos.xy(), c001, c011, c101, c111);
    return c0 * (1 - frac_pos.z) + c1 * frac_pos.z;
  }

  T Sample(const glm::vec3& inUVW) const
  {
    return this->Sample<T>(inUVW);
  }

};

struct ClampPolicy3D {
  static glm::vec3 UV(const glm::vec3& inUV)
  {
    return glm::clamp(inUV, glm::vec3(0), glm::vec3(1));
  }

  static glm::ivec3 Position(const glm::ivec3& inPos, const glm::ivec3& size) {
    return glm::clamp(inPos, glm::ivec3(0), size - 1);
  }
};

struct WrapPolicy3D {
  static glm::vec3 UV(const glm::vec3& inUV)
  {
    return glm::mod(inUV, glm::vec3(1));
  }

  static glm::ivec3 Position(const glm::ivec3& inPos, const glm::ivec3& size) {
    return glm::mod(glm::vec3(inPos), glm::vec3(size));
  }
};

bool SaveImage(const std::filesystem::path& path, const Image2D<glm::u8vec4>& image)
{
  auto s = path.string();
  auto c = s.c_str();
  const char* filename = s.c_str();

  bool result = fpng::fpng_encode_image_to_file(
    filename, reinterpret_cast<const void*>(image.GetData()),
    image.GetWidth(), image.GetHeight(), 4);
  return result;
}

struct Film
{
public:
  Film(const glm::uvec2& inSize)
    : mImage(inSize)
  {
  }
  uint32_t GetWidth() const { return mImage.GetWidth(); }
  uint32_t GetHeight() const { return mImage.GetHeight(); }
  const Image2D<glm::u8vec4>& GetImage() const { return mImage; }

  void Develop(const Image2D<glm::vec3>& inImage)
  {
    assert(mImage.GetWidth() == inImage.GetWidth());
    assert(mImage.GetHeight() == inImage.GetHeight());
    int w = static_cast<int>(mImage.GetWidth());
    int h = static_cast<int>(mImage.GetHeight());

#pragma omp parallel for schedule(static)
    for (int i = 0; i < w * h; ++i)
    {
      int x = i % w;
      int y = i / w;

      const glm::vec3& in = inImage.Load({ x, y });

      glm::vec3 c = glm::clamp(in, glm::vec3(0.0f), glm::vec3(1.0f));

      // [0,1] -> [0,255]
      glm::u8vec4 out;
      out[0] = static_cast<uint8_t>(c[0] * 255);
      out[1] = static_cast<uint8_t>(c[1] * 255);
      out[2] = static_cast<uint8_t>(c[2] * 255);
      out[3] = 255;
      mImage.Store({ x, y }, out);
    }
  }

private:
  Image2D<glm::u8vec4> mImage;
};

void Export(const std::string& inFileName, const Image2D<glm::u8vec4>& inImage)
{
  const std::filesystem::path current_path = std::filesystem::current_path();
  std::string filepath = (current_path / inFileName).string();
  SaveImage(filepath, inImage);
}

// * O’neill, M. E. (2014). PCG: A family of simple fast space-efficient statistically good algorithms for random number generation. ACM Transactions on Mathematical Software.
// * Vigna, S. (2016). An experimental exploration of Marsaglia's xorshift generators, scrambled. ACM Transactions on Mathematical Software (TOMS), 42(4), 1-23.
// * https://en.wikipedia.org/wiki/Xorshift#xorshift*
//
// 64 state bits, uint32_t output, period 2^64 - 1
struct XorShift64Star
{
  XorShift64Star() = default;
  XorShift64Star(uint64_t inSeed)
    : mState(inSeed)
  {   
  }

  uint32_t Next()
  {
    const uint64_t result = mState * MULTIPLIER;
    Advance();
    const uint32_t bit_shift = 8 * (sizeof(uint64_t) - sizeof(uint32_t));
    return result >> bit_shift;
  }

  // [0, 1]
  float NextFloat01()
  {
    constexpr uint32_t max_value = std::numeric_limits<uint32_t>::max();
    return Next() / static_cast<float>(max_value);
  }

private:
  void Advance()
  {
    mState ^= mState >> A;
    mState ^= mState << B;
    mState ^= mState >> C;
  }

  static constexpr uint32_t A = 12;
  static constexpr uint32_t B = 25;
  static constexpr uint32_t C = 27;
  static constexpr uint64_t MULTIPLIER = 0x2545F4914F6CDD1DULL;

  uint64_t mState = 1;
};

template <typename Operator, typename Resource, uint32_t Nx, uint32_t Ny, uint32_t Nz>
struct Compute
{
public:
  Compute(const Resource& inResource)
    : mOperator(Operator(inResource))
  {}

  void Exec() const
  {
    const int total = Nx * Ny * Nz;

    #pragma omp parallel for schedule(static)
    for (int idx = 0; idx < total; ++idx)
    {
      int x = idx % Nx;
      int y = (idx / Nx) % Ny;
      int z = idx / (Nx * Ny);

      mOperator.Exec({ x, y, z }, { Nx, Ny, Nz });
    }
  }

private:
  Operator mOperator;
};

struct FroxelInfo
{
public:
  constexpr FroxelInfo(
      const glm::mat4x4& inInvProj, float inFroxelNear, float inFroxelFar,
      const glm::uvec3& inResolution)
    : mResolution(inResolution)
    , mInvProj(inInvProj)
    , mFroxelNear(inFroxelNear)
    , mFroxelFar(inFroxelFar)
  {}
  
  glm::vec3 ComputeViewVector(const glm::uvec2 inPos) const
  {
    const glm::vec2 screen_pos(inPos);
    const glm::vec2 screen_res(mResolution.xy());

    // [0,1]x[0,1]
    const glm::vec2 uv = screen_pos / (screen_res - 1.0f);

    // normal device coordinates [-1,1] x [1,-1]
    glm::vec2 pos = glm::vec2(2.0f, -2.0f) * uv + glm::vec2(-1.0f, 1.0f);

    // normal device coordinates -> view coordinates(homogeneous)
    glm::vec4 homogeneous_view_vector = mInvProj * glm::vec4(pos.x, pos.y, -1, 1);

    // unnormalized view vector (vx, vy, 1)
    glm::vec3 view_vector = homogeneous_view_vector.xyz() / homogeneous_view_vector.z;
    return view_vector;
  }

  glm::vec3 ComputeViewPosition(const glm::uvec2& inScreenPos, float inDepth)
  {
    return ComputeViewVector(inScreenPos) * inDepth;
  }

  float SliceToDepth(uint32_t inSliceIndex) const
  {
    return util::SliceToDepth(inSliceIndex, mFroxelNear, mFroxelFar, mResolution.z);
  }

  uint32_t DepthToSlice(float inDepth) const
  {
    return util::DepthToSlice(inDepth, mFroxelNear, mFroxelFar, mResolution.z);
  }

  const glm::uvec3& GetResolution() const { return mResolution; }

  float GetFroxelNear() const { return mFroxelNear; }
  float GetFroxelFar() const { return mFroxelFar; }

private:
  glm::uvec3 mResolution;
  glm::mat4x4 mInvProj;
  float mFroxelNear;
  float mFroxelFar;
};

struct Frustum
{
public:
  Frustum(float inFoV, float inAspectRatio, float inNear, float inFar)
    : mFoV(inFoV), mAspectRatio(inAspectRatio), mNear(inNear), mFar(inFar)
  {
    SetViewProjMatrix(glm::perspective(inFoV, inAspectRatio, inNear, inFar));
  }

  const glm::mat4x4& GetViewProj() const { return mViewProj; }
  const glm::mat4x4& GetProjView() const { return mProjView; }

  float GetFoV() const { return mFoV; }
  float GetAspectRatio() const { return mAspectRatio; }
  float GetNear() const { return mNear; }
  float GetFar() const { return mFar; }

private:
  void SetViewProjMatrix(const glm::mat4x4& inViewProj)
  {
    mViewProj = inViewProj;
    mProjView = glm::inverse(inViewProj);
  }

  float mFoV;
  float mAspectRatio;
  float mNear;
  float mFar;

  glm::mat4x4 mViewProj;
  glm::mat4x4 mProjView;
};

struct Plane
{
  glm::vec3 position;
  glm::vec2 size;
};

glm::vec3 GenerateRandomPointOnPlane(const glm::vec2& inRandom01, const Plane& plane)
{
  float x = plane.size[0] * (inRandom01[0] - 0.5f);
  float z = plane.size[1] * (inRandom01[1] - 0.5f);
  return plane.position + glm::vec3(x, 0.0f, z);
}

// https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
glm::vec3 TonemappingACES(const glm::vec3& x)
{
  const float a = 2.51f;
  const float b = 0.03f;
  const float c = 2.43f;
  const float d = 0.59f;
  const float e = 0.14f;
  return (x * (a * x + b)) / (x * (c * x + d) + e);
}

template<typename T>
std::unique_ptr<Image3D<T>> CreateImageFromDDS(const std::filesystem::path& inFilePath)
{
  static_assert(std::is_trivially_copyable<T>::value, "T must be trivially copyable");

  uint32_t format = {};
  uint32_t texel_size_in_bytes = {};
  uint32_t width = {};
  uint32_t height = {};
  uint32_t depth = {};
  std::vector<uint8_t> data;
  bool result = LoadDDS(inFilePath, &format, &texel_size_in_bytes, &width, &height, &depth, data);
  if (!result)
  {
    return nullptr;
  }

  if (data.size() % sizeof(T) != 0)
  {
    assert(false && "Data size is not a multiple of the element size.\n");
    return nullptr;
  }

  size_t element_count = data.size() / sizeof(T);
  T* ptr = reinterpret_cast<T*>(data.data());

  if (reinterpret_cast<std::uintptr_t>(ptr) % alignof(T) != 0)
  {
    assert(false && "Data alignment is incorrect");
    return nullptr;
  }

  const std::array<uint32_t, 3> size = { width, height, depth };
  return std::make_unique<Image3D<T>>(ptr, ptr + element_count, size);
}

template<typename T, typename PositionPolicy>
std::unique_ptr<Texture3D<T, PositionPolicy>>
CreateTextureFromDDS(const std::filesystem::path& inFilePath)
{
  std::unique_ptr<Image3D<T>> img = CreateImageFromDDS<T>(inFilePath);
  return std::make_unique<Texture3D<T, PositionPolicy>>(std::move(*img));
}

} // namespace util