#pragma once
#include "util.h"

//----------------------------------------------------------------------

struct RaymarchingResource
{
  glm::mat4x4 mInvView;
  util::FroxelInfo* mFroxelInfo;

  util::Texture3D<glm::vec4, util::ClampPolicy3D>* mMaterial;
  util::Texture3D<glm::vec3, util::ClampPolicy3D>* mLight;
  util::Image3D<glm::vec3>* mOutLuminance;
  util::Image3D<glm::vec3>* mOutTransmittance;
};

class RaymarchingOperator
{
public:
  RaymarchingOperator(const RaymarchingResource& inResource)
    : mResource(inResource)
  {
  }

  struct Throughput
  {
    glm::vec3 transmittance = glm::vec3(1.0f);
    glm::vec3 luminance = glm::vec3(0.0f);
  };

  struct Material
  {
    glm::vec3 scatteringCoeff;
    glm::vec3 mieAsymmetry;
  };

  void Exec(glm::uvec3 inPos, glm::uvec3 inRes) const;

private:
  const RaymarchingResource& mResource;
};

void RaymarchingOperator::Exec(glm::uvec3 inPos, glm::uvec3 inRes) const
{
  assert(mResource.mOutLuminance->GetWidth()  == inRes.x);
  assert(mResource.mOutLuminance->GetHeight() == inRes.y);
  assert(mResource.mOutLuminance->GetDepth() > 1);
  const glm::vec3 screen_pos(inPos.x, inPos.y, inPos.z);
  const glm::vec3 screen_res(inRes.x, inRes.y, inRes.z);

  // [0,1]x[0,1]
  const glm::vec2 uv = screen_pos.xy() / (screen_res.xy() - 1.0f);

  glm::vec3 view_vector = mResource.mFroxelInfo->ComputeViewVector(inPos);
  glm::uvec3 voxel_pos = glm::uvec3(inPos.x, inPos.y, 0);

  Throughput throughput;
  Material material;
  material.mieAsymmetry = glm::vec3(0.76f);

  const int32_t slice_count = mResource.mOutLuminance->GetDepth();
  for (int32_t z = 0; z < slice_count; ++z)
  {
    voxel_pos.z = z;
    const glm::vec3 uvw = glm::vec3(uv.x, uv.y, static_cast<float>(z) / (mResource.mOutLuminance->GetDepth() - 1.0f));

    const float depth_near = mResource.mFroxelInfo->SliceToDepth(z);
    const float depth_far = mResource.mFroxelInfo->SliceToDepth(z+1);
    const float ds = depth_far - depth_near;
    assert(ds > 0);
    const glm::vec3 view_pos = view_vector * depth_near;
    const glm::vec3 world_pos = (mResource.mInvView * glm::vec4(view_pos, 1)).xyz();
    const glm::vec3 world_dir = glm::normalize((mResource.mInvView * glm::vec4(view_vector, 0)).xyz());

    glm::vec3 mu_t, luminance;
    {
      // scattering coefficients
      glm::vec4 coeffs = mResource.mMaterial->Sample(uvw);
      // rayleigh
      glm::vec3 mu_sr = coeffs.w * util::RAYLEIGH_S;
      // +mie
      glm::vec3 mu_sm = coeffs.xyz();
      glm::vec3 mu_am = coeffs.xyz() * util::MIE_ABSORB_MUL;
      // total
      mu_t = glm::max(mu_sr + mu_sm + mu_am, glm::vec3(1e-6f));

      luminance = mResource.mLight->Sample(uvw);
    }
    glm::vec3 transmittance = glm::vec3(glm::exp(-mu_t * ds));
    glm::vec3 inscatter = luminance * (1.0f - transmittance);
    throughput.transmittance *= transmittance;
    throughput.luminance += inscatter * throughput.transmittance;

    mResource.mOutLuminance->Store(voxel_pos, glm::vec4(throughput.luminance, 0));
    mResource.mOutTransmittance->Store(voxel_pos, glm::vec4(throughput.transmittance, 0));
  }
}

