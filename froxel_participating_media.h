#pragma once
#include "util.h"

//----------------------------------------------------------------------

struct ParticipatingMediaResource
{
  glm::mat4x4 mInvView;
  util::FroxelInfo* mFroxelInfo;
  float mElapsedTime;
  util::Texture3D<uint8_t, util::WrapPolicy3D>* mPerlinWorleyShape;
  util::Texture3D<uint8_t, util::WrapPolicy3D>* mPerlinWorleyErosion;
  util::Image3D<glm::vec4>* mOutImage;
};

class ParticipatingMediaOperator
{
public:
  ParticipatingMediaOperator(const ParticipatingMediaResource& inResource)
    : mResource(inResource)
  {
  }

  void Exec(glm::uvec3 inPos, glm::uvec3 inRes) const;

private:
  const ParticipatingMediaResource& mResource;
};

void ParticipatingMediaOperator::Exec(glm::uvec3 inPos, glm::uvec3 inRes) const
{
  glm::vec3 pos(inPos);
  glm::vec3 res(inRes);
  glm::vec3 uvw = pos / (res - glm::vec3(1));

  const glm::vec3 view_vector = mResource.mFroxelInfo->ComputeViewVector(inPos);
  const float depth_near = mResource.mFroxelInfo->SliceToDepth(inPos.z);
  const float depth_far = mResource.mFroxelInfo->SliceToDepth(inPos.z+1);

  glm::vec3 mie_density = glm::vec3(0);
  float rayleigh_density = 0;
  const uint32_t index = (inPos.z * inRes.y + inPos.y) * inRes.x + inPos.x;
  util::XorShift64Star random(index);
  const int sample_per_voxel = 4;
  for (int i = 0; i < sample_per_voxel; ++i)
  {
    const float t = random.NextFloat01();
    const float depth = depth_near* t + depth_far * (1.0f - t);

    const glm::vec3 view_pos = mResource.mFroxelInfo->ComputeViewPosition(inPos.xy(), depth);
    const glm::vec3 world_pos = (mResource.mInvView * glm::vec4(view_pos, 1)).xyz();

    const glm::vec3 mie_coeff = glm::vec3(0.92f, 0.93f, 1.0f);
    glm::vec3 density = glm::vec3(0);

    // ambient
    const float amb_density = 0.04f;
    const float amb_height_attenuation = glm::exp(-0.2f * glm::max(world_pos.y, 0.0f));
    rayleigh_density += 0.07f * amb_height_attenuation;
    mie_density += 0.03f * mie_coeff * amb_height_attenuation;

    // clouds
    const float noise_height_attenuation = glm::exp(-glm::max(world_pos.y + 0.0f, 0.0f));
    if (noise_height_attenuation > 0.001f)
    {
      float time = 0.6f * mResource.mElapsedTime;
      glm::vec3 sp = 0.06f * world_pos + glm::vec3(0, 0, time);
      glm::vec3 sp2 = 0.07f * world_pos + glm::vec3(0, 0, 1.4142f * time);
      float noise_high = mResource.mPerlinWorleyShape->Sample<float>(sp) / 255.0f;
      float noise_low = mResource.mPerlinWorleyErosion->Sample<float>(sp2) / 255.0f;
      float noise = pow(noise_high, 16.0f + 16.0f * pow(noise_low, 16.0f));
      mie_density += noise * mie_coeff * noise_height_attenuation;
    }
  }
  mie_density /= static_cast<float>(sample_per_voxel);
  rayleigh_density /= static_cast<float>(sample_per_voxel);

  mResource.mOutImage->Store(inPos, glm::vec4(mie_density, rayleigh_density));
}

