#pragma once
#include "util.h"

struct FinalComposeResource
{
  glm::mat4x4 mInvProj;
  glm::mat4x4 mInvView;
  float mFroxelNear;
  float mFroxelFar;

  util::Texture3D<glm::vec3, util::ClampPolicy3D>* mLuminance;
  util::Texture3D<glm::vec3, util::ClampPolicy3D>* mTransmittance;
  util::Image2D<glm::vec3>* mOutImage;
};

class FinalComposeOperator
{
public:
  FinalComposeOperator(const FinalComposeResource& inResource)
    : mResource(inResource)
  {
  }

  void Exec(glm::uvec3 inPos, glm::uvec3 inRes) const;

private:
  const FinalComposeResource& mResource;
};

void FinalComposeOperator::Exec(glm::uvec3 inPos, glm::uvec3 inRes) const
{
  const glm::vec2 screen_pos(inPos.x, inPos.y);
  const glm::vec2 screen_res(inRes.x, inRes.y);

  // [0,1]x[0,1]
  const glm::vec2 uv = screen_pos / (screen_res - 1.0f);

  // normal device coordinates [-1,1] x [1,-1]
  glm::vec2 pos = glm::vec2(2.0f, -2.0f) * uv + glm::vec2(-1.0f, 1.0f);

  const float froxel_near = mResource.mFroxelNear;
  const float froxel_far = mResource.mFroxelFar;
  const uint32_t slice_count = mResource.mLuminance->GetDepth();

  const float depth = froxel_far;
  const uint32_t index = util::DepthToSlice(depth, froxel_near, froxel_far, slice_count);

  const float z = static_cast<float>(index) / (slice_count - 1);
  glm::vec3 luminance = mResource.mLuminance->Sample(glm::vec3(uv.x, uv.y, z));
  glm::vec3 transmittance = mResource.mTransmittance->Sample(glm::vec3(uv.x, uv.y, z));

  glm::vec3 out_color = glm::vec3(0);
  out_color = out_color * transmittance + luminance;
  assert(!glm::any(glm::isnan(out_color) || glm::isinf(out_color)));
  out_color = util::TonemappingACES(out_color);
  out_color = util::LinearToSRGBApprox(out_color);

  mResource.mOutImage->Store(inPos.xy(), out_color);
}

