#pragma once
#include "util.h"

//----------------------------------------------------------------------

struct PointLight
{
  glm::vec3 position;
  glm::vec3 luminance;
};

struct LightResource
{
  glm::mat4x4 mInvView;
  util::FroxelInfo* mFroxelInfo;
  util::Texture3D<glm::vec4, util::ClampPolicy3D>* mMaterial;
  util::Texture3D<glm::vec3, util::ClampPolicy3D>* mOutImage;
  std::vector<PointLight> mPointLights;
};

class LightOperator
{
public:
  LightOperator(const LightResource& inResource)
    : mResource(inResource)
  {
  }

  void Exec(glm::uvec3 inPos, glm::uvec3 inRes) const;

private:

  const LightResource& mResource;
};

void LightOperator::Exec(glm::uvec3 inPos, glm::uvec3 inRes) const
{
  assert(inRes == mResource.mMaterial->GetDimensions());
  glm::vec3 pos(inPos);
  glm::vec3 res(inRes);

  const glm::vec3 view_vector = mResource.mFroxelInfo->ComputeViewVector(inPos);
  const float depth_near = mResource.mFroxelInfo->SliceToDepth(inPos.z);
  const float depth_far = mResource.mFroxelInfo->SliceToDepth(inPos.z + 1);

  // depth
  float depth = 0.5f * (depth_far + depth_near);
  const glm::vec3 view_pos = mResource.mFroxelInfo->ComputeViewPosition(inPos.xy(), depth);
  const glm::vec3 world_pos = (mResource.mInvView * glm::vec4(view_pos, 1)).xyz();
  const glm::vec3 world_dir = glm::normalize(mResource.mInvView * glm::vec4(view_pos, 0)).xyz();

  // mie asymmetry factor
  const glm::vec3 g = 0.80f * glm::vec3(0.92f, 0.93f, 1.0f);

  // scattering coeff
  glm::vec4 coeffs = mResource.mMaterial->Load(inPos);
  // rayleigh
  glm::vec3 mu_sr = coeffs.w * util::RAYLEIGH_S;
  // +mie
  glm::vec3 mu_sm = coeffs.xyz();
  glm::vec3 mu_am = coeffs.xyz() * util::MIE_ABSORB_MUL;
  // total
  glm::vec3 mu_t = glm::max(mu_sr + mu_sm + mu_am, glm::vec3(1e-6f));

  glm::vec3 out_luminance(0);

  // point light (diffusion approx.)
#if 1
  for (const PointLight& light : mResource.mPointLights)
  {
    if (glm::length(mu_sm) < 0.001f)
    {
      break;
    }

    glm::vec3 to_light = light.position - world_pos;
    float distance_to_light = glm::length(to_light);

    // prevent zero div.
    distance_to_light = glm::max(distance_to_light, 1e-6f);

    // effective scattering coeff.
    glm::vec3 mu_s = mu_sm + mu_sr;
    glm::vec3 mu_a = mu_am;
    glm::vec3 mu_s_eff = mu_s * (1.0f - g);

    // diffusion coeff.
    glm::vec3 D = 1.0f / (3.0f * (mu_a + mu_s_eff));

    glm::vec3 l = light.luminance / (4.0f * glm::pi<glm::vec3::value_type>() * D * D);
    l *= glm::exp(-mu_s_eff * distance_to_light) / distance_to_light;
    l = glm::clamp(l, 0.0f, 4.0f); 
    out_luminance += l;
  }
#endif

  // directional light (single scattering)
  {
    const glm::vec3 sun_lum = 1.f * glm::vec3(1, 1, 1);
    const glm::vec3 sun_dir = glm::normalize(glm::vec3(0.35f, 1.0f, -0.25f));
    const glm::vec3 view_dir = world_dir;
    const float sun_dot_view = glm::dot(sun_dir, view_dir);
    const glm::vec3 phase_func_m = util::HenyeyGreenstein(sun_dot_view, g);

    // mie
    const glm::vec3 inscatter_m = sun_lum * phase_func_m * mu_sm / mu_t;

    // rayleigh
    const float phase_func_r = util::Rayleigh(sun_dot_view);
    const glm::vec3 inscatter_r = sun_lum * phase_func_r * mu_sr / mu_t;

    out_luminance += inscatter_m + inscatter_r;
  }

  mResource.mOutImage->Store(inPos, out_luminance);
}

