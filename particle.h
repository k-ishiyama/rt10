#pragma once
#include "util.h"

struct Particle
{
  glm::vec3 position;
  glm::vec3 velocity;
  glm::vec3 emission;
  float elapsedTime;
  bool hasBuoyancy;

  glm::vec3 GetEmission() const
  {
    // fadein (elapsed time)
    float s;
    {
      float fade_start = 0.0f;
      float fade_end = 0.07f;
      s = glm::clamp((elapsedTime - fade_start) / (fade_end - fade_start), 0.0f, 1.0f);
      s = s * s * (3.0f - 2.0f * s);
    }

    float y = position.y;

    // fadein (height)
    float fi;
    {
      float fade_start = -10.0f;
      float fade_end = -13.0f;
      fi = glm::clamp((y - fade_end) / (fade_start - fade_end), 0.0f, 1.0f);
      fi = fi * fi * (3.0f - 2.0f * fi);
    }

    // fadeout (height)
    float fo;
    {
      float fade_start = 4.0f;
      float fade_end = 7.0f;
      fo = glm::clamp((fade_end - y) / (fade_end - fade_start), 0.0f, 1.0f);
    }

    float t = s * fi * fo;
    t = t * t * (3.0f - 2.0f * t);
    return t * emission;
  }

  bool isDead() const
  {
    return (glm::length(GetEmission()) < 1e-6f) && (elapsedTime > 0.07f);;
  }
};

struct ParticleSystem
{
public:
  void addParticle(
    const glm::vec3& inInitialPos,
    const glm::vec3& inInitialVel,
    const glm::vec3& inEmission,
    const bool inBuoyancy)
  {
    mParticles.emplace_back(Particle{ inInitialPos, inInitialVel, inEmission, 0.0f, inBuoyancy });
  }
  
  void Update(float inDeltaTime)
  {
    const float dt = inDeltaTime;
    const float g = 9.81f;

    for (Particle& p : mParticles)
    {
      glm::vec3 accel(0.0f, - g, 0.0f);

      if (p.hasBuoyancy)
      {
        const glm::vec3 wind_force(0.0f, 0.0f, -3.0f);

        // random force
        glm::vec3 rf(mRandom.NextFloat01(), mRandom.NextFloat01(), mRandom.NextFloat01());
        rf = 2.0f * rf - glm::vec3(1.0f);
        
        // acceleration
        float buoyancy = 0.1f * g;
        accel = glm::vec3(7.0f * rf.x, 0.01f * rf.y + buoyancy, 7.0f * rf.z) + wind_force;
      }

      // velocity verlet
      p.position += p.velocity * dt + 0.5f * accel * dt * dt;
      p.velocity += accel * dt;

      p.elapsedTime += dt;
    }

    std::erase_if(mParticles, [](const Particle& p) { return p.isDead(); });
  }

  const std::vector<Particle>& GetParticles() const
  {
    return mParticles;
  }

private:
  std::vector<Particle> mParticles;
  util::XorShift64Star mRandom = util::XorShift64Star(234567);
};