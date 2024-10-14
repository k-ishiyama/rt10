#include "dds.h"
#include "util.h"
#include "final_compose.h"
#include "froxel_participating_media.h"
#include "froxel_raymarching.h"
#include "froxel_luminance.h"
#include "timer.h"
#include "particle.h"
#include "logger.h"

int main(int argc, char** argv)
{
  Timer timer;

  namespace fs = std::filesystem;
  std::unique_ptr<util::Texture3D<uint8_t, util::WrapPolicy3D>> perlin_worley_erosion
    = util::CreateTextureFromDDS<uint8_t, util::WrapPolicy3D>(fs::current_path() / "etc" / "PerlinWorleyErosion.dds");
  std::unique_ptr<util::Texture3D<uint8_t, util::WrapPolicy3D>> perlin_worley_shape
    = util::CreateTextureFromDDS<uint8_t, util::WrapPolicy3D>(fs::current_path() / "etc" / "PerlinWorleyShape.dds");

  LightResource light_resource;
  RaymarchingResource rm_resource;
  FinalComposeResource resource;
  ParticipatingMediaResource pm_resource;
  pm_resource.mPerlinWorleyShape = perlin_worley_shape.get();
  pm_resource.mPerlinWorleyErosion = perlin_worley_erosion.get();

  constexpr float froxel_near = 1.0f;
  constexpr float froxel_far = 70.0f;
#if 1
  constexpr uint32_t film_size_x = 5 * 128;
  constexpr uint32_t film_size_y = 5 * 72;
  constexpr uint32_t froxel_size_x = film_size_x / 2;
  constexpr uint32_t froxel_size_y = film_size_y / 2;
  constexpr uint32_t froxel_size_z = 96;
  const int cMaxFrameNum = 300;
  const float fps = 30.0f; // [1/s]
  const float time_limit = 6.0f * 60.0f * 60.0f; // [s]
#else
  constexpr uint32_t film_size_x = 8 * 128;
  constexpr uint32_t film_size_y = 8 * 72;
  constexpr uint32_t froxel_size_x = film_size_x / 2;
  constexpr uint32_t froxel_size_y = film_size_y / 2;
  constexpr uint32_t froxel_size_z = 96;
  const int cMaxFrameNum = 400;
  const float fps = 40.0f; // [1/s]
  const float time_limit = 256.0f; // [s]
#endif
  const glm::uvec2 film_size = glm::uvec2(film_size_x, film_size_y);
  const glm::uvec3 froxel_size(froxel_size_x, froxel_size_y, froxel_size_z);
  const float cMaxTime = cMaxFrameNum / fps; // [s]
  const float dt = cMaxTime / cMaxFrameNum;
  float simulation_time = 0.0f;

  Log("omp_get_num_procs   = {}", omp_get_num_procs());
  Log("omp_get_max_threads = {}", omp_get_max_threads());
  Log("total time = {:.2f}[s], max frames = {}, dt = {:.4f}(FPS = {:.2f})", cMaxTime, cMaxFrameNum, dt, 1.0f/dt);
  Log("film size = {} x {}, froxel size = {} x {} x {}, near = {:.2f}, far = {:.2f}",
    film_size.x, film_size.y, froxel_size.x, froxel_size.y, froxel_size.z, froxel_near, froxel_far);
  Log("-------");

  ParticleSystem particle_system;
  particle_system.addParticle({ -60, -4, 0 }, { 50, 10, 0 }, { 200, 200, 200 }, false);

  util::XorShift64Star random(123456);
  float particle_timer = 0.0f;
  float particle_timer2 = 2.0f;
  float last_frame_time = 0.0f;
  for (int i = 0; i < cMaxFrameNum; ++i)
  {
    float start_time = timer.GetElapsedTimeSecond();

    if ((timer.GetElapsedTimeSecond() + 1.2f * last_frame_time) > time_limit)
    {
      break;
    }

    const float progress = float(i) / (cMaxFrameNum - 1);
    util::Film film(film_size);

    // frustum, camera
    glm::mat4x4 view_mat;
    float fov = 70.0f * util::PI<float> / 180.0f;
    float aspect_ratio = static_cast<float>(film.GetWidth()) / film.GetHeight();
    float near = froxel_near;
    float far = froxel_far;
    util::Frustum frustum(fov, aspect_ratio, near, far);
    float z_offset = 10.f * simulation_time;
    glm::vec3 camera_pos(-30.0f + z_offset, 0.5f * (1.0f - progress) + 6.1f * progress, 0.0f);
    glm::vec3 look_at(0.0f + z_offset, 0.5f, 0.0f);
    glm::vec3 up_dir(0.0f, 1.0f, 0.0f);
    view_mat = glm::lookAt(camera_pos, look_at, up_dir);
    const glm::mat4x4 inv_view = glm::inverse(view_mat);

    // particles
    particle_timer -= dt;
    if (particle_timer <= 0.0f)
    {
      util::Plane plane(glm::vec3(40.0f + z_offset, -12.9f, 15.0f), glm::vec2(50.0f, 30.0f));

      glm::vec3 pos
        = util::GenerateRandomPointOnPlane({ random.NextFloat01(), random.NextFloat01() }, plane);
      glm::vec3 random_vector(random.NextFloat01(), random.NextFloat01(), random.NextFloat01());
      glm::vec3 color = 30.0f * (glm::vec3(0.2f) + 0.8f * random_vector);

      float rx = 2.0f * random.NextFloat01() - 1.0f;
      float ry = random.NextFloat01();
      float rz = 2.0f * random.NextFloat01() - 1.0f;
      particle_system.addParticle(pos, { 0.01f * rx, 2.2f + 0.5f * ry, 0.01f * rz }, color, true);

      particle_timer = 0.025f;
    }

    particle_timer2 -= dt;
    if (particle_timer2 <= 0.0f)
    {
      glm::vec3 pos = glm::vec3(camera_pos.x - 30.0f, -4.0f, 2.0f * random.NextFloat01() - 1.0f);
      glm::vec3 random_vector(random.NextFloat01(), random.NextFloat01(), random.NextFloat01());
      glm::vec3 color = 200.0f * (glm::vec3(0.8f) + 0.2f * random_vector);

      float rx = random.NextFloat01();
      float ry = random.NextFloat01();
      float rz = 2.0f * random.NextFloat01() - 1.0f;
      particle_system.addParticle(pos, { 70.0f + 10.0f * rx, 10.0f * ry + 12.0f * (1.0f - ry), rz }, color, false);

      particle_timer2 = 1.91f;
    }
    particle_system.Update(dt);

#if 0
    if (i != 184)
    // if (i < 100)
    {
      // skip
      simulation_time += dt;
      continue;
    }
#endif

    // point lights
    light_resource.mPointLights.clear();
    light_resource.mPointLights.reserve(particle_system.GetParticles().size());
    for (const Particle& p : particle_system.GetParticles())
    {
      light_resource.mPointLights.emplace_back(PointLight{ p.position, p.GetEmission() });
#if 0    
      float dist = glm::length(camera_pos - p.position);
      Log("{}: pos = ({:.3f}, {:.3f}, {:.3f}), vel = ({:.3f}, {:.3f}, {:.3f}), distance = {:.3f}, emission = ({:.2f}, {:.2f}, {:.2f})",
        i, p.position.x, p.position.y, p.position.z,
        p.velocity.x, p.velocity.y, p.velocity.z,
        dist, p.emission.x, p.emission.y, p.emission.z);
#endif
    }

    resource.mInvProj = frustum.GetProjView();
    resource.mInvView = inv_view;

    util::FroxelInfo froxel_info(frustum.GetProjView(), froxel_near, froxel_far, froxel_size);
    util::Texture3D<glm::vec4, util::ClampPolicy3D> pm_volume(froxel_size);
    util::Texture3D<glm::vec3, util::ClampPolicy3D> light_volume(froxel_size);
    util::Image3D<glm::vec3> rm_luminance(froxel_size);
    util::Image3D<glm::vec3> rm_transmittance(froxel_size);
    util::Image2D<glm::vec3> render_target({ film.GetWidth(), film.GetHeight() });

    // froxel : participating media
    pm_resource.mFroxelInfo = &froxel_info;
    pm_resource.mInvView = inv_view;
    pm_resource.mElapsedTime = simulation_time;
    pm_resource.mOutImage = &pm_volume;
    util::Compute<ParticipatingMediaOperator, ParticipatingMediaResource, froxel_size_x, froxel_size_y, froxel_size_z>
      froxel_media(pm_resource);
    froxel_media.Exec();

    // froxel : lighting
    light_resource.mMaterial = &pm_volume;
    light_resource.mFroxelInfo = &froxel_info;
    light_resource.mInvView = inv_view;
    light_resource.mOutImage = &light_volume;
    util::Compute<LightOperator, LightResource, froxel_size_x, froxel_size_y, froxel_size_z>
      froxel_light(light_resource);
    froxel_light.Exec();

    // froxel : ray marching
    rm_resource.mOutLuminance = &rm_luminance;
    rm_resource.mOutTransmittance = &rm_transmittance;
    rm_resource.mMaterial = &pm_volume;
    rm_resource.mLight = &light_volume;
    rm_resource.mInvView = inv_view;
    rm_resource.mFroxelInfo = &froxel_info;
    util::Compute<RaymarchingOperator, RaymarchingResource, froxel_size_x, froxel_size_y, 1>
      froxel_raymarching(rm_resource);
    froxel_raymarching.Exec();

    // final compose
    util::Texture3D<glm::vec3, util::ClampPolicy3D> luminance_tex = std::move(rm_luminance);
    util::Texture3D<glm::vec3, util::ClampPolicy3D> transmitt_tex = std::move(rm_transmittance);
    resource.mInvProj = frustum.GetProjView();
    resource.mInvView = glm::inverse(view_mat);
    resource.mFroxelFar = froxel_far;
    resource.mFroxelNear = froxel_near;
    resource.mLuminance = &luminance_tex;
    resource.mTransmittance = &transmitt_tex;
    resource.mOutImage = &render_target;
    util::Compute<FinalComposeOperator, FinalComposeResource, film_size_x, film_size_y, 1>
      renderer(resource);
    renderer.Exec();

    // export
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(3) << i;
    std::string filename = ss.str() + ".png";
    film.Develop(render_target);
    util::Export(filename, film.GetImage());

#if 0
    // debug
    if (i == 0)
    {
      namespace fs = std::filesystem;
      SaveDDS(
        fs::current_path() / "participating_media.dds", DDS_FORMAT_R32G32B32A32_FLOAT, sizeof(glm::vec4),
        pm_volume.GetWidth(), pm_volume.GetHeight(), pm_volume.GetDepth(), pm_volume.GetData());
      SaveDDS(
        fs::current_path() / "light.dds", DDS_FORMAT_R32G32B32_FLOAT, sizeof(glm::vec3),
        light_volume.GetWidth(), light_volume.GetHeight(), light_volume.GetDepth(), light_volume.GetData());
      SaveDDS(
        fs::current_path() / "raymarching_luminance.dds", DDS_FORMAT_R32G32B32_FLOAT, sizeof(glm::vec3),
        luminance_tex.GetWidth(), luminance_tex.GetHeight(), luminance_tex.GetDepth(), luminance_tex.GetData());
      SaveDDS(
        fs::current_path() / "raymarching_transmittance.dds", DDS_FORMAT_R32G32B32_FLOAT, sizeof(glm::vec3),
        transmitt_tex.GetWidth(), transmitt_tex.GetHeight(), transmitt_tex.GetDepth(), transmitt_tex.GetData());

    }
#endif

    last_frame_time = timer.GetElapsedTimeSecond() - start_time;
    Log("[{:03d}] {:.3f}[s] : {}, simulation time = {:.3f}[s], frame time = {:.3f}[s], particles = {}",
      i, timer.GetElapsedTimeSecond(), filename, simulation_time, last_frame_time, particle_system.GetParticles().size());

    simulation_time += dt;
}
  Log("done.\n");
  return 0;
}