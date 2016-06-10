#ifndef EMITTER
#define EMITTER

#include <ngl/Vec3.h>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <ngl/Camera.h>

#include "Particle.h"
#include "AlembicExport.h"

/// @brief Emitter.h
/// Emitter for particle system
/// Sets up position and initial values for particles and contains most particle parameters
/// Controls update and rendering of particles
/// Author: Ina M. Sorensen
/// Date: 27.05.16
/// Implementation based on the code by
/// Christian Miller (2007). Realtime Explosion Simulation [online].
/// [Accessed May 2016]. Available from: <http://www.cs.utexas.edu/~ckm/explosion/>.
/// Jon Macey (2015). Simple Particles [online]
/// [Accessed May 2016] Available from: <https://github.com/NCCA/ParticleSystem>

/// @brief Emitter.h
/// Emitter for particle system
/// Sets up position and initial values for particles and contains most particle parameters
/// Controls update and rendering of particles
/// Author: Ina M. Sorensen
/// Date: 27.05.16
/// Implementation based on the code by
/// Christian Miller (2007). Realtime Explosion Simulation [online].
/// [Accessed May 2016]. Available from: <http://www.cs.utexas.edu/~ckm/explosion/>.
/// Jon Macey (2015). Simple Particles [online]
/// [Accessed May 2016] Available from: <https://github.com/NCCA/ParticleSystem>

class Emitter
{
public:
  //---------------------------------------------------------------------------------
  /// @brief Emitter ctor
  //---------------------------------------------------------------------------------
  Emitter(ngl::Vec3 _position, float _radius, int _noParticles);

  //---------------------------------------------------------------------------------
  /// @brief Emitter dtor
  //---------------------------------------------------------------------------------
  ~Emitter();


  //---------------------------------------------------------------------------------
  /// @brief Particle setup
  //---------------------------------------------------------------------------------
  void setUpParticles(float _particleMass, float _massRandomisation, ngl::Vec3 _initialParticleVelocity, float _particleRadius, float _particleDragConstant, float _particleInitialTemperature);

  //---------------------------------------------------------------------------------
  /// @brief Update particles from emitter after time step _dt
  //---------------------------------------------------------------------------------
  void update(float _dt);

  //---------------------------------------------------------------------------------
  /// @brief Render particles from emitter
  //---------------------------------------------------------------------------------
  void renderParticles(ngl::Mat4 _ModelMatrix_Camera, bool _exportAlembic);


  //---------------------------------------------------------------------------------
  /// @brief Set and call parameters for rendering of particles
  //---------------------------------------------------------------------------------
  inline void setCamera(ngl::Camera* _camera){m_camera=_camera;}
  inline ngl::Camera* getCamera() const {return m_camera;}
  inline void setShaderName(const std::string &_shaderName){m_shaderName=_shaderName;}
  inline std::string getShaderName() const {return m_shaderName;}
  inline ngl::Mat4 getModelMatrixCamera() const {return m_ModelMatrix_Camera;}
  inline float getDragConstant() const {return m_particleDragConstant;}

  //---------------------------------------------------------------------------------
  /// @brief Set and call parameters for burning of particles
  //---------------------------------------------------------------------------------
  void setParticleBurningParameters(float _burnThreshold, float _burnRate, float _thermalConductivity, float _thermalMass, float _heatRelease, float _volumeRelease, float _sootThreshold);
  inline float getBurnThreshold() const {return m_burnThreshold;}
  inline float getBurnRate() const {return m_burnRate;}
  inline float getThermalConductivity() const {return m_thermalConductivity;}
  inline float getThermalMass() const {return m_thermalMass;}
  inline float getHeatRelease() const {return m_heatRelease;}
  inline float getVolumeRelease() const {return m_volumeRelease;}
  inline float getSootThreshold() const {return m_sootThreshold;}



private:

  //---------------------------------------------------------------------------------
  /// @brief Position of the centre of the emitter
  //---------------------------------------------------------------------------------
  ngl::Vec3 m_position;
  //---------------------------------------------------------------------------------
  /// @brief Radius of disk emitter
  //---------------------------------------------------------------------------------
  float m_radius;
  //---------------------------------------------------------------------------------
  /// @brief Number of particles total
  //---------------------------------------------------------------------------------
  int m_noParticles;

  //---------------------------------------------------------------------------------
  /// @brief Particle mass
  //---------------------------------------------------------------------------------
  float m_particleMass;
  //---------------------------------------------------------------------------------
  /// @brief Particle radius
  //---------------------------------------------------------------------------------
  float m_particleRadius;
  //---------------------------------------------------------------------------------
  /// @brief Particle initial velocity
  //---------------------------------------------------------------------------------
  ngl::Vec3 m_particleInitVelocity;
  //---------------------------------------------------------------------------------
  /// @brief Particle drag coefficient. Used to calculate drag force
  //---------------------------------------------------------------------------------
  float m_particleDragConstant;

  //---------------------------------------------------------------------------------
  /// @brief Parameters for burning of particles
  //---------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------
  /// @brief Threshold temperature at which fuel particles start burning
  //---------------------------------------------------------------------------------
  float m_burnThreshold;
  //---------------------------------------------------------------------------------
  /// @brief Thermal conductivity between air and fuel
  //---------------------------------------------------------------------------------
  float m_thermalConductivity;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Rate at which fuel particles burn, ie. mass combusted per second
  //----------------------------------------------------------------------------------------------------------------------
  float m_burnRate;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Thermal mass of particle
  //----------------------------------------------------------------------------------------------------------------------
  float m_thermalMass;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Heat released per unit combusted mass of fuel
  //----------------------------------------------------------------------------------------------------------------------
  float m_heatRelease;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Gas volume released per unit combusted mass of fuel
  //----------------------------------------------------------------------------------------------------------------------
  float m_volumeRelease;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Mass threshold below which a fuel particle turns to soot
  //----------------------------------------------------------------------------------------------------------------------
  float m_sootThreshold;


  //---------------------------------------------------------------------------------
  /// @brief List of particles that have come from this emitter
  //---------------------------------------------------------------------------------
  std::vector <Particle*> m_particles;


  //---------------------------------------------------------------------------------
  /// @brief Camera
  //---------------------------------------------------------------------------------
  ngl::Camera* m_camera;
  //---------------------------------------------------------------------------------
  /// @brief Camera model matrix, based on camera rotation etc
  //---------------------------------------------------------------------------------
  ngl::Mat4 m_ModelMatrix_Camera;
  //---------------------------------------------------------------------------------
  /// @brief Name of particle shader
  //---------------------------------------------------------------------------------
  std::string m_shaderName;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Alembic exporter
  //----------------------------------------------------------------------------------------------------------------------
  AlembicExport* m_alembicExporter;
};

#endif //EMITTER
