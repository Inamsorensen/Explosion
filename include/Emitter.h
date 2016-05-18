#ifndef EMITTER
#define EMITTER

#include <ngl/Vec3.h>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <ngl/Camera.h>

#include "Particle.h"

class Emitter
{
public:
  //---------------------------------------------------------------------------------
  /// @brief Emitter ctor
  //---------------------------------------------------------------------------------
  Emitter(ngl::Vec3 _position, float _radius, int _noParticles, int _emissionRate);

  //---------------------------------------------------------------------------------
  /// @brief Emitter dtor
  //---------------------------------------------------------------------------------
  ~Emitter();


  //---------------------------------------------------------------------------------
  /// @brief Particle setup
  //---------------------------------------------------------------------------------
  void setUpParticles(float _particleMass, ngl::Vec3 _initialParticleVelocity, float _particleLifeTime, float _particleRadius, float _particleDragConstant, float _particleInitialTemperature);

  //---------------------------------------------------------------------------------
  /// @brief Update particles from emitter after time step _dt
  //---------------------------------------------------------------------------------
  void update(float _dt);

  //---------------------------------------------------------------------------------
  /// @brief Render particles from emitter
  //---------------------------------------------------------------------------------
  void renderParticles(ngl::Mat4 _ModelMatrix_Camera);


  //---------------------------------------------------------------------------------
  /// @brief Set and call parameters for rendering of particles
  //---------------------------------------------------------------------------------
  inline void setCamera(ngl::Camera* _camera){m_camera=_camera;}
  inline ngl::Camera* getCamera() const {return m_camera;}
  inline void setShaderName(const std::string &_shaderName){m_shaderName=_shaderName;}
  inline std::string getShaderName() const {return m_shaderName;}
  inline ngl::Mat4 getModelMatrixCamera() const {return m_ModelMatrix_Camera;}
  inline float getDragConstant() const {return m_particleDragConstant;}


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
  /// @brief Number of particles emitted each time step
  //---------------------------------------------------------------------------------
  int m_emissionRate;

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
  /// @brief Particle average lifetime
  //---------------------------------------------------------------------------------
  float m_particleLifeTime;
  //---------------------------------------------------------------------------------
  /// @brief Particle drag coefficient. Used to calculate drag force
  //---------------------------------------------------------------------------------
  float m_particleDragConstant;


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
};

#endif //EMITTER
