#ifndef PARTICLE
#define PARTICLE

#include <ngl/Vec3.h>
#include <ngl/Mat4.h>
#include "Grid.h"
#include "mathFunction.h"

/// @brief Particle.h
/// Single particle in simulation
/// Contains position, velocity, temperature and burn state data
/// Author: Ina M. Sorensen
/// Date: 27.05.16
/// Implementation based on the code by
/// Christian Miller (2007). Realtime Explosion Simulation [online].
/// [Accessed May 2016]. Available from: <http://www.cs.utexas.edu/~ckm/explosion/>.

enum particleStates
{
  Unburnt,
  Burning,
  Soot
};

class Emitter;

class Particle
{
public:
  //---------------------------------------------------------------------------------
  /// @brief Particle ctor
  //---------------------------------------------------------------------------------
  Particle(ngl::Vec3 _position, ngl::Vec3 _velocity, float _mass, float _radius, float _initialTemperature, Emitter *_emitter);

  //---------------------------------------------------------------------------------
  /// @brief Update particle position and velocity for time step _dt
  //---------------------------------------------------------------------------------
  void update(float _dt);

  //---------------------------------------------------------------------------------
  /// @brief Render out the particle
  //---------------------------------------------------------------------------------
  void render() const;

  //---------------------------------------------------------------------------------
  /// @brief Get particle position
  //---------------------------------------------------------------------------------
  inline ngl::Vec3 getPosition(){return m_position;}

  //---------------------------------------------------------------------------------
  /// @brief Is particle active, if it is will be rendered
  //---------------------------------------------------------------------------------
  bool m_active;

private:
  //---------------------------------------------------------------------------------
  /// @brief Update particle position from gravity and gas field velocity
  //---------------------------------------------------------------------------------
  void updatePosition(float _dt);
  //---------------------------------------------------------------------------------
  /// @brief Update particle temperature from gas field
  //---------------------------------------------------------------------------------
  void updateTemperature(float _dt);
  //---------------------------------------------------------------------------------
  /// @brief Burn particle
  //---------------------------------------------------------------------------------
  void burnParticle(float _dt);

  //---------------------------------------------------------------------------------
  /// @brief Properties of one particle
  //---------------------------------------------------------------------------------
  ngl::Vec3 m_position;
  ngl::Vec3 m_velocity;

  //---------------------------------------------------------------------------------
  /// @brief Size and mass of particle
  //---------------------------------------------------------------------------------
  float m_radius;
  float m_mass;

  //---------------------------------------------------------------------------------
  /// @brief Temperature of particle
  //---------------------------------------------------------------------------------
  float m_temperature;

  //---------------------------------------------------------------------------------
  /// @brief Particle state and soot accumulator
  //---------------------------------------------------------------------------------
  particleStates m_burnState;
  float m_sootAccumulator;

  //---------------------------------------------------------------------------------
  /// @brief Emitter that particle belongs to
  //---------------------------------------------------------------------------------
  const Emitter* m_emitter;

  //---------------------------------------------------------------------------------
  /// @brief Boundaries for particles, ie. size of bounding box
  //---------------------------------------------------------------------------------
  float m_xMin;
  float m_xMax;
  float m_yMin;
  float m_yMax;
  float m_zMin;
  float m_zMax;

};

#endif //PARTICLE
