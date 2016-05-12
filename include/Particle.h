#ifndef PARTICLE
#define PARTICLE

#include <ngl/Vec3.h>
#include "Grid.h"

class Emitter;

class Particle
{
public:
  //---------------------------------------------------------------------------------
  /// @brief Particle ctor
  //---------------------------------------------------------------------------------
  Particle(ngl::Vec3 _position, ngl::Vec3 _velocity, float _mass, float _lifeTime, Emitter *_emitter);

  //---------------------------------------------------------------------------------
  /// @brief Update particle position and velocity for time step _dt
  //---------------------------------------------------------------------------------
  void update(float _dt);

  //---------------------------------------------------------------------------------
  /// @brief Render out the particle
  //---------------------------------------------------------------------------------
  void render() const;

   bool m_active;

private:
  //---------------------------------------------------------------------------------
  /// @brief Properties of one particle
  //---------------------------------------------------------------------------------
  float m_mass;
  ngl::Vec3 m_position;
  ngl::Vec3 m_velocity;

  ngl::Vec3 m_origin;
  ngl::Vec3 m_initVelocity;

  //---------------------------------------------------------------------------------
  /// @brief Particle lifetime parameters
  //---------------------------------------------------------------------------------
  float m_lifeTime;
  float m_currLife;

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
