#include "Particle.h"
#include "Emitter.h"
#include "Grid.h"

#include <ngl/VAOPrimitives.h>
#include <ngl/Transformation.h>
#include <ngl/ShaderLib.h>
#include <ngl/Mat4.h>
#include <ngl/Mat3.h>

Particle::Particle(ngl::Vec3 _position, ngl::Vec3 _velocity, float _mass, float _lifeTime, Emitter* _emitter)
{
  m_mass=_mass;
  m_position=_position;
  m_velocity=_velocity;

  m_origin=_position;
  m_initVelocity=_velocity;

  m_lifeTime=_lifeTime;
  m_currLife=0.0;
  m_active=0;

  m_emitter=_emitter;

  //Set boundaries for particles
  Grid* grid=Grid::getGrid();
  ngl::Vec3 gridOrigin=grid->getPosition();
  float gridSize=grid->getGridSize();

  //Make sure particles cannot move into the outer cells in the grid.
  ///Could optionally not let the trilinear interpolation sample from cells that do not exist.
  m_xMin=gridOrigin.m_x;
  m_yMin=gridOrigin.m_y;
  m_zMin=gridOrigin.m_z;
  m_xMax=gridOrigin.m_x+gridSize;
  m_yMax=gridOrigin.m_y+gridSize;
  m_zMax=gridOrigin.m_z+gridSize;

}

void Particle::update(float _dt)
{
  ///To do: Take velocity from grid

  m_currLife+=_dt;


//Reset particle when reaches outside bounding box

  if (m_currLife>m_lifeTime || m_position.m_x<m_xMin || m_position.m_x>m_xMax || m_position.m_y<m_yMin || m_position.m_y>m_yMax || m_position.m_z<m_zMin || m_position.m_z>m_zMax)
  {
    m_position=m_origin;
    m_velocity=m_initVelocity;
    m_currLife=0.0;
    m_active=0;

  }
  else
  {
    ///To do: Should velocity be updated before or after position
    ///       Velocity should probably be drag force, not direct velocity
    Grid* grid=Grid::getGrid();

    ngl::Vec3 velocityFromField=grid->getVelocityFromField(m_position);

    m_position+=_dt*m_velocity;
    m_velocity+=velocityFromField;
  }
}


void Particle::render() const
{
  ngl::VAOPrimitives *prim=ngl::VAOPrimitives::instance();
  ngl::Transformation transformation;
  ngl::ShaderLib *shader=ngl::ShaderLib::instance();
  //Pass in shader name as parameter
  shader->use(m_emitter->getShaderName());
  transformation.setPosition(m_position);

  ngl::Mat4 MV;
  ngl::Mat4 MVP;
  ngl::Mat3 normalMatrix;
  ngl::Mat4 M;
  M=transformation.getMatrix();
  MV=M*m_emitter->getCamera()->getViewMatrix();
  MVP=MV*m_emitter->getCamera()->getProjectionMatrix();
  normalMatrix=MV;
  normalMatrix.inverse();

  shader->setShaderParamFromMat4("MV",MV);
  shader->setShaderParamFromMat4("MVP",MVP);
  shader->setShaderParamFromMat3("normalMatrix",normalMatrix);
  shader->setShaderParamFromMat4("M",M);

  prim->draw("sphere");
}
