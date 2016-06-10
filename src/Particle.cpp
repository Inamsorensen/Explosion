#include "Particle.h"
#include "Emitter.h"
#include "Grid.h"

#include <ngl/VAOPrimitives.h>
#include <ngl/Transformation.h>
#include <ngl/ShaderLib.h>
#include <ngl/Mat4.h>
#include <ngl/Mat3.h>

Particle::Particle(ngl::Vec3 _position, ngl::Vec3 _velocity, float _mass, float _radius, float _initialTemperature, Emitter* _emitter)
{
  m_position=_position;
  m_velocity=_velocity;

  m_active=true;

  m_mass=_mass;
  m_radius=_radius;

  m_temperature=_initialTemperature;

  m_burnState=Unburnt;
  m_sootAccumulator=0.0;

  m_emitter=_emitter;

  //Set boundaries for particles
  Grid* grid=Grid::getGrid();
  ngl::Vec3 gridOrigin=grid->getPosition();
  float gridSize=grid->getGridSize();
  m_xMin=gridOrigin.m_x;
  m_yMin=gridOrigin.m_y;
  m_zMin=gridOrigin.m_z;
  m_xMax=gridOrigin.m_x+gridSize;
  m_yMax=gridOrigin.m_y+gridSize;
  m_zMax=gridOrigin.m_z+gridSize;

}

void Particle::update(float _dt)
{

//Stop updating and rendering particle when reaches outside bounding box

  if (m_position.m_x<m_xMin || m_position.m_x>m_xMax || m_position.m_y<m_yMin || m_position.m_y>m_yMax || m_position.m_z<m_zMin || m_position.m_z>m_zMax)
  {
    m_active=false;

  }
  else
  {

    //Update position and temperature, and burn particle
    updatePosition(_dt);
    updateTemperature(_dt);
    burnParticle(_dt);
  }

}

void Particle::updatePosition(float _dt)
{
  Grid* grid=Grid::getGrid();

    //Find velocity of grid at particle position
  ngl::Vec3 velocityFromField=grid->getVelocityFromField(m_position);

  //Only use drag force on fuel particles, not on soot
  if (m_burnState!=Soot)
  {
    ngl::Vec3 velocityDifference=velocityFromField-m_velocity;

    float dragArea=pow(m_radius,2);
    float density=1.0/m_mass;

    ngl::Vec3 dragForce=(m_emitter->getDragConstant()*density*dragArea)*velocityDifference*(velocityDifference.length());
//    ngl::Vec3 dragForce=(m_emitter->getDragConstant()*density*dragArea)*velocityDifference;

    //If statement inserted to make sure the particles don't speed up a lot when mass decreases
    if (std::abs(dragForce.m_x)<std::abs(velocityDifference.m_x) && std::abs(dragForce.m_y)<std::abs(velocityDifference.m_y) && std::abs(dragForce.m_z)<std::abs(velocityDifference.m_z))
    {
      m_velocity+=_dt*dragForce;
    }
    else
    {
      m_velocity=velocityFromField;
    }
  }
  else
  {
    m_velocity=velocityFromField;
  }

  //Calculate new position
  m_position+=_dt*m_velocity;

}

void Particle::updateTemperature(float _dt)
{
  //Get temperature of air at particle position
  Grid* grid=Grid::getGrid();
  float gridTemp=grid->getTemperatureFromField(m_position);

  //Use to calculate temperature change
  float tempChange=(m_emitter->getThermalConductivity()/m_emitter->getThermalMass())*(gridTemp-m_temperature);

  //Calculate new temperature
  m_temperature+=_dt*tempChange;
}

void Particle::burnParticle(float _dt)
{
  //If unburnt, check temperature to see if ignited
  if (m_burnState==Unburnt)
  {
    if (m_temperature>=m_emitter->getBurnThreshold())
    {
      m_burnState=Burning;
    }
  }
  else if (m_burnState==Burning)
  {
    //Get burn rate
    float burnRate=m_emitter->getBurnRate();

    //Get position of particle in field
    Grid* grid=Grid::getGrid();
    ngl::Vec3 indexParticle=(1/grid->getCellSize())*(m_position-grid->getPosition());
    int index0_X;
    int index0_Y;
    int index0_Z;
    index0_X=floor(indexParticle.m_x);
    index0_Y=floor(indexParticle.m_y);
    index0_Z=floor(indexParticle.m_z);
    int gridIndex=grid->getVectorIndex(index0_X, index0_Y, index0_Z);

    //Add temperature contribution
    float tempContribution=_dt*m_emitter->getHeatRelease()*burnRate;
    grid->addTemperatureFromParticle(gridIndex, tempContribution);

    //Add divergence contribution
    float divergenceContribution=_dt*m_emitter->getVolumeRelease()*burnRate;
    grid->addDivergenceFromParticle(gridIndex, divergenceContribution);

    //Burn off mass
    m_mass-=_dt*m_emitter->getBurnRate();

    //Change to soot if remaining mass is above threshold
    if (m_mass<m_emitter->getSootThreshold())
    {
      m_burnState=Soot;
    }

  }


}

void Particle::render() const
{
  //Get instance so can draw sphere
  ngl::VAOPrimitives *prim=ngl::VAOPrimitives::instance();

  //Set shader depending on burn state
  ngl::ShaderLib *shader=ngl::ShaderLib::instance();

  if (m_burnState==Soot)
  {
    shader->use("PhongSilver");
  }
  else if (m_burnState==Unburnt)
  {
    shader->use("PhongCopper");
  }
  else
  {
    shader->use(m_emitter->getShaderName());
  }


  //Calculate rotated positions - in case camera is rotated
  ngl::Vec4 position;
  position.m_x=m_position.m_x;
  position.m_y=m_position.m_y;
  position.m_z=m_position.m_z;
  position.m_w=1.0;
  ngl::Vec4 rotatedPosition=position*m_emitter->getModelMatrixCamera();

  //Set MVP matrix
  ngl::Transformation transformation;
  transformation.setPosition(rotatedPosition.m_x, rotatedPosition.m_y, rotatedPosition.m_z);
  transformation.setScale(0.05,0.05,0.05);

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
