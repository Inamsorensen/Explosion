#include "Particle.h"
#include "Emitter.h"
#include "Grid.h"

#include <ngl/VAOPrimitives.h>
#include <ngl/Transformation.h>
#include <ngl/ShaderLib.h>
#include <ngl/Mat4.h>
#include <ngl/Mat3.h>

Particle::Particle(ngl::Vec3 _position, ngl::Vec3 _velocity, float _mass, float _radius, float _lifeTime, float _initialTemperature, Emitter* _emitter)
{
  m_position=_position;
  m_velocity=_velocity;

  m_origin=_position;
  m_initVelocity=_velocity;

  m_lifeTime=_lifeTime;
  m_currLife=0.0;
  m_active=0;

  m_mass=_mass;
  m_radius=_radius;

  m_initialTemperature=_initialTemperature;
  m_temperature=_initialTemperature;

  m_burnState=Unburnt;
  m_sootAccumulator=0.0;

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
  ///To do: When reset particle, should temperature be initTemp or current temp in field at init position?

  m_currLife+=_dt;


//Reset particle when reaches outside bounding box

//  if (m_currLife>m_lifeTime || m_position.m_x<m_xMin || m_position.m_x>m_xMax || m_position.m_y<m_yMin || m_position.m_y>m_yMax || m_position.m_z<m_zMin || m_position.m_z>m_zMax)
//  {
//    m_position=m_origin;
//    m_velocity=m_initVelocity;
//    m_currLife=0.0;
//    m_active=0;
//    m_temperature=m_initialTemperature;

//  }
//  else
//  {
//    ///To do: Should velocity be updated before or after position
//    ///       Velocity should probably be drag force, not direct velocity
//    ///       Update temperature

//    //Update position and temperature
//    updatePosition(_dt);
//    updateTemperature(_dt);
//    burnParticle(_dt);
//  }

  updatePosition(_dt);
  updateTemperature(_dt);
  burnParticle(_dt);

}

void Particle::updatePosition(float _dt)
{
  //Calculate drag force

  Grid* grid=Grid::getGrid();
  //Find velocity of grid at particle position
  ngl::Vec3 velocityFromField=grid->getVelocityFromField(m_position);

//    std::cout<<"Position: ["<<m_position.m_x<<" "<<m_position.m_y<<" "<<m_position.m_z<<"] Velocity: ["<<velocityFromField.m_x<<" "<<velocityFromField.m_y<<" "<<velocityFromField.m_z<<"]\n";

  ngl::Vec3 velocityDifference=velocityFromField-m_velocity;
//    ngl::Vec3 velocityDiffSquared=velocityDifference.dot(velocityDifference);

  float dragArea=pow(m_radius,2);
  float density=1/m_mass;

  ngl::Vec3 dragForce=(m_emitter->getDragConstant()*density*dragArea)*velocityDifference;


  //Calculate new velocity using gravity and dragForce
  ngl::Vec3 gravity=ngl::Vec3(0.0,-9.81,0.0);
//  m_velocity+=_dt*gravity;
  m_velocity+=_dt*dragForce;

  //Calculate new position
  m_position+=_dt*m_velocity;
//    m_position+=_dt*velocityFromField;
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
//  m_temperature=gridTemp;
}

void Particle::burnParticle(float _dt)
{
  //Check if unburnt
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
    index0_X=trunc(indexParticle.m_x);
    index0_Y=trunc(indexParticle.m_y);
    index0_Z=trunc(indexParticle.m_z);
    int gridIndex=grid->getVectorIndex(index0_X, index0_Y, index0_Z);

    //Add temperature contribution
    float tempContribution=_dt*m_emitter->getHeatRelease()*burnRate;
    grid->addTemperatureFromParticle(gridIndex, tempContribution);

    //Add divergence contribution
    float divergenceContribution=_dt*m_emitter->getVolumeRelease()*burnRate;
    grid->addDivergenceFromParticle(gridIndex, divergenceContribution);

    //Burn off mass and accumulate soot
    m_mass-=_dt*m_emitter->getBurnRate();
    m_sootAccumulator+=_dt*m_emitter->getSootCreation()*burnRate;

    //Change to soot if soot accumulation is above threshold
    if (m_sootAccumulator>=m_emitter->getSootThreshold())
    {
      m_burnState=Soot;
    }
  }


}

void Particle::render() const
{
  //Get instance so can draw sphere
  ngl::VAOPrimitives *prim=ngl::VAOPrimitives::instance();

  //Set up shader with shader name from emitter
  ngl::ShaderLib *shader=ngl::ShaderLib::instance();
  shader->use(m_emitter->getShaderName());

  //Calculate rotated positions - in case camera is rotated
  ngl::Vec4 position;
  position.m_x=m_position.m_x;
  position.m_y=m_position.m_y;
  position.m_z=m_position.m_z;
  position.m_w=1;
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
