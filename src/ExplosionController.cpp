#include <ExplosionController.h>

ExplosionController* ExplosionController::m_instance=nullptr;


ExplosionController::ExplosionController()
{
  //Timer setup
  m_simTimeStep=0.01;
  m_simElapsedTime=0.0;
  m_noFrames=0;

  //Grid setup
  m_gridPosition=ngl::Vec3(-2.5,-2.5,-2.5);
  m_gridSize=5.0;
  m_noCells=16;

  m_noiseConstant=0.2;
  m_vorticityConstant=4.0;

  m_ambientTemp=293.0; //20C=293K
  m_thermalConductivity=0.01;
  m_coolingConstant=1550.0;
  m_buoyancyConstant=0.003;

  m_grid=Grid::createGrid(m_gridPosition, m_gridSize, m_noCells);
  m_grid->setVelocityFieldConstants(m_noiseConstant, m_vorticityConstant);
  m_grid->setTemperatureFieldConstants(m_ambientTemp, m_thermalConductivity, m_coolingConstant, m_buoyancyConstant);

  //Explosion setup
  m_explosionOrigin=ngl::Vec3(0.0,0.0,0.0);
//  m_explosionRadius=0.3;
  m_explosionRadius=1.0;
  m_explosionTemperature=10000;
//  m_explosionDivergence=0.01;
  m_explosionDivergence=10.0;

  m_grid->setExplosion(m_explosionOrigin, m_explosionRadius, m_explosionTemperature, m_explosionDivergence);

  //Particle setup
  m_particleMass=0.1;
  m_particleRadius=0.1;
  m_particleInitVelocity=ngl::Vec3(0.0,0.0,0.0);
  m_particleLifeTime=10;
  m_particleDragConstant=950.0;
  m_particleInitTemperature=m_ambientTemp;

  //Emitter setup
  m_emitterPosition=m_explosionOrigin;
  m_emitterRadius=m_explosionRadius;
  m_noParticles=2000;
//  m_noParticles=1;
  m_emissionRate=10;

  m_emitter=new Emitter(m_emitterPosition, m_emitterRadius, m_noParticles, m_emissionRate);
  m_emitter->setUpParticles(m_particleMass, m_particleInitVelocity, m_particleLifeTime, m_particleRadius, m_particleDragConstant, m_particleInitTemperature);

}

ExplosionController* ExplosionController::instance()
{
  if (m_instance==nullptr)
  {
    m_instance=new ExplosionController();
  }

  return m_instance;
}

void ExplosionController::setRenderVariables(ngl::Camera *_camera, std::string _shaderName)
{
  m_camera=_camera;
  m_shaderName=_shaderName;

  m_emitter->setCamera(m_camera);
  m_emitter->setShaderName(m_shaderName);
}

void ExplosionController::update()
{
  //Increase elapsed time
  m_simElapsedTime+=m_simTimeStep;

  //Update grid values
  m_grid->update(m_simTimeStep);

  //Update particle values
  m_emitter->update(m_simTimeStep);

}

void ExplosionController::render(ngl::Mat4 _ModelMatrix_Camera)
{
  m_noFrames+=1;
  m_emitter->renderParticles(_ModelMatrix_Camera);
}
