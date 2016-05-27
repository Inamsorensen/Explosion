#include <ExplosionController.h>

ExplosionController* ExplosionController::m_instance=nullptr;


ExplosionController::ExplosionController()
{
  //Timer setup
  m_simTimeStep=0.025;
  m_simElapsedTime=0.0;
  m_noFrames=0;

  //Grid setup
  m_gridPosition=ngl::Vec3(-0.5,-0.5,-0.5);
  m_gridSize=1.0;
  m_noCells=16;
//  m_noCells=32;

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
  m_explosionOrigin=ngl::Vec3(0.0,-0.3,0.0);
  m_explosionRadius=0.1;
  m_explosionTemperature=4000;
  m_explosionDivergence=3.0;

  m_grid->setExplosion(m_explosionOrigin, m_explosionRadius, m_explosionTemperature, m_explosionDivergence);


  //Particle setup
  m_particleMass=0.7;
  m_particleMassRandomisation=2.0;
  m_particleRadius=0.1;
  m_particleInitVelocity=ngl::Vec3(0.0,0.0,0.0);
  m_particleDragConstant=950.0;
  m_particleInitTemperature=m_ambientTemp;


  //Emitter setup
  m_emitterPosition=m_explosionOrigin; ///Increasing this so definitely inside sphere would be good too
  m_emitterRadius=m_explosionRadius; ///Decreasing this below the explosion radius will decrease the number of particles going sideways
  m_noParticles=16000;

  m_emitter=new Emitter(m_emitterPosition, m_emitterRadius, m_noParticles);
  m_emitter->setUpParticles(m_particleMass, m_particleMassRandomisation, m_particleInitVelocity, m_particleRadius, m_particleDragConstant, m_particleInitTemperature);


  //Fuel particle burn set up
  m_particleBurnThreshold=500.0;
  m_particleBurnRate=0.67;
  m_particleHeatRelease=40.0*(65536.0/(float)m_noParticles);
  m_particleVolumeRelease=0.001*(65536.0/(float)m_noParticles);
  m_particleThermalMass=1.6;
  m_particleThermalConductivity=1.0;
  m_particleSootThreshold=0.01;

  m_emitter->setParticleBurningParameters(m_particleBurnThreshold, m_particleBurnRate, m_particleThermalConductivity, m_particleThermalMass, m_particleHeatRelease, m_particleVolumeRelease, m_particleSootThreshold);

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
