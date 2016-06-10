#include "Emitter.h"

#include <ngl/Random.h>
#include <cmath>
#include <math.h>

Emitter::Emitter(ngl::Vec3 _position, float _radius, int _noParticles)
{
  m_position=_position;
  m_radius=_radius;

  m_noParticles=_noParticles;

  m_alembicExporter=new AlembicExport();



}

Emitter::~Emitter()
{
  int noParticlesCurr=m_particles.size();

  if(noParticlesCurr!=0)
  {
    for (int i=0; i<noParticlesCurr; i++)
    {
      delete m_particles[i];
    }
  }
  m_particles.clear();
}

void Emitter::setUpParticles(float _particleMass, float _massRandomisation, ngl::Vec3 _initialParticleVelocity, float _particleRadius, float _particleDragConstant, float _particleInitialTemperature)
{
  m_particleMass=_particleMass;
  m_particleInitVelocity=_initialParticleVelocity;
  m_particleRadius=_particleRadius;
  m_particleDragConstant=_particleDragConstant;

  for (int i=0; i<m_noParticles; i++)
  {
    //Create random position within sphere
    ngl::Random* rand=ngl::Random::instance();
    rand->setSeed(i);
    float randRadius=rand->randomPositiveNumber(m_radius);
    float randTheta=rand->randomPositiveNumber(2*M_PI);
    float randPhi=rand->randomPositiveNumber(M_PI);
    ngl::Vec3 pos=m_position;
    pos.m_x+=randRadius*cos(randTheta)*sin(randPhi);
    pos.m_y+=randRadius*cos(randPhi);
    pos.m_z+=randRadius*sin(randTheta)*sin(randPhi);

    //Randomise mass slightly
    float massRand=rand->randomNumber(m_particleMass);
    float mass=m_particleMass+_massRandomisation*massRand;

    Particle* newParticle=new Particle(pos, m_particleInitVelocity, mass, m_particleRadius, _particleInitialTemperature, this);
    m_particles.push_back(newParticle);
  }
}

void Emitter::setParticleBurningParameters(float _burnThreshold, float _burnRate, float _thermalConductivity, float _thermalMass, float _heatRelease, float _volumeRelease, float _sootThreshold)
{
  m_burnThreshold=_burnThreshold;
  m_burnRate=_burnRate;
  m_thermalConductivity=_thermalConductivity;
  m_thermalMass=_thermalMass;
  m_heatRelease=_heatRelease;
  m_volumeRelease=_volumeRelease;
  m_sootThreshold=_sootThreshold;
}

void Emitter::update(float _dt)
{
  for (int i=0; i<m_noParticles; i++)
  {
    if (m_particles[i]->m_active==true)
    {
      m_particles[i]->update(_dt);
    }
  }
}

void Emitter::renderParticles(ngl::Mat4 _ModelMatrix_Camera, bool _exportAlembic)
{

  m_ModelMatrix_Camera=_ModelMatrix_Camera;
  for (int i=0; i<m_noParticles; i++)
  {
    if (m_particles[i]->m_active==true)
    {
     m_particles[i]->render();
    }
  }

  //Export frame to alembic file if _exportAlembic is true
  if (_exportAlembic==true)
  {
    //Create list of particle positions
    std::vector <ngl::Vec3> particlePositions;
    for (int i=0; i<m_noParticles; i++)
    {
      particlePositions.push_back(m_particles.at(i)->getPosition());
    }

    m_alembicExporter->exportFrame(&particlePositions, m_noParticles);
  }
//  else
//  {
//    std::cout<<"Alembic export stopped\n";
//  }

}


