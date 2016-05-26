#include "Emitter.h"

#include <ngl/Random.h>
#include <cmath>
#include <math.h>

Emitter::Emitter(ngl::Vec3 _position, float _radius, int _noParticles, int _emissionRate)
{
  m_position=_position;
  m_radius=_radius;

  m_noParticles=_noParticles;
  m_emissionRate=_emissionRate;


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

void Emitter::setUpParticles(float _particleMass, ngl::Vec3 _initialParticleVelocity, float _particleLifeTime, float _particleRadius, float _particleDragConstant, float _particleInitialTemperature)
{
  m_particleMass=_particleMass;
  m_particleInitVelocity=_initialParticleVelocity;
  m_particleLifeTime=_particleLifeTime;
  m_particleRadius=_particleRadius;
  m_particleDragConstant=_particleDragConstant;

  for (int i=0; i<m_noParticles; i++)
  {
    //Create random position within disk
    ngl::Random* rand=ngl::Random::instance();
    rand->setSeed(i);
    float randRadius=rand->randomPositiveNumber(m_radius);
    float randTheta=rand->randomPositiveNumber(2*M_PI);
    float randPhi=rand->randomPositiveNumber(M_PI);
    ngl::Vec3 pos=m_position;
    pos.m_x+=randRadius*cos(randTheta)*sin(randPhi);
    pos.m_y+=randRadius*cos(randPhi);
    pos.m_z+=randRadius*sin(randTheta)*sin(randPhi);

//    ngl::Random* rand=ngl::Random::instance();
//    rand->setSeed(i);
//    ngl::Vec3 pos;
//    pos.m_x=rand->randomNumber(2.5);
//    rand->setSeed(i+100);
//    pos.m_y=rand->randomNumber(2.5);
//    rand->setSeed(i+200);
//    pos.m_z=rand->randomNumber(2.5);

    ///To do: Randomise velocity slightly?
    //Initial velocity noise
    ngl::Vec3 velocity=m_particleInitVelocity;
//    float velNoise=rand->randomNumber(5.0);
//    velocity.m_x+=velNoise;
//    velocity.m_z+=velNoise;

    //Give random lifetime
    float lifeTime=m_particleLifeTime;
    lifeTime=rand->randomPositiveNumber(lifeTime);

    Particle* newParticle=new Particle(pos, velocity, m_particleMass, m_particleRadius, lifeTime, _particleInitialTemperature, this);
    m_particles.push_back(newParticle);
  }
}

void Emitter::setParticleBurningParameters(float _burnThreshold, float _burnRate, float _thermalConductivity, float _thermalMass, float _heatRelease, float _volumeRelease, float _sootCreationConstant, float _sootThreshold)
{
  m_burnThreshold=_burnThreshold;
  m_burnRate=_burnRate;
  m_thermalConductivity=_thermalConductivity;
  m_thermalMass=_thermalMass;
  m_heatRelease=_heatRelease;
  m_volumeRelease=_volumeRelease;
  m_sootCreation=_sootCreationConstant;
  m_sootThreshold=_sootThreshold;
}

void Emitter::update(float _dt)
{
//  int countActive=0;
//  for (int i=0; i<m_noParticles; i++)
//  {
//    if (m_particles[i]->m_active==true)
//    {
//      m_particles[i]->update(_dt);
//    }
//    else if(countActive<m_emissionRate)
//    {
//      m_particles[i]->m_active=true;
//      countActive+=1;
//    }

//  }
  for (int i=0; i<m_noParticles; i++)
  {
    m_particles[i]->update(_dt);
  }


}

void Emitter::renderParticles(ngl::Mat4 _ModelMatrix_Camera)
{
  m_ModelMatrix_Camera=_ModelMatrix_Camera;
//  for (int i=0; i<m_noParticles; i++)
//  {
//    if (m_particles[i]->m_active==true)
//    {
//     m_particles[i]->render();
//    }
//  }
  for (int i=0; i<m_noParticles; i++)
  {
    m_particles[i]->render();
  }

}


