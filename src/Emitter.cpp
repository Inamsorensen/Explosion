#include "Emitter.h"

#include <ngl/Random.h>
#include <cmath>
#include <math.h>

namespace AbcG=Alembic::AbcGeom;

Emitter::Emitter(ngl::Vec3 _position, float _radius, int _noParticles)
{
  m_position=_position;
  m_radius=_radius;

  m_noParticles=_noParticles;

//  m_alembicExporter=new AlembicExport();
  m_alembicExporter.reset(new AlembicExport());

////  m_alembicArchive=AbcG::OArchive(Alembic::AbcCoreOgawa::WriteArchive(),"Explosion.abc");
//  m_alembicArchive.reset(new AbcG::OArchive(Alembic::AbcCoreOgawa::WriteArchive(),"Explosion.abc"));

//  AbcG::OObject archiveTop=m_alembicArchive->getTop();
////  AbcG::OObject archiveTop(m_alembicArchive.get(),AbcG::kTop);

//  const Alembic::Abc::chrono_t dt=1.0/25.0;
//  AbcG::TimeSampling timeSample(dt, 0.0);

//  Alembic::Util::uint32_t tsidx=archiveTop.getArchive().addTimeSampling(timeSample);

////  m_alembicPoints=AbcG::OPoints(archiveTop, "particles", tsidx);
//  m_alembicPoints.reset(new AbcG::OPoints(archiveTop, "particles", tsidx));



}

Emitter::~Emitter()
{
  int noParticlesCurr=m_particles.size();

  if(noParticlesCurr!=0)
  {
    std::cout<<"Deleting particles\n";
    for (int i=0; i<noParticlesCurr; i++)
    {
      delete m_particles[i];
    }
  }
  m_particles.clear();

  if (m_alembicExporter!=nullptr)
  {
    std::cout<<"Removing alembic export pointer\n";
//    delete m_alembicExporter;
  }

//  if (m_alembicArchive!=nullptr)
//  {
//    std::cout<<"Removing alembic archive pointer\n";
//    delete m_alembicArchive;
//  }
//  if (m_alembicPoints!=nullptr)
//  {
//    std::cout<<"Removing alembic point object pointer\n";
//    delete m_alembicPoints;
//  }

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
    std::vector <float> particleTemperatures;
    for (int i=0; i<m_noParticles; i++)
    {
      particlePositions.push_back(m_particles.at(i)->getPosition());
      particleTemperatures.push_back(m_particles.at(i)->getTemperature());
    }

    m_alembicExporter->exportFrame(&particlePositions, &particleTemperatures, m_noParticles);

//    std::vector<Imath::V3f> positions;
//    std::vector<Alembic::Util::uint64_t> ids;

//    unsigned int unsignedNoParticles=m_noParticles;

//    for (unsigned int i=0; i<unsignedNoParticles; i++)
//    {
//      ngl::Vec3 particlePos=m_particles[i]->getPosition();
//      positions.push_back(Imath::V3f(particlePos.m_x, particlePos.m_y, particlePos.m_z));
//      ids.push_back(i);
//    }

//    AbcG::V3fArraySample pos(positions);
//    AbcG::UInt64ArraySample id(ids);
//    AbcG::OPointsSchema::Sample particleSample(pos, id);

//    AbcG::OPointsSchema particleSchema=m_alembicPoints->getSchema();
//    particleSchema.set(particleSample);

  }
//  else
//  {
//    std::cout<<"Alembic export stopped, \n";
//  }

}


