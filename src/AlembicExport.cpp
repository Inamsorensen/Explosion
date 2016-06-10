#include "AlembicExport.h"

namespace AbcG=Alembic::AbcGeom;


AlembicExport::AlembicExport()
{
  //Create alembic geometry output archive called explosion.abc
  m_archive.reset(new AbcG::OArchive(Alembic::AbcCoreOgawa::WriteArchive(),"explosion.abc"));

  //Set up timesampling to be for 25fps
  AbcG::TimeSampling timeSamp(1.0f/25.0f,0.0f);

  //Get top of archive
  AbcG::OObject topObj(*m_archive.get(), AbcG::kTop);

  //Add in time sampling
  Alembic::Util::uint32_t timeSampAdd=topObj.getArchive().addTimeSampling(timeSamp);

  //Reset particle position output
  m_pointSamples.reset(new AbcG::OPoints(topObj,"simpleParticles", timeSampAdd));

  ///To do: Add temperature property

}

void AlembicExport::exportFrame(std::vector<ngl::Vec3> *_particlePositions, int _noParticles)
{
  ///To do: Set temperature properties as well

  //Set up Imath vectors for position and id data
  std::vector<Imath::V3f> positions;
  std::vector<Alembic::Util::uint64_t> id;

  unsigned int noParticles=_noParticles;

  for (unsigned int i=0; i<noParticles; i++)
  {
    float posX=_particlePositions->at(i).m_x;
    float posY=_particlePositions->at(i).m_y;
    float posZ=_particlePositions->at(i).m_z;
    positions.push_back(Imath::V3f(posX, posY, posZ));
    id.push_back(i);
  }

  //Need to create samples from data
  AbcG::V3fArraySample samplePos(positions);
  AbcG::UInt64ArraySample sampleId(id);
  AbcG::OPointsSchema::Sample pointsSample(samplePos, sampleId);

  m_pointSamples->getSchema().set(pointsSample);

}
