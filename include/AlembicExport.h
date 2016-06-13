#ifndef ALEMBICEXPORT
#define ALEMBICEXPORT

#include <Alembic/AbcGeom/All.h>
#include <Alembic/AbcCoreOgawa/All.h>

#include <ngl/Vec3.h>

class AlembicExport
{
public:
  //---------------------------------------------------------------------------------
  /// @brief Alembic export ctor
  //---------------------------------------------------------------------------------
  AlembicExport();
  //---------------------------------------------------------------------------------
  /// @brief Alembic export dtor
  //---------------------------------------------------------------------------------
  ~AlembicExport(){;}
  //---------------------------------------------------------------------------------
  /// @brief Write frame to alembic
  //---------------------------------------------------------------------------------
  void exportFrame(std::vector<ngl::Vec3> *_particlePositions, std::vector<float> *_particleTemperatures, int _noParticles);

private:
  //---------------------------------------------------------------------------------
  /// @brief Archive for the alembic data
  //---------------------------------------------------------------------------------
  std::unique_ptr <Alembic::AbcGeom::OArchive> m_archive;
  //---------------------------------------------------------------------------------
  /// @brief Point samples given by simulation
  //---------------------------------------------------------------------------------
  std::unique_ptr <Alembic::AbcGeom::OPoints> m_pointSamples;
  //---------------------------------------------------------------------------------
  /// @brief Array of temperatures for each point
  //---------------------------------------------------------------------------------
  std::unique_ptr <Alembic::AbcGeom::OFloatArrayProperty> m_temperatures;
};

#endif // ALEMBICEXPORT

