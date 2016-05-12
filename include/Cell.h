#ifndef CELL
#define CELL

#include <ngl/Vec3.h>

struct Cell
{
  //---------------------------------------------------------------------------------
  /// @brief True if cell is solid boundary, false if fluid
  //---------------------------------------------------------------------------------
  bool m_isSolid;
  //---------------------------------------------------------------------------------
  /// @brief Old velocity field vector
  //---------------------------------------------------------------------------------
  ngl::Vec3 m_oldVelocity;
  //---------------------------------------------------------------------------------
  /// @brief Updated velocity field vector
  //---------------------------------------------------------------------------------
  ngl::Vec3 m_newVelocity;
  //---------------------------------------------------------------------------------
  /// @brief External force field vector
  //---------------------------------------------------------------------------------
  ngl::Vec3 m_externalForce;
  //---------------------------------------------------------------------------------
  /// @brief Cell pressure
  //---------------------------------------------------------------------------------
  float m_pressure;

  //---------------------------------------------------------------------------------
  /// @brief Stores temporary vector
  //---------------------------------------------------------------------------------
  ngl::Vec3 m_tempStoreVec3;
  //---------------------------------------------------------------------------------
  /// @brief Store temporary float
  //---------------------------------------------------------------------------------
  float m_tempStoreFloat;

};

#endif //CELL
