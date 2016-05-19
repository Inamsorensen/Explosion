#ifndef MATHFUNCTION
#define MATHFUNCTION

#include <vector>

#include <ngl/Vec3.h>

#include "Grid.h"

class mathFunction
{
public:
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get vector index from coordinates (i,j,k)
  //----------------------------------------------------------------------------------------------------------------------
  static int getVectorIndex(int i, int j, int k, int _noCells);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Runge Kutta 2nd order integration for vec3
  //----------------------------------------------------------------------------------------------------------------------
  static ngl::Vec3 RK2_integrator(ngl::Vec3 _u, ngl::Vec3 _du, float _dt);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Linear interpolation
  //----------------------------------------------------------------------------------------------------------------------
  static float linearInterp(std::vector<float> *_function, float _x);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Trilinear interpolation for function
  //----------------------------------------------------------------------------------------------------------------------
  static ngl::Vec3 trilinearInterpVec3(std::vector<ngl::Vec3> *_function, int _index0_X, int _index0_Y, int _index0_Z, int _index1_X, int _index1_Y, int _index1_Z, ngl::Vec3 _indexActual, int _noCells);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Trilinear interpolation for function
  //----------------------------------------------------------------------------------------------------------------------
  static float trilinearInterpFloat(std::vector<float> *_function, int _index0_X, int _index0_Y, int _index0_Z, int _index1_X, int _index1_Y, int _index1_Z, ngl::Vec3 _indexActual, int _noCells);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate divergence of a Vec3 field at index (i,j,k)
  //----------------------------------------------------------------------------------------------------------------------
  static float calcDivergenceVec3(std::vector<ngl::Vec3> *_field, int _index_X, int _index_Y, int _index_Z, float _cellSize, int _noCells);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate divergence of a float field at index (i,j,k)
  //----------------------------------------------------------------------------------------------------------------------
  static ngl::Vec3 calcDivergenceFloat(std::vector<float> *_field, int _index_X, int _index_Y, int _index_Z, float _cellSize, int _noCells);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Calculate curl of a vector. In this case, used to get curl of velocity.
  //----------------------------------------------------------------------------------------------------------------------
  static void calculateCurl(std::vector<ngl::Vec3> *_field, std::vector<ngl::Vec3> *curlField, std::vector<float> *curlMagnitude, int _noCells, float _cellSize);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Linear system solver to calculate resulting float field
  //----------------------------------------------------------------------------------------------------------------------
  static void linearSystemSolveFloat(std::vector<float> *result, std::vector<float> *_initField, std::vector<float> *_b, float _Aii, float _Aij, int _iterations, int _noCells, float _setMinimumValue);
};

#endif // MATHFUNCTION

