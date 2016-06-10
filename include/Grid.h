#ifndef GRID
#define GRID

#include <vector>
#include <ngl/Vec3.h>
#include <iostream>

#include "mathFunction.h"

/// @brief Grid.h
/// Grid providing a Eulerian Grid calculation of Navier-Stokes
/// Contains fields for velocity, pressure, forces and temperature
/// Updates values through advection, pressure projection and diffusion
/// Author: Ina M. Sorensen
/// Date: 27.05.16
/// Implementation based on the code by
/// Christian Miller (2007). Realtime Explosion Simulation [online].
/// [Accessed May 2016]. Available from: <http://www.cs.utexas.edu/~ckm/explosion/>.

class Grid
{
public:
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Grid create instance
  //----------------------------------------------------------------------------------------------------------------------
  static Grid* createGrid(ngl::Vec3 _origin, float _gridSize, int _noCells);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Grid instance
  //----------------------------------------------------------------------------------------------------------------------
  static Grid* getGrid();

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get grid position
  //----------------------------------------------------------------------------------------------------------------------
  inline ngl::Vec3 getPosition() const {return m_origin;}
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get grid size
  //----------------------------------------------------------------------------------------------------------------------
  inline float getGridSize() const {return m_gridSize;}
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get cell size
  //----------------------------------------------------------------------------------------------------------------------
  inline float getCellSize() const {return m_cellSize;}

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Setup explosion by giving position of half sphere, its radius and the temperature and divergence addition
  //----------------------------------------------------------------------------------------------------------------------
  void setExplosion(ngl::Vec3 _originExplosion, float _radius, float _temperature, float _divAddition);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set temperature field constants
  //----------------------------------------------------------------------------------------------------------------------
  void setTemperatureFieldConstants(float _ambientTemperature, float _thermalConductivity, float _coolingConstant, float _buoyancyConstant);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set velocity field constants
  //----------------------------------------------------------------------------------------------------------------------
  void setVelocityFieldConstants(float _noise, float _vorticityConstant);


  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Update grid after time step _dt
  //----------------------------------------------------------------------------------------------------------------------
  void update(float _dt);

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get vector index from coordinates (i,j,k)
  //----------------------------------------------------------------------------------------------------------------------
  int getVectorIndex(int i, int j, int k);


  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get velocity from field based on input position
  //----------------------------------------------------------------------------------------------------------------------
  ngl::Vec3 getVelocityFromField(ngl::Vec3 _particlePosition);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Get temperature from field based on input position
  //----------------------------------------------------------------------------------------------------------------------
  float getTemperatureFromField(ngl::Vec3 _particlePosition);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Add temperature to field from particle
  //----------------------------------------------------------------------------------------------------------------------
  void addTemperatureFromParticle(int _index, float _temperature);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Add pressure from combustion gas created by fuel particle
  //----------------------------------------------------------------------------------------------------------------------
  void addDivergenceFromParticle(int _index, float _divergence);


private:
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Grid constructor
  //----------------------------------------------------------------------------------------------------------------------
  Grid(ngl::Vec3 _origin, float _gridSize, int _noCells);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Grid instance pointer
  //----------------------------------------------------------------------------------------------------------------------
  static Grid* m_instance;


  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Origin (set as back lower corner)
  //----------------------------------------------------------------------------------------------------------------------
  ngl::Vec3 m_origin;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Length of one side of the grid (same for all directions as cube)
  //----------------------------------------------------------------------------------------------------------------------
  float m_gridSize;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Number of cells along one side (same in all direction, total number of cells=noCells^3)
  //----------------------------------------------------------------------------------------------------------------------
  int m_noCells;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Length of one grid cell (again cubic so same in all directions)
  //----------------------------------------------------------------------------------------------------------------------
  float m_cellSize;


  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Vorticity constant
  //----------------------------------------------------------------------------------------------------------------------
  float m_vorticityConstant;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Noise constant
  //----------------------------------------------------------------------------------------------------------------------
  float m_noiseConstant;

  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Ambient temperature. Set to 20C=293K
  //----------------------------------------------------------------------------------------------------------------------
  float m_ambientTemp;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Maximum temperature. Set through explosion setup
  //----------------------------------------------------------------------------------------------------------------------
  float m_maxTemp;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Thermal conductivity. Should give how much diffusion of temperature from/through the fluid
  //----------------------------------------------------------------------------------------------------------------------
  float m_thermalConductivity;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Cooling constant. Heat given off to the environment
  //----------------------------------------------------------------------------------------------------------------------
  float m_coolingConstant;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Buoyancy constant. The amount of effect of temperature on fluid velocity
  //----------------------------------------------------------------------------------------------------------------------
  float m_buoyancyConstant;


  //----------------------------------------------------------------------------------------------------------------------
  /// @brief New velocity field
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<ngl::Vec3> m_newVelocityField;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Force field
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<ngl::Vec3> m_forceField;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Pressure field
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<float> m_pressureField;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief New temperature field
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<float> m_newTemperatureField;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Temporarily store Vec3 field
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<ngl::Vec3> m_storeFieldVec3;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Temporarily store float field
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<float> m_storeFieldFloat;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Store zero vector field
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<ngl::Vec3> m_storeZeroFieldVec3;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Store zero float field
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<float> m_storeZeroFieldFloat;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Store b for calculating projection
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<float> m_storeB;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Store curl vectors and curl magnitudes
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<ngl::Vec3> m_curlVectorField;
  std::vector<float> m_curlMagnitudes;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Store divergence field
  ///        Calculates interaction of explosion with fluid/air by adding this to b in projection
  //----------------------------------------------------------------------------------------------------------------------
  std::vector<float> m_divergenceField;


  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set boundary values for _vec3Field equal to the nearest neighbour cell
  //----------------------------------------------------------------------------------------------------------------------
  void setBoundaryValuesVec3(std::vector<ngl::Vec3> *_vec3Field);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set boundary values for _floatField equal to the nearest neighbour cell
  //----------------------------------------------------------------------------------------------------------------------
  void setBoundaryValuesFloat(std::vector<float> *_floatField);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set boundary velocity values to zero in direction of boundary normal
  //----------------------------------------------------------------------------------------------------------------------
  void setBoundaryVelocity();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Set boundary temperature equal to ambient temperature. Needs to be changed if should be solid temperature
  //----------------------------------------------------------------------------------------------------------------------
  void setBoundaryTemperature();


  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Update velocity field
  //----------------------------------------------------------------------------------------------------------------------
  void updateVelocityField(float _dt);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Step1: Add force
  //----------------------------------------------------------------------------------------------------------------------
  void addForce(float _dt);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Step2: Advect velocity
  //----------------------------------------------------------------------------------------------------------------------
  void advectVelocity(float _dt);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Step4: Project
  //----------------------------------------------------------------------------------------------------------------------
  void project(float _dt);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Function to calculate vorticity forces
  //----------------------------------------------------------------------------------------------------------------------
  void calculateVorticity();
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Function to calculate buoyancy forces from temperature. Done before temperature update so will use
  /// temperature field from previous step
  //----------------------------------------------------------------------------------------------------------------------
  void calculateBuoyancy();


  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Update temperature field
  //----------------------------------------------------------------------------------------------------------------------
  void updateTemperatureField(float _dt);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Step1: Advect temperature
  //----------------------------------------------------------------------------------------------------------------------
  void advectTemperature(float _dt);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Step2: Diffuse temperature
  //----------------------------------------------------------------------------------------------------------------------
  void diffuseTemperature(float _dt);
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Step3: Cool down through dissipation to environment
  //----------------------------------------------------------------------------------------------------------------------
  void dissipateTemperature();
  
};

#endif //GRID
