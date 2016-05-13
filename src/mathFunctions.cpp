#include "Grid.h"

#include <vector>

ngl::Vec3 Grid::RK2_integrator(ngl::Vec3 _u, ngl::Vec3 _du, float _dt)
{
  ngl::Vec3 result;

  ///Option 1:
  /// Think this is actually midpoint rule/modified Euler
//  //k1=h*f(xn,yn)
//  ngl::Vec3 k1=_du;

//  //k2=h*f(xn+h/2,yn+k1/2)
//  ngl::Vec3 k2=_du+(_dt/2.0)*k1;

//  //yn+1=yn+k2
//  result=_u+(_dt*k2);


  ///Option 2:
  /// Not sure if different. Think the first one is midpoint rule/modified Euler, and the one below is actually RK2

  //k1=f(tn,un)
  ngl::Vec3 k1=_du;

  //k2=f(tn+h,un+hk1)
  ngl::Vec3 k2=_du+(_dt*k1);

  //un+1=un+h/2*(k1+k2)
  result=_u+((_dt/2)*(k1+k2));


  return result;

}

float Grid::linearInterp(std::vector<float> *_function, float _x)
{
  float result=0;

  return result;

}

ngl::Vec3 Grid::trilinearInterpVec3(std::string _fieldType, int _index0_X, int _index0_Y, int _index0_Z, int _index1_X, int _index1_Y, int _index1_Z, ngl::Vec3 _indexActual)
{
  ngl::Vec3 result=ngl::Vec3(0.0,0.0,0.0);

  //Calculate differences to be used in the interpolation
  float x_diff=(_indexActual.m_x-_index0_X)/(_index1_X-_index0_X);
  float y_diff=(_indexActual.m_y-_index0_Y)/(_index1_Y-_index0_Y);
  float z_diff=(_indexActual.m_z-_index0_Z)/(_index1_Z-_index0_Z);

  //Find indices to be used when finding values of the function
  int index_000=getVectorIndex(_index0_X, _index0_Y, _index0_Z);
  int index_100=getVectorIndex(_index1_X, _index0_Y, _index0_Z);
  int index_001=getVectorIndex(_index0_X, _index0_Y, _index1_Z);
  int index_101=getVectorIndex(_index1_X, _index0_Y, _index1_Z);
  int index_110=getVectorIndex(_index1_X, _index1_Y, _index0_Z);
  int index_010=getVectorIndex(_index0_X, _index1_Y, _index0_Z);
  int index_011=getVectorIndex(_index0_X, _index1_Y, _index1_Z);
  int index_111=getVectorIndex(_index1_X, _index1_Y, _index1_Z);

  //Find function values for the various indices
//  ngl::Vec3 function_000=_function->at(index_000);
//  ngl::Vec3 function_100=_function->at(index_100);
//  ngl::Vec3 function_001=_function->at(index_001);
//  ngl::Vec3 function_101=_function->at(index_101);
//  ngl::Vec3 function_110=_function->at(index_110);
//  ngl::Vec3 function_010=_function->at(index_010);
//  ngl::Vec3 function_011=_function->at(index_011);
//  ngl::Vec3 function_111=_function->at(index_111);

  //Find function values for the various indices
  ngl::Vec3 function_000;
  ngl::Vec3 function_100;
  ngl::Vec3 function_001;
  ngl::Vec3 function_101;
  ngl::Vec3 function_110;
  ngl::Vec3 function_010;
  ngl::Vec3 function_011;
  ngl::Vec3 function_111;

  if (_fieldType=="velocity")
  {
    function_000=m_listCells.at(index_000)->m_newVelocity;
    function_100=m_listCells.at(index_100)->m_newVelocity;
    function_001=m_listCells.at(index_001)->m_newVelocity;
    function_101=m_listCells.at(index_101)->m_newVelocity;
    function_110=m_listCells.at(index_110)->m_newVelocity;
    function_010=m_listCells.at(index_010)->m_newVelocity;
    function_011=m_listCells.at(index_011)->m_newVelocity;
    function_111=m_listCells.at(index_111)->m_newVelocity;
  }

  //Interpolate along x
  ngl::Vec3 interpX_00=((1-x_diff)*function_000)+(x_diff*function_100);
  ngl::Vec3 interpX_01=((1-x_diff)*function_001)+(x_diff*function_101);
  ngl::Vec3 interpX_10=((1-x_diff)*function_010)+(x_diff*function_110);
  ngl::Vec3 interpX_11=((1-x_diff)*function_011)+(x_diff*function_111);

  //Interpolate along y
  ngl::Vec3 interpY_0=((1-y_diff)*interpX_00)+(y_diff*interpX_10);
  ngl::Vec3 interpY_1=((1-y_diff)*interpX_01)+(y_diff*interpX_11);

  //Interpolate along z
  ngl::Vec3 interpZ=((1-z_diff)*interpY_0)+(z_diff*interpY_1);

  //Return result
  result=interpZ;
  return result;
}

float Grid::trilinearInterpFloat(std::vector<float> *_function, int _index0_X, int _index0_Y, int _index0_Z, int _index1_X, int _index1_Y, int _index1_Z, ngl::Vec3 _indexActual)
{
  float result=0;

  //Calculate differences to be used in the interpolation
  float x_diff=(_indexActual.m_x-_index0_X)/(_index1_X-_index0_X);
  float y_diff=(_indexActual.m_y-_index0_Y)/(_index1_Y-_index0_Y);
  float z_diff=(_indexActual.m_z-_index0_Z)/(_index1_Z-_index0_Z);

  //Find indices to be used when finding values of the function
  int index_000=getVectorIndex(_index0_X, _index0_Y, _index0_Z);
  int index_100=getVectorIndex(_index1_X, _index0_Y, _index0_Z);
  int index_001=getVectorIndex(_index0_X, _index0_Y, _index1_Z);
  int index_101=getVectorIndex(_index1_X, _index0_Y, _index1_Z);
  int index_110=getVectorIndex(_index1_X, _index1_Y, _index0_Z);
  int index_010=getVectorIndex(_index0_X, _index1_Y, _index0_Z);
  int index_011=getVectorIndex(_index0_X, _index1_Y, _index1_Z);
  int index_111=getVectorIndex(_index1_X, _index1_Y, _index1_Z);

  //Find function values for the various indices
  float function_000=_function->at(index_000);
  float function_100=_function->at(index_100);
  float function_001=_function->at(index_001);
  float function_101=_function->at(index_101);
  float function_110=_function->at(index_110);
  float function_010=_function->at(index_010);
  float function_011=_function->at(index_011);
  float function_111=_function->at(index_111);

  //Interpolate along x
  float interpX_00=((1-x_diff)*function_000)+(x_diff*function_100);
  float interpX_01=((1-x_diff)*function_001)+(x_diff*function_101);
  float interpX_10=((1-x_diff)*function_010)+(x_diff*function_110);
  float interpX_11=((1-x_diff)*function_011)+(x_diff*function_111);

  //Interpolate along y
  float interpY_0=((1-y_diff)*interpX_00)+(y_diff*interpX_10);
  float interpY_1=((1-y_diff)*interpX_01)+(y_diff*interpX_11);

  //Interpolate along z
  float interpZ=((1-z_diff)*interpY_0)+(z_diff*interpY_1);

  //Return result
  result=interpZ;
  return result;
}


void Grid::linearSystemSolve(std::vector<float> *result, std::vector<float> *_A, std::vector<float> *_b)
{

}

ngl::Vec3 Grid::calcDivergenceVec3(std::string _fieldType, int _index_X, int _index_Y, int _index_Z)
{
  ngl::Vec3 result;

  if (_fieldType=="velocity")
  {
    ngl::Vec3 difference_X=m_listCells.at(getVectorIndex((_index_X+1), _index_Y, _index_Z))-m_listCells.at(getVectorIndex((_index_X-1), _index_Y, _index_Z));
    ngl::Vec3 difference_Y=m_listCells.at(getVectorIndex(_index_X, (_index_Y+1), _index_Z))-m_listCells.at(getVectorIndex(_index_X, (_index_Y-1), _index_Z));
    ngl::Vec3 difference_Z=m_listCells.at(getVectorIndex(_index_X+1, _index_Y, (_index_Z+1)))-m_listCells.at(getVectorIndex(_index_X, _index_Y, (_index_Z-1)));

    float constant=0.5/m_cellSize;

    result=constant*(difference_X + difference_Y + difference_Z);
  }

  return result;
}
