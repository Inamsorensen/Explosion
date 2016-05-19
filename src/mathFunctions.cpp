#include "mathFunction.h"

#include <vector>
#include <cmath>

int mathFunction::getVectorIndex(int i, int j, int k, int _noCells)
{
  int vectorIndex=i+(_noCells*j)+(_noCells*_noCells*k);
  return vectorIndex;
}

ngl::Vec3 mathFunction::RK2_integrator(ngl::Vec3 _u, ngl::Vec3 _du, float _dt)
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

float mathFunction::linearInterp(std::vector<float> *_function, float _x)
{
  float result=0;

  return result;

}

ngl::Vec3 mathFunction::trilinearInterpVec3(std::vector<ngl::Vec3> *_function, int _index0_X, int _index0_Y, int _index0_Z, int _index1_X, int _index1_Y, int _index1_Z, ngl::Vec3 _indexActual, int _noCells)
{
  ngl::Vec3 result=ngl::Vec3(0.0,0.0,0.0);

  //Calculate differences to be used in the interpolation
  float x_diff=(_indexActual.m_x-_index0_X)/(_index1_X-_index0_X);
  float y_diff=(_indexActual.m_y-_index0_Y)/(_index1_Y-_index0_Y);
  float z_diff=(_indexActual.m_z-_index0_Z)/(_index1_Z-_index0_Z);

  //Need to set difference values to positive, so 1-diff gives a ratio, ie. 1-diff<1, not 1+diff.
  x_diff=std::abs(x_diff);
  y_diff=std::abs(y_diff);
  z_diff=std::abs(z_diff);

  //Find indices to be used when finding values of the function
  int index_000=getVectorIndex(_index0_X, _index0_Y, _index0_Z, _noCells);
  int index_100=getVectorIndex(_index1_X, _index0_Y, _index0_Z, _noCells);
  int index_001=getVectorIndex(_index0_X, _index0_Y, _index1_Z, _noCells);
  int index_101=getVectorIndex(_index1_X, _index0_Y, _index1_Z, _noCells);
  int index_110=getVectorIndex(_index1_X, _index1_Y, _index0_Z, _noCells);
  int index_010=getVectorIndex(_index0_X, _index1_Y, _index0_Z, _noCells);
  int index_011=getVectorIndex(_index0_X, _index1_Y, _index1_Z, _noCells);
  int index_111=getVectorIndex(_index1_X, _index1_Y, _index1_Z, _noCells);

  //Find function values for the various indices
  ngl::Vec3 function_000=_function->at(index_000);
  ngl::Vec3 function_100=_function->at(index_100);
  ngl::Vec3 function_001=_function->at(index_001);
  ngl::Vec3 function_101=_function->at(index_101);
  ngl::Vec3 function_110=_function->at(index_110);
  ngl::Vec3 function_010=_function->at(index_010);
  ngl::Vec3 function_011=_function->at(index_011);
  ngl::Vec3 function_111=_function->at(index_111);

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

//  std::cout<<

  //Return result
  result=interpZ;
  return result;
}

float mathFunction::trilinearInterpFloat(std::vector<float> *_function, int _index0_X, int _index0_Y, int _index0_Z, int _index1_X, int _index1_Y, int _index1_Z, ngl::Vec3 _indexActual, int _noCells)
{
  float result=0;

  //Calculate differences to be used in the interpolation
  float x_diff=(_indexActual.m_x-_index0_X)/(_index1_X-_index0_X);
  float y_diff=(_indexActual.m_y-_index0_Y)/(_index1_Y-_index0_Y);
  float z_diff=(_indexActual.m_z-_index0_Z)/(_index1_Z-_index0_Z);

  //Need to set difference values to positive, so 1-diff gives a ratio, ie. 1-diff<1, not 1+diff.
  x_diff=std::abs(x_diff);
  y_diff=std::abs(y_diff);
  z_diff=std::abs(z_diff);

  //Find indices to be used when finding values of the function
  int index_000=getVectorIndex(_index0_X, _index0_Y, _index0_Z, _noCells);
  int index_100=getVectorIndex(_index1_X, _index0_Y, _index0_Z, _noCells);
  int index_001=getVectorIndex(_index0_X, _index0_Y, _index1_Z, _noCells);
  int index_101=getVectorIndex(_index1_X, _index0_Y, _index1_Z, _noCells);
  int index_110=getVectorIndex(_index1_X, _index1_Y, _index0_Z, _noCells);
  int index_010=getVectorIndex(_index0_X, _index1_Y, _index0_Z, _noCells);
  int index_011=getVectorIndex(_index0_X, _index1_Y, _index1_Z, _noCells);
  int index_111=getVectorIndex(_index1_X, _index1_Y, _index1_Z, _noCells);

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


void mathFunction::linearSystemSolveFloat(std::vector<float> *result, std::vector<float> *_initField, std::vector<float> *_b, float _Aii, float _Aij, int _iterations, int _noCells, float _setMinimumValue)
{
  //Set up variables
  float p_ijk;
  float p_i1jk;
  float p_i0jk;
  float p_ij1k;
  float p_ij0k;
  float p_ijk1;
  float p_ijk0;

  int index;

  float sumNeighbours;
  float newValue;


  //Loop over iterations
  for (int iterationCount=0; iterationCount<_iterations; iterationCount++)
  {
    //Loop over i,j,k
    for (int k=1; k<(_noCells-1); k++)
    {
      for (int j=1; j<(_noCells-1); j++)
      {
        for (int i=1; i<(_noCells-1); i++)
        {
          //Get index for this cell (i,j,k)
          index=getVectorIndex(i,j,k, _noCells);

          //Find initField of itself and 6 nearest neighbours
          p_ijk=_initField->at(index);
          p_i1jk=_initField->at(getVectorIndex(i+1,j,k, _noCells));
          p_i0jk=_initField->at(getVectorIndex(i-1,j,k, _noCells));
          p_ij1k=_initField->at(getVectorIndex(i,j+1,k, _noCells));
          p_ij0k=_initField->at(getVectorIndex(i,j-1,k, _noCells));
          p_ijk1=_initField->at(getVectorIndex(i,j,k+1, _noCells));
          p_ijk0=_initField->at(getVectorIndex(i,j,k-1, _noCells));

          //Add nearest neighbour pressures
          sumNeighbours=_Aij*(p_i1jk + p_i0jk + p_ij1k + p_ij0k + p_ijk1 + p_ijk0);

          //Subtract from b to find new pressure at (i,j,k)
          newValue=(1.0/_Aii)*(_b->at(index)-sumNeighbours);
//          newValue=(1.0/_Aii)*(sumNeighbours+p_ijk);

          //Set result at (i,j,k) to this new pressure value
          ///Tried to set to minimum value if below this. Not sure if correct
          if (newValue>=_setMinimumValue)
          {
            result->at(index)=newValue;
          }
          else
          {
            result->at(index)=_setMinimumValue;
          }
          ///This is probably what should be used
//          result->at(index)=newValue;

//          result->at(index)=std::abs(newValue);

        }

      }
    }

    //Check for convergence?

    //Set values found to start for next iteration
    _initField=result;
  }

}

float mathFunction::calcDivergenceVec3(std::vector<ngl::Vec3> *_field, int _index_X, int _index_Y, int _index_Z, float _cellSize, int _noCells)
{
  float result;

  ngl::Vec3 fieldValue_i1jk=_field->at(getVectorIndex((_index_X+1), _index_Y, _index_Z, _noCells));
  ngl::Vec3 fieldValue_i0jk=_field->at(getVectorIndex((_index_X-1), _index_Y, _index_Z, _noCells));
  ngl::Vec3 fieldValue_ij1k=_field->at(getVectorIndex(_index_X, (_index_Y+1), _index_Z, _noCells));
  ngl::Vec3 fieldValue_ij0k=_field->at(getVectorIndex(_index_X, (_index_Y-1), _index_Z, _noCells));
  ngl::Vec3 fieldValue_ijk1=_field->at(getVectorIndex(_index_X, _index_Y, (_index_Z+1), _noCells));
  ngl::Vec3 fieldValue_ijk0=_field->at(getVectorIndex(_index_X, _index_Y, (_index_Z-1), _noCells));

  float dfxdx=fieldValue_i1jk.m_x-fieldValue_i0jk.m_x;
  float dfydy=fieldValue_ij1k.m_y-fieldValue_ij0k.m_y;
  float dfzdz=fieldValue_ijk1.m_z-fieldValue_ijk0.m_z;

  float constant=0.5/_cellSize;

  result=constant*(dfxdx+dfydy+dfzdz);

  return result;
}

ngl::Vec3 mathFunction::calcDivergenceFloat(std::vector<float> *_field, int _index_X, int _index_Y, int _index_Z, float _cellSize, int _noCells)
{
  ngl::Vec3 result;

  float fieldValue_i1jk=_field->at(getVectorIndex((_index_X+1), _index_Y, _index_Z, _noCells));
  float fieldValue_i0jk=_field->at(getVectorIndex((_index_X-1), _index_Y, _index_Z, _noCells));
  float fieldValue_ij1k=_field->at(getVectorIndex(_index_X, (_index_Y+1), _index_Z, _noCells));
  float fieldValue_ij0k=_field->at(getVectorIndex(_index_X, (_index_Y-1), _index_Z, _noCells));
  float fieldValue_ijk1=_field->at(getVectorIndex(_index_X, _index_Y, (_index_Z+1), _noCells));
  float fieldValue_ijk0=_field->at(getVectorIndex(_index_X, _index_Y, (_index_Z-1), _noCells));

  float dfdx=fieldValue_i1jk-fieldValue_i0jk;
  float dfdy=fieldValue_ij1k-fieldValue_ij0k;
  float dfdz=fieldValue_ijk1-fieldValue_ijk0;

  float constant=0.5/_cellSize;

  result=ngl::Vec3(dfdx, dfdy, dfdz);
  result=constant*result;

  return result;

}

void mathFunction::calculateCurl(std::vector<ngl::Vec3> *_field, std::vector<ngl::Vec3> *curlField, std::vector<float> *curlMagnitude, int _noCells, float _cellSize)
{
  //Set up needed variables
  ngl::Vec3 field_i1jk;
  ngl::Vec3 field_i0jk;
  ngl::Vec3 field_ij1k;
  ngl::Vec3 field_ij0k;
  ngl::Vec3 field_ijk1;
  ngl::Vec3 field_ijk0;

  float dfxdy;
  float dfxdz;
  float dfydx;
  float dfydz;
  float dfzdx;
  float dfzdy;

  ngl::Vec3 resultCurl;

  for (int k=1; k<(_noCells-1); k++)
  {
    for (int j=1; j<(_noCells-1); j++)
    {
      for (int i=1; i<(_noCells-1); i++)
      {
        //Find values of field at (i+1,j,k), (i-1,j,k) etc.
        field_i1jk=_field->at(getVectorIndex(i+1, j, k, _noCells));
        field_i0jk=_field->at(getVectorIndex(i-1, j, k, _noCells));
        field_ij1k=_field->at(getVectorIndex(i, j+1, k, _noCells));
        field_ij0k=_field->at(getVectorIndex(i, j-1, k, _noCells));
        field_ijk1=_field->at(getVectorIndex(i, j, k+1, _noCells));
        field_ijk0=_field->at(getVectorIndex(i, j, k-1, _noCells));

        //Calculate approximate value for df_x/dy etc
        //df_x/dy=(1/2*cellSize)*(field(i+1,j,k)-field(i-1,j,k))
        dfxdy=(0.5/_cellSize)*(field_ij1k.m_x-field_ij0k.m_x);
        dfxdz=(0.5/_cellSize)*(field_ijk1.m_x-field_ijk0.m_x);
        dfydx=(0.5/_cellSize)*(field_i1jk.m_y-field_i0jk.m_y);
        dfydz=(0.5/_cellSize)*(field_ijk1.m_y-field_ijk0.m_y);
        dfzdx=(0.5/_cellSize)*(field_i1jk.m_z-field_i0jk.m_z);
        dfzdy=(0.5/_cellSize)*(field_ij1k.m_z-field_ij0k.m_z);

        //Curl is nabla crossed with velocity
        //Curl_x=df_z/dy - df_y/dz
        resultCurl.m_x=(dfzdy-dfydz);
        resultCurl.m_y=(dfxdz-dfzdx);
        resultCurl.m_z=(dfydx-dfxdy);

        curlField->at(getVectorIndex(i,j,k, _noCells))=resultCurl;

        //Calculate magnitude of curl
        curlMagnitude->at(getVectorIndex(i,j,k, _noCells))=resultCurl.length();
      }
    }
  }
}
