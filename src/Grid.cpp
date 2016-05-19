#include "Grid.h"

#include <stdexcept>
#include <cmath>
#include <math.h>
#include <iostream>

#include <ngl/Random.h>

Grid* Grid::m_instance=nullptr;

Grid::Grid(ngl::Vec3 _origin, float _gridSize, int _noCells)
{
  //Set initial variables for grid and calculate cell size
  m_origin=_origin;
  m_gridSize=_gridSize;
  m_noCells=_noCells;
  m_cellSize=m_gridSize/m_noCells;

  //Set up space and initial values for fields
  for (int i=0; i<pow(m_noCells,3); i++)
  {
    //Zero fields for pressure and force
    m_pressureField.push_back(0.0);
    m_forceField.push_back(ngl::Vec3(0.0,0.0,0.0));

    m_oldVelocityField.push_back(ngl::Vec3(0.0,0.0,0.0));


    //Set ambient temperature everywhere
    m_oldTemperatureField.push_back(0.0);

    //Set space for store fields
    m_storeFieldVec3.push_back(ngl::Vec3(0.0,0.0,0.0));
    m_storeFieldFloat.push_back(0.0);
    m_storeZeroFieldVec3.push_back(ngl::Vec3(0.0,0.0,0.0));

    //Set space for curl fields
    m_curlVectorField.push_back(ngl::Vec3(0.0,0.0,0.0));
    m_curlMagnitudes.push_back(0.0);

    //Set space for divergence field
    m_divergenceField.push_back(0.0);

    //Set space for storeB
    m_storeB.push_back(0.0);

  }

  m_newVelocityField=m_oldVelocityField;
  m_newTemperatureField=m_oldTemperatureField;

}


Grid* Grid::createGrid(ngl::Vec3 _origin, float _gridSize, int _noCells)
{
  if(m_instance==nullptr)
  {
    m_instance=new Grid(_origin, _gridSize, _noCells);
  }

  return m_instance;
}


Grid* Grid::getGrid()
{
  if (m_instance==nullptr)
  {
    throw std::invalid_argument("You need to create the grid first.");
  }
  return m_instance;
}

void Grid::setVelocityFieldConstants(float _noise, float _vorticityConstant)
{
  m_noiseConstant=_noise;
  m_vorticityConstant=_vorticityConstant;

  //Get instance for random number generator
  ngl::Random* rand=ngl::Random::instance();
  float noiseValX;
  float noiseValY;
  float noiseValZ;
  ngl::Vec3 velocity;

  for (int i=0; i<pow(m_noCells,3); i++)
  {
    //Noise in velocity field
    rand->setSeed(i);
    noiseValX=rand->randomNumber(m_noiseConstant);
    rand->setSeed(i+1000);
    noiseValY=rand->randomNumber(m_noiseConstant);
    rand->setSeed(i+2000);
    noiseValZ=rand->randomNumber(m_noiseConstant);
  //    rand->setSeed(i);
  //    noiseValX=rand->randomPositiveNumber(_noise);
  //    rand->setSeed(i+1000);
  //    noiseValY=rand->randomPositiveNumber(_noise);
  //    rand->setSeed(i+2000);
  //    noiseValZ=rand->randomPositiveNumber(_noise);
    velocity.m_x=noiseValX;
    velocity.m_y=noiseValY;
    velocity.m_z=noiseValZ;

    //Test velocity field effect
  //    velocity.m_x=0.0;
  //    velocity.m_y=1.0;
  //    velocity.m_z=0.0;

    m_newVelocityField.at(i)=velocity;
  }

  setBoundaryVelocity();

}

void Grid::setTemperatureFieldConstants(float _ambientTemperature, float _thermalConductivity, float _coolingConstant, float _buoyancyConstant)
{
  m_ambientTemp=_ambientTemperature;
  m_thermalConductivity=_thermalConductivity;
  m_coolingConstant=_coolingConstant;
  m_buoyancyConstant=_buoyancyConstant;

  for (int i=0; i<pow(m_noCells,3); i++)
  {
    m_newTemperatureField.at(i)=m_ambientTemp;
  }

}

void Grid::setExplosion(ngl::Vec3 _originExplosion, float _radius, float _temperature, float _divAddition)
{
  m_maxTemp=_temperature;

  //Find cells that will contain explosive sphere
    //Find cell where the origin is
  ngl::Vec3 indexExplosionOrigin=(1/m_cellSize)*(_originExplosion-m_origin);
  int indexOrigin_X;
  int indexOrigin_Y;
  int indexOrigin_Z;
  indexOrigin_X=trunc(indexExplosionOrigin.m_x);
  indexOrigin_Y=trunc(indexExplosionOrigin.m_y);
  indexOrigin_Z=trunc(indexExplosionOrigin.m_z);

    //Check that index is inside grid, if not throw error
  if (indexOrigin_X<1 || indexOrigin_Y<1 || indexOrigin_Z<1 || indexOrigin_X>(m_noCells-1) || indexOrigin_Y>(m_noCells-1) || indexOrigin_Z>(m_noCells-1))
  {
    std::cout<<"ExplosionOrigin: ["<<indexOrigin_X<<" "<<indexOrigin_Y<<" "<<indexOrigin_Z<<"]\n";
    throw std::invalid_argument("Origin of explosion is in a boundary cell or outside the grid");
  }

    //Check radius against cellSize to see how many cells will be changed in one direction
  int noCellsInRadius=trunc(_radius/m_cellSize);

  //Check that explosion isn't too big compared to grid
  if (_radius>(((m_noCells-2)/2)*m_cellSize))
  {
    std::cout<<"No cells in radius: "<<noCellsInRadius<<"\n";
    throw std::invalid_argument("Explosion takes up the whole grid or more.");
  }

  //Loop over cells that are likely to be within the radius of explosion
  for (int k=(indexOrigin_Z-noCellsInRadius); k<=(indexOrigin_Z+noCellsInRadius); k++)
  {
    for (int j=(indexOrigin_Y-noCellsInRadius); j<=(indexOrigin_Y+noCellsInRadius); j++)
    {
      for (int i=(indexOrigin_X-noCellsInRadius); i<=(indexOrigin_X+noCellsInRadius); i++)
      {

        //Calculate position of current grid point relative to explosion origin
        float positionX=(i-indexOrigin_X)*m_cellSize;
        float positionY=(j-indexOrigin_Y)*m_cellSize;
        float positionZ=(k-indexOrigin_Z)*m_cellSize;

        //Check that index is inside radius of sphere
        if (pow(_radius,2)>=(pow(positionX,2)+pow(positionY,2)+pow(positionZ,2)))
        {

          //Check that cell isn't boundary cell
          if (i<1 || j<1 || k<1 || i>(m_noCells-1) || j>(m_noCells-1) || k>(m_noCells-1))
          {
            //If outside - do nothing
          }
          else
          {
          int indexCurrentCell=getVectorIndex(i,j,k);

            //Set temperature of cells found to _temperature
          m_newTemperatureField.at(indexCurrentCell)=m_maxTemp;

            //Add divAddition.
          m_divergenceField.at(indexCurrentCell)+=_divAddition;

          }
        }
      }
    }
  }
}

void Grid::update(float _dt)
{
  //Set old velocity field equal to new velocity field so new velocity field can be given new values. Do same for temperature field
  swap();

  updateVelocityField(_dt);

  updateTemperatureField(_dt);

}

int Grid::getVectorIndex(int i, int j, int k)
{
  int vectorIndex=i+(m_noCells*j)+(m_noCells*m_noCells*k);
  return vectorIndex;
}

void Grid::setBoundaryValuesVec3(std::vector<ngl::Vec3> *_vec3Field)
{
  for (int i=0; i<m_noCells; i++)
  {
    for (int j=0; j<m_noCells; j++)
    {
      _vec3Field->at(getVectorIndex(0,i,j))=_vec3Field->at(getVectorIndex(1,i,j));
      _vec3Field->at(getVectorIndex(m_noCells-1,i,j))=_vec3Field->at(getVectorIndex(m_noCells-2,i,j));
      _vec3Field->at(getVectorIndex(i,0,j))=_vec3Field->at(getVectorIndex(i,1,j));
      _vec3Field->at(getVectorIndex(i,m_noCells-1,j))=_vec3Field->at(getVectorIndex(i,m_noCells-2,j));
      _vec3Field->at(getVectorIndex(i,j,0))=_vec3Field->at(getVectorIndex(i,j,1));
      _vec3Field->at(getVectorIndex(i,j,m_noCells-1))=_vec3Field->at(getVectorIndex(i,j,m_noCells-2));

    }
  }
}

void Grid::setBoundaryValuesFloat(std::vector<float> *_floatField)
{
  for (int i=0; i<m_noCells; i++)
  {
    for (int j=0; j<m_noCells; j++)
    {
      _floatField->at(getVectorIndex(0,i,j))=_floatField->at(getVectorIndex(1,i,j));
      _floatField->at(getVectorIndex(m_noCells-1,i,j))=_floatField->at(getVectorIndex(m_noCells-2,i,j));
      _floatField->at(getVectorIndex(i,0,j))=_floatField->at(getVectorIndex(i,1,j));
      _floatField->at(getVectorIndex(i,m_noCells-1,j))=_floatField->at(getVectorIndex(i,m_noCells-2,j));
      _floatField->at(getVectorIndex(i,j,0))=_floatField->at(getVectorIndex(i,j,1));
      _floatField->at(getVectorIndex(i,j,m_noCells-1))=_floatField->at(getVectorIndex(i,j,m_noCells-2));
    }
  }
}


void Grid::setBoundaryVelocity()
{
  for (int i=0; i<m_noCells; i++)
  {
    for (int j=0; j<m_noCells; j++)
    {
      m_newVelocityField.at(getVectorIndex(0,i,j)).m_x=0;
      m_newVelocityField.at(getVectorIndex(m_noCells-1,i,j)).m_x=0;
      m_newVelocityField.at(getVectorIndex(i,0,j)).m_y=0;
      m_newVelocityField.at(getVectorIndex(i,m_noCells-1,j)).m_y=0;
      m_newVelocityField.at(getVectorIndex(i,j,0)).m_z=0;
      m_newVelocityField.at(getVectorIndex(i,j,m_noCells-1)).m_z=0;
    }
  }
}

void Grid::setBoundaryTemperature()
{
  for (int i=0; i<m_noCells; i++)
  {
    for (int j=0; j<m_noCells; j++)
    {
      m_newTemperatureField.at(getVectorIndex(0,i,j))=m_ambientTemp;
      m_newTemperatureField.at(getVectorIndex(m_noCells-1,i,j))=m_ambientTemp;
      m_newTemperatureField.at(getVectorIndex(i,0,j))=m_ambientTemp;
      m_newTemperatureField.at(getVectorIndex(i,m_noCells-1,j))=m_ambientTemp;
      m_newTemperatureField.at(getVectorIndex(i,j,0))=m_ambientTemp;
      m_newTemperatureField.at(getVectorIndex(i,j,m_noCells-1))=m_ambientTemp;
    }
  }
}


void Grid::swap()
{
  m_oldVelocityField=m_newVelocityField;

  m_oldTemperatureField=m_newTemperatureField;
}


void Grid::updateVelocityField(float _dt)
{
  //Add external forces. w1=w0+dt*f
  addForce(_dt);

  ///Code has another project in it. Not sure why but makes the result better?
//  project();

  //Advect. w2=w1(p(x,-dt))
  advectVelocity(_dt);

  //Diffuse (I-vdtD2)w3=w2

  //Project D2p=D.w3, w4=w3-Dp
  project();

  //Debug check for final velocity. Prints index and velocity
//  for (int k=0; k<m_noCells; k++)
//  {
//    for (int j=0; j<m_noCells; j++)
//    {
//      for (int i=0; i<m_noCells; i++)
//      {
//        std::cout<<"Cell: ["<<i<<" "<<j<<" "<<k<<"] Final velocity: ["<<m_newVelocityField.at(getVectorIndex(i,j,k)).m_x<<" "<<m_newVelocityField.at(getVectorIndex(i,j,k)).m_y<<" "<<m_newVelocityField.at(getVectorIndex(i,j,k)).m_z<<"]\n";
//      }
//    }
//  }
}



void Grid::addForce(float _dt)
{
  //Set force field to zero since all force calculations adds onto this field
  m_forceField=m_storeZeroFieldVec3;

  //Calculate vorticity forces
  calculateVorticity();

  //Calculate buoyancy forces
  calculateBuoyancy();

  //Add forces to velocity field
  for (int i=0; i<pow(m_noCells,3); i++)
  {
    m_newVelocityField.at(i)=m_oldVelocityField.at(i)+(m_forceField.at(i)*_dt);
  }

  //Make sure velocities at boundaries are correct
  setBoundaryVelocity();

}



void Grid::advectVelocity(float _dt)
{
  //Need to loop over 3D space i,j,k
  for (int k=0; k<m_noCells; k++)
  {
    for (int j=0; j<m_noCells; j++)
    {
      for (int i=0; i<m_noCells; i++)
      {
        //Find previous position for fictive particle at (i,j,k)
        //Get old_velocity at (i,j,k)
        int index=getVectorIndex(i,j,k);

        ngl::Vec3 oldVelocity=m_oldVelocityField.at(index);

        //Get actual position for (i,j,k)
        ngl::Vec3 currPosition;
        currPosition.m_x=(i+0.5)*m_cellSize+m_origin.m_x;
        currPosition.m_y=(j+0.5)*m_cellSize+m_origin.m_y;
        currPosition.m_z=(k+0.5)*m_cellSize+m_origin.m_z;

        //Use these values in RK2 to find x0
        ngl::Vec3 prevPosition=mathFunction::RK2_integrator(currPosition,oldVelocity,(-1)*_dt);

        //Determine the cell that is in
        //Is cell inside boundaries? If not find nearest boundary cell and use that value for advection.
        ngl::Vec3 advectVelocity=getVelocityFromField(prevPosition);

        //Need to store new values temporarily so don't interfer with the w1 we're collecting values from
        m_storeFieldVec3.at(index)=advectVelocity;

//        std::cout<<"ijk: ["<<i<<" "<<j<<" "<<k<<"] vel: ["<<advectVelocity.m_x<<" "<<advectVelocity.m_y<<" "<<advectVelocity.m_z<<"]\n";

      }
    }
  }

  //When all advection velocities calculated, set the temporarily stored field to the newVelocity field.
  m_newVelocityField=m_storeFieldVec3;

  setBoundaryVelocity();

}


void Grid::diffuse(float _dt)
{

}


void Grid::project()
{
  ///To do:
  /// Add check for convergence in loop


  //Set number of iterations
  int iterations=30;
//  int iterations=100;

  //Set up b. Same for each iteration so only do this once
  for (int k=1; k<(m_noCells-1); k++)
  {
    for (int j=1; j<(m_noCells-1); j++)
    {
      for (int i=1; i<(m_noCells-1); i++)
      {
        int index=getVectorIndex(i,j,k);

        //Set up b=div(w3)
        m_storeB.at(index)=mathFunction::calcDivergenceVec3(&m_newVelocityField,i,j,k, m_cellSize, m_noCells);

        //Add divergence from explosion setup or from particle interaction
        m_storeB.at(index)+=m_divergenceField.at(index);
      }
    }
  }

  //Set up A
//  float Aii=((-1.0)*6.0)/(pow(m_cellSize,2));
  float Aii=(6.0)/(pow(m_cellSize,2));
//  float Aii=6.0;
  float Aij=1.0/(pow(m_cellSize,2));
//  float Aij=1.0;

  //Initial guess is zero so initial field is b=div(w3)
//  m_storeFieldFloat=m_storeB;
  m_storeFieldFloat=m_pressureField;


  //Set boundary values for this guess
  setBoundaryValuesFloat(&m_storeFieldFloat);

  //Send to linear solver to find p
  mathFunction::linearSystemSolveFloat(&m_pressureField, &m_storeFieldFloat, &m_storeB, Aii, Aij, iterations, m_noCells, 0.0);


  //From new pressure field calculate new velocity
  ngl::Vec3 divP;

  for (int k=1; k<(m_noCells-1); k++)
  {
    for (int j=1; j<(m_noCells-1); j++)
    {
      for (int i=1; i<(m_noCells-1); i++)
      {
        int index=getVectorIndex(i,j,k);

        //Calulate div(p)
        divP=mathFunction::calcDivergenceFloat(&m_pressureField, i, j, k, m_cellSize, m_noCells);

        //Calculate new velocity w4=w3-div(p)
        m_newVelocityField.at(index)=m_newVelocityField.at(index)-divP;

      }
    }
  }

  setBoundaryVelocity();

}


void Grid::calculateVorticity()
{
  //Set up variables
  ngl::Vec3 curlGradient;
  float curlGradientLength;
  ngl::Vec3 curlNormal;
  ngl::Vec3 vorticityForce;

  //Calculate curl vectors, W, and magnitudes, magW
  mathFunction::calculateCurl(&m_newVelocityField, &m_curlVectorField, &m_curlMagnitudes, m_noCells, m_cellSize);
  setBoundaryValuesVec3(&m_curlVectorField);
  setBoundaryValuesFloat(&m_curlMagnitudes);

  for (int k=1; k<(m_noCells-1); k++)
  {
    for (int j=1; j<(m_noCells-1); j++)
    {
      for (int i=1; i<(m_noCells-1); i++)
      {
        int index=getVectorIndex(i,j,k);

        //Calculate gradient of curl, div(magW)
        curlGradient=mathFunction::calcDivergenceFloat(&m_curlMagnitudes, i, j, k, m_cellSize, m_noCells);

        //Normalise to find gradient normal, N
        curlGradientLength=curlGradient.length() + 0.0001; //NB! Add 0.0001 to avoid issues when dividing by zero
        curlNormal=(1/curlGradientLength)*curlGradient;

        //Calulate vorticity force: f=constant*cellSize*(NxW)
        vorticityForce.m_x=(curlNormal.m_y*curlGradient.m_z)-(curlNormal.m_z-curlGradient.m_y);
        vorticityForce.m_y=(curlNormal.m_z*curlGradient.m_x)-(curlNormal.m_x-curlGradient.m_z);
        vorticityForce.m_x=(curlNormal.m_x*curlGradient.m_y)-(curlNormal.m_y-curlGradient.m_x);

        vorticityForce*=(m_vorticityConstant*m_cellSize);

        //Add to force field
        m_forceField.at(index)+=vorticityForce;

      }
    }
  }
}

void Grid::calculateBuoyancy()
{
  //Applies force only in up direction, ie. y-direction
  //Force=buoyancyConst*(T_at_index-T_ambient)

  float forceY;

  for (int k=1; k<(m_noCells-1); k++)
  {
    for (int j=1; j<(m_noCells-1); j++)
    {
      for (int i=1; i<(m_noCells-1); i++)
      {
        int index=getVectorIndex(i,j,k);
        forceY=m_buoyancyConstant*(m_newTemperatureField.at(index)-m_ambientTemp);
        m_forceField.at(index).m_y+=forceY;
      }
    }
  }
}

void Grid::updateTemperatureField(float _dt)
{
  //Advect temperature
  advectTemperature(_dt);

  //Diffuse
  diffuseTemperature(_dt);

  //Set boundaries
  setBoundaryTemperature();

  //Cool particles
  dissipateTemperature();


}


void Grid::advectTemperature(float _dt)
{
  //Need to loop over 3D space i,j,k
  for (int k=0; k<m_noCells; k++)
  {
    for (int j=0; j<m_noCells; j++)
    {
      for (int i=0; i<m_noCells; i++)
      {
        //Find previous position for fictive particle at (i,j,k)
        //Get old_velocity at (i,j,k)
        int index=getVectorIndex(i,j,k);

        ngl::Vec3 oldVelocity=m_oldVelocityField.at(index);

        //Get actual position for (i,j,k)
        ngl::Vec3 currPosition;
        currPosition.m_x=(i+0.5)*m_cellSize+m_origin.m_x;
        currPosition.m_y=(j+0.5)*m_cellSize+m_origin.m_y;
        currPosition.m_z=(k+0.5)*m_cellSize+m_origin.m_z;

        //Use these values in RK2 to find x0
        ngl::Vec3 prevPosition=mathFunction::RK2_integrator(currPosition,oldVelocity,(-1.0)*_dt);

        //Determine the cell that is in
        //Is cell inside boundaries? If not find nearest boundary cell and use that value for advection.
        float advectTemperature=getTemperatureFromField(prevPosition);

        //Need to store new values temporarily so don't interfer with the w1 we're collecting values from
        m_storeFieldFloat.at(index)=advectTemperature;

      }
    }
  }

  //When all advection velocities calculated, set the temporarily stored field to the newVelocity field.
  m_newTemperatureField=m_storeFieldFloat;

  setBoundaryTemperature();

}

void Grid::diffuseTemperature(float _dt)
{
  //Set number of iterations
  int iterations=30;

  //Calculate Aij and Aii
  float Aij;
  float Aii;

  ///Might need to divide by 1/m_cellSize^2, like velocity projection
  Aij=_dt*m_thermalConductivity*pow(m_cellSize,2);

  Aii=1.0+(6.0*Aij);

  //Set up b in AT=b
  m_storeB=m_newTemperatureField;

  //Set initial guess to previous temperature field?
  m_storeFieldFloat=m_oldTemperatureField;

  //Use linear solver to find new temp
  mathFunction::linearSystemSolveFloat(&m_newTemperatureField, &m_storeFieldFloat, &m_storeB, Aii, Aij, iterations, m_noCells, m_ambientTemp);

  //Set boundary temp
  setBoundaryTemperature();
}

void Grid::dissipateTemperature()
{
  ///To do: m_tempMax should change as fluid cools down. Would need a function to find largest value in vector

  float overallDifference=m_maxTemp-m_ambientTemp;
  float currentDifference;
  float ratioTempDifference;

  float coolingContribution;

  for (int i=0; i<pow(m_noCells,3); i++)
  {
    currentDifference=m_newTemperatureField.at(i)-m_ambientTemp;

    //RatioTempDiff=(T_at_index-T_ambient)/(T_max-T_ambient)
    ratioTempDifference=currentDifference/overallDifference;

    //Calculate cooling contribution=coolConst*(ratioTempDiff)^4
    coolingContribution=m_coolingConstant*(pow(ratioTempDifference,4));

    //Calculate new temperature after cooling
    m_newTemperatureField.at(i)-=coolingContribution;


    //Need to make sure that temp does not go below ambient temp
    if (m_newTemperatureField.at(i)<m_ambientTemp)
    {
      m_newTemperatureField.at(i)=m_ambientTemp;
    }
  }
}

ngl::Vec3 Grid::getVelocityFromField(ngl::Vec3 _particlePosition)
{

  ///To do: Choose field to take interpolation from?

  ngl::Vec3 fieldVelocity;

  //Find grid indices from particle position.
  ngl::Vec3 indexParticle=(1/m_cellSize)*(_particlePosition-m_origin);

  //Find which cell the particle is in and set this to index0. Are there any issues with using round?
  int index0_X;
  int index0_Y;
  int index0_Z;
  index0_X=trunc(indexParticle.m_x);
  index0_Y=trunc(indexParticle.m_y);
  index0_Z=trunc(indexParticle.m_z);


  //Check if particle inside boundaries
  if (index0_X>0 && index0_X<(m_noCells-1) && index0_Y>0 && index0_Y<(m_noCells-1) && index0_Z>0 && index0_Z<(m_noCells-1))
  {
    //Find the index of the nearest neighbour cell in i,j,k directions
    int index1_X=index0_X-1;
    int index1_Y=index0_Y-1;
    int index1_Z=index0_Z-1;

    if ((indexParticle.m_x-index0_X)>=0.5)
    {
      index1_X=index0_X+1;
    }
    if ((indexParticle.m_y-index0_Y)>=0.5)
    {
      index1_Y=index0_Y+1;
    }
    if ((indexParticle.m_z-index0_Z)>=0.5)
    {
      index1_Z=index0_Z+1;
    }

    //  //Test routine for finding indices
    //  std::cout<<"Position particle: ["<<_particlePosition.m_x<<","<<_particlePosition.m_y<<","<<_particlePosition.m_z<<"]\n";
    //  std::cout<<"Index particle: ["<<indexParticle.m_x<<","<<indexParticle.m_y<<","<<indexParticle.m_z<<"]\n";
    //  std::cout<<"Index in Cell: ["<<index0_X<<","<<index0_Y<<","<<index0_Z<<"]\n";
    //  std::cout<<"Index neigbour cell: ["<<index1_X<<","<<index1_Y<<","<<index1_Z<<"]\n";

    //Use trilinear interpolation to get field velocity based on position
    fieldVelocity=mathFunction::trilinearInterpVec3(&m_newVelocityField, index0_X, index0_Y, index0_Z, index1_X, index1_Y, index1_Z, indexParticle, m_noCells);
  }

  //If either in boundary cell or outside bounding box
  else
  {
    //If outside bounding box by i,j or k being smaller than 0 or greater than m_noCells-1, set them to 0 or m_noCells-1.
    //If inside boundary cell then will not change indices and instead just find velocity inside this cell.

    //If index0_X<0 set it to 0
    if (index0_X<0)
    {
      index0_X=0;
    }
    //If index0_X>m_noCells-1, set it to m_noCells-1
    if (index0_X>(m_noCells-1))
    {
      index0_X=(m_noCells-1);
    }
    //Do same for Y and Z
    if (index0_Y<0)
    {
      index0_Y=0;
    }
    if (index0_Y>(m_noCells-1))
    {
      index0_Y=(m_noCells-1);
    }

    if (index0_Z<0)
    {
      index0_Z=0;
    }
    //If index0_X>m_noCells-1, set it to m_noCells-1
    if (index0_Z>(m_noCells-1))
    {
      index0_Z=(m_noCells-1);
    }


    //Find velocity at new index
    fieldVelocity=m_newVelocityField.at(getVectorIndex(index0_X,index0_Y,index0_Z));

  }


  //Return velocity
  return fieldVelocity;
}

float Grid::getTemperatureFromField(ngl::Vec3 _particlePosition)
{

  ///To do: Choose field to take interpolation from?

  float fieldTemperature;

  //Find grid indices from particle position.
  ngl::Vec3 indexParticle=(1/m_cellSize)*(_particlePosition-m_origin);

  //Find which cell the particle is in and set this to index0. Are there any issues with using round?
  int index0_X;
  int index0_Y;
  int index0_Z;
  index0_X=trunc(indexParticle.m_x);
  index0_Y=trunc(indexParticle.m_y);
  index0_Z=trunc(indexParticle.m_z);


  //Check if particle inside boundaries
  if (index0_X>0 && index0_X<(m_noCells-1) && index0_Y>0 && index0_Y<(m_noCells-1) && index0_Z>0 && index0_Z<(m_noCells-1))
  {
    //Find the index of the nearest neighbour cell in i,j,k directions
    int index1_X=index0_X-1;
    int index1_Y=index0_Y-1;
    int index1_Z=index0_Z-1;

    if ((indexParticle.m_x-index0_X)>=0.5)
    {
      index1_X=index0_X+1;
    }
    if ((indexParticle.m_y-index0_Y)>=0.5)
    {
      index1_Y=index0_Y+1;
    }
    if ((indexParticle.m_z-index0_Z)>=0.5)
    {
      index1_Z=index0_Z+1;
    }

    //  //Test routine for finding indices
    //  std::cout<<"Position particle: ["<<_particlePosition.m_x<<","<<_particlePosition.m_y<<","<<_particlePosition.m_z<<"]\n";
    //  std::cout<<"Index particle: ["<<indexParticle.m_x<<","<<indexParticle.m_y<<","<<indexParticle.m_z<<"]\n";
    //  std::cout<<"Index in Cell: ["<<index0_X<<","<<index0_Y<<","<<index0_Z<<"]\n";
    //  std::cout<<"Index neigbour cell: ["<<index1_X<<","<<index1_Y<<","<<index1_Z<<"]\n";

    //Use trilinear interpolation to get field velocity based on position
    fieldTemperature=mathFunction::trilinearInterpFloat(&m_newTemperatureField, index0_X, index0_Y, index0_Z, index1_X, index1_Y, index1_Z, indexParticle, m_noCells);
  }

  //If either in boundary cell or outside bounding box
  else
  {
    //If outside bounding box by i,j or k being smaller than 0 or greater than m_noCells-1, set them to 0 or m_noCells-1.
    //If inside boundary cell then will not change indices and instead just find velocity inside this cell.

    //If index0_X<0 set it to 0
    if (index0_X<0)
    {
      index0_X=0;
    }
    //If index0_X>m_noCells-1, set it to m_noCells-1
    if (index0_X>(m_noCells-1))
    {
      index0_X=(m_noCells-1);
    }
    //Do same for Y and Z
    if (index0_Y<0)
    {
      index0_Y=0;
    }
    if (index0_Y>(m_noCells-1))
    {
      index0_Y=(m_noCells-1);
    }

    if (index0_Z<0)
    {
      index0_Z=0;
    }
    //If index0_X>m_noCells-1, set it to m_noCells-1
    if (index0_Z>(m_noCells-1))
    {
      index0_Z=(m_noCells-1);
    }


    //Find velocity at new index
    fieldTemperature=m_newTemperatureField.at(getVectorIndex(index0_X,index0_Y,index0_Z));

  }


  //Return velocity
  return fieldTemperature;
}
