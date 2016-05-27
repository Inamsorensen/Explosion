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

  //Set space for fields
  for (int i=0; i<pow(m_noCells,3); i++)
  {
    m_pressureField.push_back(0.0);
    m_forceField.push_back(ngl::Vec3(0.0,0.0,0.0));

    m_newVelocityField.push_back(ngl::Vec3(0.0,0.0,0.0));

    m_newTemperatureField.push_back(0.0);

    m_storeFieldVec3.push_back(ngl::Vec3(0.0,0.0,0.0));
    m_storeFieldFloat.push_back(0.0);
    m_storeZeroFieldVec3.push_back(ngl::Vec3(0.0,0.0,0.0));
    m_storeZeroFieldFloat.push_back(0.0);

    m_curlVectorField.push_back(ngl::Vec3(0.0,0.0,0.0));
    m_curlMagnitudes.push_back(0.0);

    m_divergenceField.push_back(0.0);
    m_storeB.push_back(0.0);

  }

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

  //Create noise in velocity with random values
  ngl::Random* rand=ngl::Random::instance();
  float noiseValX;
  float noiseValY;
  float noiseValZ;
  ngl::Vec3 velocity;

  for (int i=0; i<pow(m_noCells,3); i++)
  {
    rand->setSeed(i);
    noiseValX=rand->randomNumber(m_noiseConstant);
    rand->setSeed(i+1000);
    noiseValY=rand->randomNumber(m_noiseConstant);
    rand->setSeed(i+2000);
    noiseValZ=rand->randomNumber(m_noiseConstant);
    velocity.m_x=noiseValX;
    velocity.m_y=noiseValY;
    velocity.m_z=noiseValZ;

    m_newVelocityField.at(i)=velocity;
  }

  //For solid boundary, set normal velocity to zero
//  setBoundaryVelocity();

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

    //Check that the index is inside grid, if not throw error
  if (indexOrigin_X<1 || indexOrigin_Y<1 || indexOrigin_Z<1 || indexOrigin_X>(m_noCells-1) || indexOrigin_Y>(m_noCells-1) || indexOrigin_Z>(m_noCells-1))
  {
    std::cout<<"ExplosionOrigin: ["<<indexOrigin_X<<" "<<indexOrigin_Y<<" "<<indexOrigin_Z<<"]\n";
    throw std::invalid_argument("Origin of explosion is in a boundary cell or outside the grid");
  }

    //Check radius against cellSize to see how many cells will be changed in one direction
  int noCellsInRadius=trunc(_radius/m_cellSize);

  //Check that explosion isn't too big compared to grid
  if (_radius>(((m_noCells-2.0)/2.0)*m_cellSize))
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


int Grid::getVectorIndex(int i, int j, int k)
{
  int vectorIndex=i+(m_noCells*j)+(m_noCells*m_noCells*k);
  return vectorIndex;
}

void Grid::setBoundaryValuesVec3(std::vector<ngl::Vec3> *_vec3Field)
{
  for (int j=0; j<m_noCells; j++)
  {
    for (int i=0; i<m_noCells; i++)
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
  for (int j=0; j<m_noCells; j++)
  {
    for (int i=0; i<m_noCells; i++)
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
  for (int j=0; j<m_noCells; j++)
  {
    for (int i=0; i<m_noCells; i++)
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
  for (int j=0; j<m_noCells; j++)
  {
    for (int i=0; i<m_noCells; i++)
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

void Grid::update(float _dt)
{

  updateVelocityField(_dt);

  updateTemperatureField(_dt);

}


void Grid::updateVelocityField(float _dt)
{
  //Add external forces. w1=w0+dt*f
  addForce(_dt);

  //Project Div^2p=Div.w3, w4=w3-Dp
  ///Reference code has another project in it. Makes the results slightly better, but can be commented out
  project(_dt);

  //Advect. w2=w1(backTrace(x,-dt))
  advectVelocity(_dt);

  //Project Div^2p=Div.w3, w4=w3-Divp
  project(_dt);

  //Set divergence field to zero before new update
  m_divergenceField=m_storeZeroFieldFloat;

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
  //Forces only affect the inside cells, so boundary values are unchanged
  for (int i=0; i<pow(m_noCells,3); i++)
  {
    m_newVelocityField.at(i)+=_dt*m_forceField.at(i);
  }

}

void Grid::calculateVorticity()
{
  //Set up variables
  ngl::Vec3 curlVector;
  ngl::Vec3 curlGradient;
  float curlGradientLength;
  ngl::Vec3 curlNormal;
  ngl::Vec3 vorticityForce;

  //Calculate curl vectors, W, and magnitudes, magW
  mathFunction::calculateCurl(&m_newVelocityField, &m_curlVectorField, &m_curlMagnitudes, m_noCells, m_cellSize);

//  setBoundaryValuesVec3(&m_curlVectorField);
//  setBoundaryValuesFloat(&m_curlMagnitudes);

  for (int k=1; k<(m_noCells-1); k++)
  {
    for (int j=1; j<(m_noCells-1); j++)
    {
      for (int i=1; i<(m_noCells-1); i++)
      {
        int index=getVectorIndex(i,j,k);

        //Calculate gradient of curl magnitude, div(magW)
        curlGradient=mathFunction::calcDivergenceFloat(&m_curlMagnitudes, i, j, k, m_cellSize, m_noCells);

        //Normalise to find gradient normal, N
        curlGradientLength=curlGradient.length() + 0.0001; //NB! Add 0.0001 to avoid issues when dividing by zero
        curlNormal=(1.0/curlGradientLength)*curlGradient;

        //Calulate vorticity force: f=constant*cellSize*(NxW)
        curlVector=m_curlVectorField.at(index);
        vorticityForce.m_x=(curlNormal.m_y*curlVector.m_z)-(curlNormal.m_z*curlVector.m_y);
        vorticityForce.m_y=(curlNormal.m_z*curlVector.m_x)-(curlNormal.m_x*curlVector.m_z);
        vorticityForce.m_x=(curlNormal.m_x*curlVector.m_y)-(curlNormal.m_y*curlVector.m_x);

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


void Grid::advectVelocity(float _dt)
{

///Method only advects non-boundary cells and uses current velocity field for back trace.

  //Set field to store calculated velocity values
  m_storeFieldVec3=m_newVelocityField;

  //Need to loop over non-boundary cells
  for (int k=1; k<(m_noCells-1); k++)
  {
    for (int j=1; j<(m_noCells-1); j++)
    {
      for (int i=1; i<(m_noCells-1); i++)
      {
        //Find previous position for fictive particle at (i,j,k)
          //Get velocity at (i,j,k)
        int index=getVectorIndex(i,j,k);
        ngl::Vec3 velocity=m_newVelocityField.at(index);

          //Get actual position for (i,j,k)
        ngl::Vec3 currPosition;
        float fi=(float)i;
        float fj=(float)j;
        float fk=(float)k;
        currPosition.m_x=(fi*m_cellSize)+m_origin.m_x;
        currPosition.m_y=(fj*m_cellSize)+m_origin.m_y;
        currPosition.m_z=(fk*m_cellSize)+m_origin.m_z;

          //Use these values in RK2 to find previous position
        ngl::Vec3 prevPosition=mathFunction::RK2_integrator(currPosition,velocity,(-1.0)*_dt);

        //Get velocity at that position.
        //getVelocityFromField checks if within boundaries and sets what to do if not
        ngl::Vec3 advectVelocity=getVelocityFromField(prevPosition);

        //Need to store new values temporarily so don't interfer with the w1 we're collecting values from
        m_storeFieldVec3.at(index)=advectVelocity;

      }
    }
  }

  //When all advection velocities calculated, set the temporarily stored field to the newVelocity field.
  m_newVelocityField=m_storeFieldVec3;

}


void Grid::project(float _dt)
{
  ///To do:
  /// Add check for convergence in loop
  /// Set proper boundary values for p
  ///
  /// Unsure about whether the values used for A and b in Ap=b are correct. Other possibilities commented out.
  /// Current choice based on values used by J. Stam [1999]

  //Set number of iterations
  int iterations=30;
//  int iterations=20;

  //Set up b. Same for each iteration so only do this once
  for (int k=1; k<(m_noCells-1); k++)
  {
    for (int j=1; j<(m_noCells-1); j++)
    {
      for (int i=1; i<(m_noCells-1); i++)
      {
        int index=getVectorIndex(i,j,k);

        //Set up b=div(w3)
        /// Not sure which of the following should be used, the one in use is likely the correct one
//        m_storeB.at(index)=(-1.0)*mathFunction::calcDivergenceVec3(&m_newVelocityField,i,j,k, m_cellSize, m_noCells);
//        m_storeB.at(index)=(1.0/_dt)*mathFunction::calcDivergenceVec3(&m_newVelocityField,i,j,k, m_cellSize, m_noCells);
         m_storeB.at(index)=mathFunction::calcDivergenceVec3(&m_newVelocityField,i,j,k, m_cellSize, m_noCells);

        //Add divergence from explosion setup or from particle interaction
//        m_storeB.at(index)+=m_divergenceField.at(index);
        m_storeB.at(index)-=m_divergenceField.at(index);
      }
    }
  }

  //Set up A
  float Aii=((-1.0)*6.0)/(pow(m_cellSize,2)); //Weird
//  float Aii=(6.0)/(pow(m_cellSize,2)); //Good
//  float Aii=6.0;
  float Aij=1.0/(pow(m_cellSize,2)); //Good
//  float Aij=-1.0/(pow(m_cellSize,2)); //Weird
//  float Aij=1.0;

  //Set initial guess for pressure to previous pressure values
  m_storeFieldFloat=m_pressureField;

  //Set boundary values
  setBoundaryValuesFloat(&m_storeFieldFloat);
  setBoundaryValuesFloat(&m_storeB);

  //Send to linear solver to find p
  mathFunction::linearSystemSolveFloat(&m_pressureField, &m_storeFieldFloat, &m_storeB, Aii, Aij, iterations, m_noCells, 0.0);

  //Set boundary values for pressure field equal to adjacent cell
  setBoundaryValuesFloat(&m_pressureField);

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

        //Calculate new velocity w4=w3-dt*div(p)
        m_newVelocityField.at(index)=m_newVelocityField.at(index)-(_dt*divP);

      }
    }
  }

//  setBoundaryValuesVec3(&m_newVelocityField);
  setBoundaryVelocity();

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
  //Need to loop over non-boundary cells
  for (int k=1; k<(m_noCells-1); k++)
  {
    for (int j=1; j<(m_noCells-1); j++)
    {
      for (int i=1; i<(m_noCells-1); i++)
      {
        //Find previous position for fictive particle at (i,j,k)
        //Get velocity at (i,j,k)
        int index=getVectorIndex(i,j,k);
        ngl::Vec3 velocity=m_newVelocityField.at(index);

        //Get actual position for (i,j,k)
        ngl::Vec3 currPosition;
        currPosition.m_x=(((float)i)*m_cellSize)+m_origin.m_x;
        currPosition.m_y=(((float)j)*m_cellSize)+m_origin.m_y;
        currPosition.m_z=(((float)k)*m_cellSize)+m_origin.m_z;

        //Use these values in RK2 to find previous position
        ngl::Vec3 prevPosition=mathFunction::RK2_integrator(currPosition,velocity,(-1.0)*_dt);

        //Get temperature at this previous position
        float advectTemperature=getTemperatureFromField(prevPosition);

        //Need to store new values temporarily so don't interfer with the w1 we're collecting values from
        m_storeFieldFloat.at(index)=advectTemperature;

      }
    }
  }

  //When all advection velocities calculated, set the temporarily stored field to the newVelocity field.
  m_newTemperatureField=m_storeFieldFloat;

//  setBoundaryTemperature();

}


void Grid::diffuseTemperature(float _dt)
{
  //Set number of iterations
  int iterations=30;

  //Calculate Aij and Aii
  float Aij;
  float Aii;

  Aij=(-1.0)*_dt*m_thermalConductivity/pow(m_cellSize,2);

  Aii=1.0+(6.0*(-1.0)*Aij);

  //Set up b in AT=b
  m_storeB=m_newTemperatureField;

  //Set initial guess to current temperature field
  m_storeFieldFloat=m_newTemperatureField;

  //Use linear solver to find new temp
  mathFunction::linearSystemSolveFloat(&m_newTemperatureField, &m_storeFieldFloat, &m_storeB, Aii, Aij, iterations, m_noCells, m_ambientTemp);

}


void Grid::dissipateTemperature()
{
  ///To do:
  /// m_tempMax should change as fluid cools down. Would need a function to find largest value in vector

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

  ngl::Vec3 fieldVelocity;

  //Find grid indices from particle position.
  ngl::Vec3 indexParticle=(1.0/m_cellSize)*(_particlePosition-m_origin);

  //Find which cell the particle is in and set this to index0. Are there any issues with using round?
  int index0_X;
  int index0_Y;
  int index0_Z;
  index0_X=floor(indexParticle.m_x);
  index0_Y=floor(indexParticle.m_y);
  index0_Z=floor(indexParticle.m_z);


  //Check if particle inside boundaries
  if (index0_X>0 && index0_X<(m_noCells-1) && index0_Y>0 && index0_Y<(m_noCells-1) && index0_Z>0 && index0_Z<(m_noCells-1))
  {

    ///Method for index setup so matches MAC grid
    /// index1 is always the cell with i,j,k one higher than cell at index0
    int index1_X=index0_X+1;
    int index1_Y=index0_Y+1;
    int index1_Z=index0_Z+1;

    //Use trilinear interpolation to get field velocity based on position
    fieldVelocity=mathFunction::trilinearInterpVec3(&m_newVelocityField, index0_X, index0_Y, index0_Z, index1_X, index1_Y, index1_Z, indexParticle, m_noCells);
  }

  //If either in boundary cell or outside bounding box, return velocity of closest boundary cell
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

  ///To do:
  /// Set return from boundary/outside to ambient temperature

  float fieldTemperature;

  //Find grid indices from particle position.
  ngl::Vec3 indexParticle=(1.0/m_cellSize)*(_particlePosition-m_origin);

  //Find which cell the particle is in and set this to index0.
  int index0_X;
  int index0_Y;
  int index0_Z;
  index0_X=floor(indexParticle.m_x);
  index0_Y=floor(indexParticle.m_y);
  index0_Z=floor(indexParticle.m_z);


  //Check if particle inside boundaries
  if (index0_X>0 && index0_X<(m_noCells-1) && index0_Y>0 && index0_Y<(m_noCells-1) && index0_Z>0 && index0_Z<(m_noCells-1))
  {
    ///Method for index setup so matches MAC grid
    /// index1 is always the cell with i,j,k one higher than cell at index0
    int index1_X=index0_X+1;
    int index1_Y=index0_Y+1;
    int index1_Z=index0_Z+1;

    //Use trilinear interpolation to get field velocity based on position
    fieldTemperature=mathFunction::trilinearInterpFloat(&m_newTemperatureField, index0_X, index0_Y, index0_Z, index1_X, index1_Y, index1_Z, indexParticle, m_noCells);
  }

  //If either in boundary cell or outside bounding box
  else
  {
    //If outside bounding box by i,j or k being smaller than 0 or greater than m_noCells-1, set them to 0 or m_noCells-1.
    //If inside boundary cell then will not change indices and instead just find temperature inside this cell.

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


void Grid::addTemperatureFromParticle(int _index, float _temperature)
{
  m_newTemperatureField.at(_index)+=_temperature;
}


void Grid::addDivergenceFromParticle(int _index, float _divergence)
{
  m_divergenceField.at(_index)+=_divergence;
}
