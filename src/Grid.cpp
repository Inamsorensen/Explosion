#include "Grid.h"

#include <stdexcept>
#include <cmath>
#include <math.h>
#include <iostream>

Grid* Grid::m_instance=nullptr;

Grid::Grid(ngl::Vec3 _origin, float _gridSize, int _noCells)
{
  m_origin=_origin;
  m_gridSize=_gridSize;
  m_noCells=_noCells;

  m_cellSize=m_gridSize/m_noCells;

  //Set up zero fields for velocity, pressure and force
  for (int i=0; i<pow(m_noCells,3); i++)
  {
    m_oldVelocityField.push_back(ngl::Vec3(0.0,0.0,0.0));
    m_pressureField.push_back(0.0);
    m_forceField.push_back(ngl::Vec3(0.0,0.0,0.0));

    //Set space for store fields
    m_storeFieldVec3.push_back(ngl::Vec3(0.0,0.0,0.0));
    m_storeFieldFloat.push_back(0.0);

    //Set space for storeB and storeA. They only contain values for cells that aren't solid. A=[(m_noCells-2)^3,(m_noCells-2)^3]
    m_storeB.push_back(0.0);

    for (int j=0; j<pow(m_noCells,3); j++)
    {
      m_storeA.push_back(0.0);
    }
  }
  m_newVelocityField=m_oldVelocityField;


  ///Need to set boundary values if change from zero

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


void Grid::setVelocityField(std::vector<ngl::Vec3> _velocityField)
{
  m_newVelocityField=_velocityField;
  setBoundaryValuesVec3(&m_newVelocityField);
  setBoundaryVelocity();
}

void Grid::setForceField(std::vector<ngl::Vec3> _forceField)
{
  m_forceField=_forceField;
  setBoundaryValuesVec3(&m_forceField);
}

void Grid::setPressureField(std::vector<float> _pressureField)
{
  m_pressureField=_pressureField;
  setBoundaryValuesFloat(&m_pressureField);
}

void Grid::update(float _dt)
{
  updateVelocityField(_dt);

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




void Grid::swap()
{
  m_oldVelocityField=m_newVelocityField;
}


void Grid::updateVelocityField(float _dt)
{
  //Set old velocity field equal to new velocity field so new velocity field can be given new values
  swap();

  //Add external forces. w1=w0+dt*f
  addForce(_dt);

  //Advect. w2=w1(p(x,-dt))
  advect(_dt);

  //Diffuse (I-vdtD2)w3=w2

  //Project D2p=D.w3, w4=w3-Dp
  project();
}


void Grid::addForce(float _dt)
{
  for (int i=0; i<pow(m_noCells,3); i++)
  {
    m_newVelocityField.at(i)=m_oldVelocityField.at(i)+(m_forceField.at(i)*_dt);
  }

  setBoundaryValuesVec3(&m_newVelocityField);
  setBoundaryVelocity();

}



void Grid::advect(float _dt)
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
        ngl::Vec3 prevPosition=RK2_integrator(currPosition,oldVelocity,(-1)*_dt);

        //Determine the cell that is in
        //Is cell inside boundaries? If not find nearest boundary cell and use that value for advection.
        ngl::Vec3 advectVelocity=getVelocityFromField(prevPosition);

        //Need to store new values temporarily so don't interfer with the w1 we're collecting values from
        m_storeFieldVec3.at(index)=advectVelocity;

      }
    }
  }

  //When all advection velocities calculated, set the temporarily stored field to the newVelocity field.
  m_newVelocityField=m_storeFieldVec3;

  setBoundaryValuesVec3(&m_newVelocityField);
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

  //Set up b. Same for each iteration so only do this once
  for (int k=1; k<(m_noCells-1); k++)
  {
    for (int j=1; j<(m_noCells-1); j++)
    {
      for (int i=1; i<(m_noCells-1); i++)
      {
        int index=getVectorIndex(i,j,k);

        //Set up b=div(w3)
        m_storeB.at(index)=calcDivergenceVec3(&m_newVelocityField,i,j,k);
      }
    }
  }

  //Set up A
  float Aii=((-1)*6.0)/(pow(m_cellSize,2));
  float Aij=1.0/(pow(m_cellSize,2));

  //Initial guess is zero so initial field is b=div(w3)
  m_storeFieldFloat=m_storeB;

  //Set boundary values for this guess
  setBoundaryValuesFloat(&m_storeFieldFloat);

  //Send to linear solver to find p
  linearSystemSolve(&m_pressureField, &m_storeFieldFloat, &m_storeB, Aii, Aij, iterations);


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
        divP=calcDivergenceFloat(&m_pressureField, i, j, k);

        //Calculate new velocity w4=w3-div(p)
        m_newVelocityField.at(index)=m_newVelocityField.at(index)-divP;

      }
    }
  }

  setBoundaryVelocity();

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
    fieldVelocity=trilinearInterpVec3(&m_newVelocityField, index0_X, index0_Y, index0_Z, index1_X, index1_Y, index1_Z, indexParticle);
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
