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

//  //Set up zero fields for velocity, pressure and force
//  for (int i=0; i<pow(m_noCells,3); i++)
//  {
//    m_oldVelocityField.push_back(ngl::Vec3(0.0,0.0,0.0));
//    m_pressureField.push_back(0.0);
//    m_forceField.push_back(ngl::Vec3(0.0,0.0,0.0));

//    //Set all cells to not be boundary cells
//    m_solidBoundary.push_back(false);
//  }
//  m_newVelocityField=m_oldVelocityField;

//  //Set solid boundary cells
//  for (int i=0; i<m_noCells; i++)
//  {
//    for (int j=0; j<m_noCells; j++)
//    {
//      m_solidBoundary.at(getVectorIndex(0,i,j))=true;
//      m_solidBoundary.at(getVectorIndex(m_noCells-1,i,j))=true;
//      m_solidBoundary.at(getVectorIndex(i,0,j))=true;
//      m_solidBoundary.at(getVectorIndex(i,m_noCells-1,j))=true;
//      m_solidBoundary.at(getVectorIndex(i,j,0))=true;
//      m_solidBoundary.at(getVectorIndex(i,j,m_noCells-1))=true;
//    }
//  }

  //Create cells and set zero fields for velocity, pressure and force. Also set all cells to fluid cells
  for (int i=0; i<pow(m_noCells,3); i++)
  {
    Cell* newCell=new Cell;
    newCell->m_isSolid=false;
    newCell->m_oldVelocity=ngl::Vec3(0.0,0.0,0.0);
    newCell->m_newVelocity=newCell->m_oldVelocity;
    newCell->m_externalForce=ngl::Vec3(0.0,0.0,0.0);
    newCell->m_pressure=0.0;
    m_listCells.push_back(newCell);
  }
  //Set solid boundary cells
  for (int i=0; i<m_noCells; i++)
  {
    for (int j=0; j<m_noCells; j++)
    {
      m_listCells.at(getVectorIndex(0,i,j))->m_isSolid=true;
      m_listCells.at(getVectorIndex(m_noCells-1,i,j))->m_isSolid=true;
      m_listCells.at(getVectorIndex(i,0,j))->m_isSolid=true;
      m_listCells.at(getVectorIndex(i,m_noCells-1,j))->m_isSolid=true;
      m_listCells.at(getVectorIndex(i,j,0))->m_isSolid=true;
      m_listCells.at(getVectorIndex(i,j,m_noCells-1))->m_isSolid=true;
    }
  }

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
//  m_newVelocityField=_velocityField;
//  setBoundaryValuesVec3(&m_newVelocityField);
//  setBoundaryVelocity();

  for (int i=0; i<pow(m_noCells,3); i++)
  {
    m_listCells.at(i)->m_newVelocity=_velocityField.at(i);
  }
  setBoundaryValuesVec3("velocity");
  setBoundaryVelocity();
}

void Grid::setForceField(std::vector<ngl::Vec3> _forceField)
{
//  m_forceField=_forceField;
//  setBoundaryValuesVec3(&m_forceField);

  for (int i=0; i<pow(m_noCells,3); i++)
  {
    m_listCells.at(i)->m_externalForce=_forceField.at(i);
  }
  setBoundaryValuesVec3("force");
  setBoundaryVelocity();
}

void Grid::setPressureField(std::vector<float> _pressureField)
{
//  m_pressureField=_pressureField;
//  setBoundaryValuesFloat(&m_pressureField);

  for (int i=0; i<pow(m_noCells,3); i++)
  {
    m_listCells.at(i)->m_pressure=_pressureField.at(i);
  }
  setBoundaryValuesFloat("pressure");
  setBoundaryVelocity();
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

void Grid::setBoundaryValuesVec3(std::string _fieldType)
{
  for (int i=0; i<m_noCells; i++)
  {
    for (int j=0; j<m_noCells; j++)
    {
//      _vec3Field->at(getVectorIndex(0,i,j))=_vec3Field->at(getVectorIndex(1,i,j));
//      _vec3Field->at(getVectorIndex(m_noCells-1,i,j))=_vec3Field->at(getVectorIndex(m_noCells-2,i,j));
//      _vec3Field->at(getVectorIndex(i,0,j))=_vec3Field->at(getVectorIndex(i,1,j));
//      _vec3Field->at(getVectorIndex(i,m_noCells-1,j))=_vec3Field->at(getVectorIndex(i,m_noCells-2,j));
//      _vec3Field->at(getVectorIndex(i,j,0))=_vec3Field->at(getVectorIndex(i,j,1));
//      _vec3Field->at(getVectorIndex(i,j,m_noCells-1))=_vec3Field->at(getVectorIndex(i,j,m_noCells-2));

      if(_fieldType=="velocity")
      {
        m_listCells.at(getVectorIndex(0,i,j))->m_newVelocity=m_listCells.at(getVectorIndex(1,i,j))->m_newVelocity;
        m_listCells.at(getVectorIndex(m_noCells-1,i,j))->m_newVelocity=m_listCells.at(getVectorIndex(m_noCells-2,i,j))->m_newVelocity;
        m_listCells.at(getVectorIndex(i,0,j))->m_newVelocity=m_listCells.at(getVectorIndex(i,1,j))->m_newVelocity;
        m_listCells.at(getVectorIndex(i,m_noCells-1,j))->m_newVelocity=m_listCells.at(getVectorIndex(i,m_noCells-2,j))->m_newVelocity;
        m_listCells.at(getVectorIndex(i,j,0))->m_newVelocity=m_listCells.at(getVectorIndex(i,j,1))->m_newVelocity;
        m_listCells.at(getVectorIndex(i,j,m_noCells-1))->m_newVelocity=m_listCells.at(getVectorIndex(i,j,m_noCells-2))->m_newVelocity;

      }
      else if (_fieldType=="force")
      {
        m_listCells.at(getVectorIndex(0,i,j))->m_externalForce=m_listCells.at(getVectorIndex(1,i,j))->m_externalForce;
        m_listCells.at(getVectorIndex(m_noCells-1,i,j))->m_externalForce=m_listCells.at(getVectorIndex(m_noCells-2,i,j))->m_externalForce;
        m_listCells.at(getVectorIndex(i,0,j))->m_externalForce=m_listCells.at(getVectorIndex(i,1,j))->m_externalForce;
        m_listCells.at(getVectorIndex(i,m_noCells-1,j))->m_externalForce=m_listCells.at(getVectorIndex(i,m_noCells-2,j))->m_externalForce;
        m_listCells.at(getVectorIndex(i,j,0))->m_externalForce=m_listCells.at(getVectorIndex(i,j,1))->m_externalForce;
        m_listCells.at(getVectorIndex(i,j,m_noCells-1))->m_externalForce=m_listCells.at(getVectorIndex(i,j,m_noCells-2))->m_externalForce;
      }


    }
  }
}

void Grid::setBoundaryValuesFloat(std::string _fieldType)
{
  for (int i=0; i<m_noCells; i++)
  {
    for (int j=0; j<m_noCells; j++)
    {
//      _floatField->at(getVectorIndex(0,i,j))=_floatField->at(getVectorIndex(1,i,j));
//      _floatField->at(getVectorIndex(m_noCells-1,i,j))=_floatField->at(getVectorIndex(m_noCells-2,i,j));
//      _floatField->at(getVectorIndex(i,0,j))=_floatField->at(getVectorIndex(i,1,j));
//      _floatField->at(getVectorIndex(i,m_noCells-1,j))=_floatField->at(getVectorIndex(i,m_noCells-2,j));
//      _floatField->at(getVectorIndex(i,j,0))=_floatField->at(getVectorIndex(i,j,1));
//      _floatField->at(getVectorIndex(i,j,m_noCells-1))=_floatField->at(getVectorIndex(i,j,m_noCells-2));

      if (_fieldType=="pressure")
      {
        m_listCells.at(getVectorIndex(0,i,j))->m_pressure=m_listCells.at(getVectorIndex(1,i,j))->m_pressure;
        m_listCells.at(getVectorIndex(m_noCells-1,i,j))->m_pressure=m_listCells.at(getVectorIndex(m_noCells-2,i,j))->m_pressure;
        m_listCells.at(getVectorIndex(i,0,j))->m_pressure=m_listCells.at(getVectorIndex(i,1,j))->m_pressure;
        m_listCells.at(getVectorIndex(i,m_noCells-1,j))->m_pressure=m_listCells.at(getVectorIndex(i,m_noCells-2,j))->m_pressure;
        m_listCells.at(getVectorIndex(i,j,0))->m_pressure=m_listCells.at(getVectorIndex(i,j,1))->m_pressure;
        m_listCells.at(getVectorIndex(i,j,m_noCells-1))->m_pressure=m_listCells.at(getVectorIndex(i,j,m_noCells-2))->m_pressure;
      }

    }
  }
}

void Grid::setBoundaryVelocity()
{
  for (int i=0; i<m_noCells; i++)
  {
    for (int j=0; j<m_noCells; j++)
    {
//      m_newVelocityField.at(getVectorIndex(0,i,j)).m_x=0;
//      m_newVelocityField.at(getVectorIndex(m_noCells-1,i,j)).m_x=0;
//      m_newVelocityField.at(getVectorIndex(i,0,j)).m_y=0;
//      m_newVelocityField.at(getVectorIndex(i,m_noCells-1,j)).m_y=0;
//      m_newVelocityField.at(getVectorIndex(i,j,0)).m_z=0;
//      m_newVelocityField.at(getVectorIndex(i,j,m_noCells-1)).m_z=0;
      m_listCells.at(getVectorIndex(0,i,j))->m_newVelocity.m_x=0;
      m_listCells.at(getVectorIndex(m_noCells-1,i,j))->m_newVelocity.m_x=0;
      m_listCells.at(getVectorIndex(i,0,j))->m_newVelocity.m_y=0;
      m_listCells.at(getVectorIndex(i,m_noCells-1,j))->m_newVelocity.m_y=0;
      m_listCells.at(getVectorIndex(i,j,0))->m_newVelocity.m_z=0;
      m_listCells.at(getVectorIndex(i,j,m_noCells-1))->m_newVelocity.m_z=0;

    }
  }
}



void Grid::swap()
{
//  m_oldVelocityField=m_newVelocityField;

  for (int i=0; i<pow(m_noCells,3); i++)
  {
    m_listCells.at(i)->m_oldVelocity=m_listCells.at(i)->m_newVelocity;
  }
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
}

void Grid::addForce(float _dt)
{
  for (int i=0; i<pow(m_noCells,3); i++)
  {
//    if (m_solidBoundary.at(i)==false)
//    {
//      m_newVelocityField.at(i)=m_oldVelocityField.at(i)+(m_forceField.at(i)*_dt);
//    }

    if (m_listCells.at(i)->m_isSolid==false)
    {
      m_listCells.at(i)->m_newVelocity=m_listCells.at(i)->m_oldVelocity+(m_listCells.at(i)->m_externalForce*_dt);
    }
  }

  setBoundaryValuesVec3("velocity");
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

        ngl::Vec3 oldVelocity=m_listCells.at(index)->m_oldVelocity;

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
        m_listCells.at(index)->m_tempStoreVec3=advectVelocity;

      }
    }
  }

  //When all advection velocities calculated, set the temporarily stored field to the newVelocity field.
  for (int i=0; i<pow(m_noCells,3); i++)
  {
    m_listCells.at(i)->m_newVelocity=m_listCells.at(i)->m_tempStoreVec3;
  }

  setBoundaryValuesVec3("velocity");
  setBoundaryVelocity();

}

void Grid::diffuse(float _dt)
{

}

void Grid::project(float _dt)
{

//  int sizeOfVectors=(int)(pow((m_noCells-2),3));
//  std::vector<ngl::Vec3> b[sizeOfVectors];
//  std::vector<float> p[sizeOfVectors];
//  std::vector<float> A[1];
  //May be worth checking if boundary cell since don't have to calculate for these

  for (int k=1; k<(m_noCells-1); k++)
  {
    for (int j=1; j<(m_noCells-1); j++)
    {
      for (int i=1; i<(m_noCells-1); i++)
      {

        //Set up b=div(w3)





      }
    }
  }
  //Set up b

  //Set up A

  //Send to linear solver to find p

  //Calulate div(p)

  //Calculate new velocity w4=w3-div(p)


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
  //  fieldVelocity=trilinearInterpVec3(&m_newVelocityField, index0_X, index0_Y, index0_Z, index1_X, index1_Y, index1_Z, indexParticle);
    fieldVelocity=trilinearInterpVec3("velocity", index0_X, index0_Y, index0_Z, index1_X, index1_Y, index1_Z, indexParticle);
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
    fieldVelocity=m_listCells.at(getVectorIndex(index0_X,index0_Y,index0_Z))->m_newVelocity;

  }


  //Return velocity
  return fieldVelocity;
}
