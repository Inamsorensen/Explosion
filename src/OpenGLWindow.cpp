#include <QMouseEvent>
#include <QGuiApplication>

#include <ngl/Camera.h>
#include <ngl/Light.h>
#include <ngl/Material.h>
#include <ngl/NGLInit.h>
#include <ngl/VAOPrimitives.h>
#include <ngl/ShaderLib.h>
#include <ngl/Transformation.h>

#include "OpenGLWindow.h"

//----------------------------------------------------------------------------------------------------------------------
/// @brief the increment for x/y translation with mouse movement
//----------------------------------------------------------------------------------------------------------------------
const static float INCREMENT=0.01;
//----------------------------------------------------------------------------------------------------------------------
/// @brief the increment for the wheel zoom
//----------------------------------------------------------------------------------------------------------------------
const static float ZOOM=0.1;

OpenGLWindow::OpenGLWindow()
{
  //Set rotation of camera to zero for initialisation
  m_rotate=false;
  m_spinXFace=0;
  m_spinYFace=0;

  setTitle("Explosion");

  m_simElapsedTime=0.0;
  m_noFrames=0;

}


OpenGLWindow::~OpenGLWindow()
{
  std::cout<<"Shutting down NGL, removing VAO's and Shaders\n";
}

void OpenGLWindow::resizeGL(QResizeEvent *_event )
{
  m_width=_event->size().width()*devicePixelRatio();
  m_height=_event->size().height()*devicePixelRatio();
  
  // Reset camera with new width and height
  m_camera.setShape(45.0f,(float)width()/height(),0.05f,350.0f);
}

void OpenGLWindow::resizeGL(int _w , int _h)
{
  //Calculate new width and height
  m_camera.setShape(45.0f,(float)_w/_h,0.05f,350.0f);
  m_width=_w*devicePixelRatio();
  m_height=_h*devicePixelRatio();
}

void OpenGLWindow::initializeGL()
{
  ngl::NGLInit::instance();

  // Grey Background
  glClearColor(0.4f, 0.4f, 0.4f, 1.0f);

  // enable depth testing for drawing
  glEnable(GL_DEPTH_TEST);
  // enable multisampling for smoother drawing
  glEnable(GL_MULTISAMPLE);

  //Camera setup
  ngl::Vec3 from(0,0,15);
  ngl::Vec3 to(0,0,0);
  ngl::Vec3 up(0,1,0);
  m_camera.set(from,to,up);
  // set the shape using FOV 45 Aspect Ratio based on Width and Height
  // The final two are near and far clipping planes of 0.5 and 10
  m_camera.setShape(60,(float)720.0/576.0,0.5,150);

  //Shader setup
  //Setup Phong shader
  ngl::ShaderLib *shader=ngl::ShaderLib::instance();
  shader->createShaderProgram("Phong");
  shader->attachShader("PhongVertex",ngl::ShaderType::VERTEX);
  shader->attachShader("PhongFragment",ngl::ShaderType::FRAGMENT);
  shader->loadShaderSource("PhongVertex","shaders/Phong.vs");
  shader->loadShaderSource("PhongFragment","shaders/Phong.fs");
  shader->compileShader("PhongVertex");
  shader->compileShader("PhongFragment");
  shader->attachShaderToProgram("Phong","PhongVertex");
  shader->attachShaderToProgram("Phong","PhongFragment");
  shader->bindAttribute("Phong",0,"inVert");
  shader->bindAttribute("Phong",1,"inUV");
  shader->bindAttribute("Phong",2,"inNormal");
  shader->linkProgramObject("Phong");
  (*shader)["Phong"]->use();
  shader->setShaderParam1i("Normalize",1);

  // now pass the modelView and projection values to the shader
  // the shader will use the currently active material and light0 so set them
  ngl::Material m(ngl::STDMAT::GOLD);
  m.loadToShader("material");
  ngl::Light light(ngl::Vec3(2,2,20),ngl::Colour(1,1,1,1),ngl::Colour(1,1,1,1),ngl::LightModes::POINTLIGHT);

  //Create light
  ngl::Mat4 iv=m_camera.getViewMatrix();
  iv.transpose();
  light.setTransform(iv);
  light.setAttenuation(1,0,0);
  light.enable();
  light.loadToShader("light");
  ngl::VAOPrimitives *prim=ngl::VAOPrimitives::instance();
  prim->createSphere("sphere", 0.1,10);

  //Setup colour shader
  shader->createShaderProgram("Colour");
  shader->attachShader("ColourVertex", ngl::ShaderType::VERTEX);
  shader->attachShader("ColourFragment", ngl::ShaderType::FRAGMENT);
  shader->loadShaderSource("ColourVertex","shaders/colour.vs");
  shader->loadShaderSource("ColourFragment","shaders/colour.fs");
  shader->compileShader("ColourVertex");
  shader->compileShader("ColourFragment");
  shader->attachShaderToProgram("Colour","ColourVertex");
  shader->attachShaderToProgram("Colour","ColourFragment");
  shader->bindAttribute("Colour",0,"inVertex");
  shader->bindAttribute("Colour",1,"inUV");
  shader->bindAttribute("Colour",2,"inNormal");
  shader->linkProgramObject("Colour");
  (*shader)["Colour"]->use();
  shader->setShaderParam1i("Normalize",1);



  //Setup simulation
  simulationSetup();

  //Setup VAO for bounding box
  buildVAO();

  //Set camera
  m_emitter->setCamera(&m_camera);
  m_emitter->setShaderName("Phong");



  // as re-size is not explicitly called we need to do this.
  glViewport(0,0,width(),height());

  //Start timer.
  startTimer(10);


}

void OpenGLWindow::simulationSetup()
{
  //Time step setup
  m_simTimeStep=0.01;


  //Grid setup
  m_gridPosition=ngl::Vec3(-2.5,-2.5,-2.5);
//  m_gridPosition=ngl::Vec3(-5,-5,-5);
  m_gridSize=5.0;
  m_noCells=16;

  float noise=2.0;
  float vorticityConstant=30.0;

  Grid* grid=Grid::createGrid(m_gridPosition,m_gridSize,m_noCells, noise, vorticityConstant);

  std::vector<ngl::Vec3> zeroVectorField;
  std::vector<float> zeroFloatField;
  std::vector<ngl::Vec3> gravity;
  for (int i=0; i<pow(m_noCells,3); i++)
  {
    zeroVectorField.push_back(ngl::Vec3(0.0,0.0,0.0));
    zeroFloatField.push_back(0.0);
    gravity.push_back(ngl::Vec3(0.0,-9.81,0.0));
  }

  //grid->setVelocityField(zeroVectorField);
//  grid->setPressureField(zeroFloatField);
  //grid->setForceField(gravity);


  //Particle setup
  float particleMass=0.1;
  float particleRadius=0.1;
  ngl::Vec3 initialVelocity=ngl::Vec3(0.0,0.0,0.0);
  float particleLifeTime=10;
  float particleDrag=1;
  float particleInitialTemperature=293.0;


  //Emitter setup
  m_emitterPosition=ngl::Vec3(0.0,0.0,0.0);
  m_emitterRadius=0.3;
  m_noParticles=2000;
//  m_noParticles=1;
  m_emissionRate=10;

  m_emitter=new Emitter(m_emitterPosition,m_emitterRadius,m_noParticles,m_emissionRate);
  m_emitter->setUpParticles(particleMass, initialVelocity, particleLifeTime, particleRadius, particleDrag, particleInitialTemperature);

}

void OpenGLWindow::buildVAO()
{
  m_vao.reset( ngl::VertexArrayObject::createVOA(GL_LINES));
  m_vao->bind();

//Line indices
  const static GLubyte indices[]=  {
                                      0,1,1,2,2,3,3,0, //top
                                      0,4,4,5,5,1,1,0, //back
                                      0,4,4,7,7,3,3,0, //left
                                      3,2,2,6,6,7,7,3, //front
                                      7,6,6,5,5,4,4,7, //bottom
                                   };
//Vertices of lines
   GLfloat vertices[] = {0,1,0,
                         1,1,0,
                         1,1,1,
                         0,1,1,
                         0,0,0,
                         1,0,0,
                         1,0,1,
                         0,0,1
                        };
//Colour of each vertex
   GLfloat colours[] = {1,0,0,
                        1,0,0,
                        1,0,0,
                        1,0,0,
                        1,0,0,
                        1,0,0,
                        1,0,0,
                        1,0,0
                       };


//Feed data into the VAO
   m_vao->setIndexedData(24*sizeof(GLfloat),vertices[0],sizeof(indices),&indices[0],GL_UNSIGNED_BYTE,GL_STATIC_DRAW);
   m_vao->setVertexAttributePointer(0,3,GL_FLOAT,0,0);
   m_vao->setIndexedData(24*sizeof(GLfloat),colours[0],sizeof(indices),&indices[0],GL_UNSIGNED_BYTE,GL_STATIC_DRAW);
   m_vao->setVertexAttributePointer(1,3,GL_FLOAT,0,0);
   m_vao->setNumIndices(sizeof(indices));
   m_vao->unbind();

}



void OpenGLWindow::paintGL()
{
  // clear the screen and depth buffer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0,0,m_width,m_height);

  //Set ModelMatrix for camera orientation
  ngl::Mat4 ModelMatrix_Camera;
  //Create rotation matrices
  ngl::Mat4 rotationX;
  ngl::Mat4 rotationY;
  rotationX.rotateX(m_spinXFace);
  rotationY.rotateY(m_spinYFace);

  //Multiply the rotation matrices and add translation
  ModelMatrix_Camera=rotationX*rotationY;
  ModelMatrix_Camera.m_30=m_modelPosition.m_x;
  ModelMatrix_Camera.m_31=m_modelPosition.m_y;
  ModelMatrix_Camera.m_32=m_modelPosition.m_z;


  //Draw bounding box
  ngl::ShaderLib *shader=ngl::ShaderLib::instance();
  shader->use("Colour");

  ngl::Mat4 MVP;
  ngl::Transformation ModelMatrix_BBox;
  ModelMatrix_BBox.setPosition(m_gridPosition);
  ModelMatrix_BBox.setScale(m_gridSize,m_gridSize,m_gridSize);


  //MVP calculated from multiplying ModelMatrix_BBox with ModelMatrix_Camera?
  ngl::Mat4 ModelMatrix;
  ModelMatrix=ModelMatrix_BBox.getMatrix()*ModelMatrix_Camera;

  MVP=ModelMatrix*m_camera.getVPMatrix();
//  MVP=ModelMatrix_Camera*m_camera.getVPMatrix();
  shader->setShaderParamFromMat4("MVP",MVP);

  m_vao->bind();
  m_vao->draw();
  m_vao->unbind();


  //Draw particles
  m_emitter->renderParticles(ModelMatrix_Camera);

  //Increase frame count
  m_noFrames+=1;
 }


//----------------------------------------------------------------------------------------------------------------------
void OpenGLWindow::mouseMoveEvent (QMouseEvent * _event)
{
  // note the method buttons() is the button state when event was called
  // this is different from button() which is used to check which button was
  // pressed when the mousePress/Release event is generated
  if(m_rotate && _event->buttons() == Qt::LeftButton)
  {
    int diffx=_event->x()-m_origX;
    int diffy=_event->y()-m_origY;
    m_spinXFace += (float) 0.5f * diffy;
    m_spinYFace += (float) 0.5f * diffx;
    m_origX = _event->x();
    m_origY = _event->y();
    update();

  }
        // right mouse translate code
  else if(m_translate && _event->buttons() == Qt::RightButton)
  {
    int diffX = (int)(_event->x() - m_origXPos);
    int diffY = (int)(_event->y() - m_origYPos);
    m_origXPos=_event->x();
    m_origYPos=_event->y();
    m_modelPosition.m_x += INCREMENT * diffX;
    m_modelPosition.m_y -= INCREMENT * diffY;
    update();

   }
}



//----------------------------------------------------------------------------------------------------------------------
void OpenGLWindow::mousePressEvent ( QMouseEvent * _event)
{
  // this method is called when the mouse button is pressed in this case we
  // store the value where the maouse was clicked (x,y) and set the Rotate flag to true
  if(_event->button() == Qt::LeftButton)
  {
    m_origX = _event->x();
    m_origY = _event->y();
    m_rotate =true;
  }
  // right mouse translate mode
  else if(_event->button() == Qt::RightButton)
  {
    m_origXPos = _event->x();
    m_origYPos = _event->y();
    m_translate=true;
  }

}


//----------------------------------------------------------------------------------------------------------------------
void OpenGLWindow::mouseReleaseEvent ( QMouseEvent * _event )
{
  // this event is called when the mouse button is released
  // we then set Rotate to false
  if (_event->button() == Qt::LeftButton)
  {
    m_rotate=false;
  }
        // right mouse translate mode
  if (_event->button() == Qt::RightButton)
  {
    m_translate=false;
  }

}


//----------------------------------------------------------------------------------------------------------------------
void OpenGLWindow::wheelEvent(QWheelEvent *_event)
{

  // check the diff of the wheel position (0 means no change)
  if(_event->delta() > 0)
  {
    m_modelPosition.m_z+=ZOOM;
  }
  else if(_event->delta() <0 )
  {
    m_modelPosition.m_z-=ZOOM;
  }
  update();

}

//----------------------------------------------------------------------------------------------------------------------

void OpenGLWindow::keyPressEvent(QKeyEvent *_event)
{
  // this method is called every time the main window recives a key event.
  // we then switch on the key value and set the camera in the GLWindow
  switch (_event->key())
  {
  // escape key to quite
  case Qt::Key_Escape : QGuiApplication::exit(EXIT_SUCCESS); break;
  // turn on wirframe rendering
  case Qt::Key_W : glPolygonMode(GL_FRONT_AND_BACK,GL_LINE); break;
  // turn off wire frame
  case Qt::Key_S : glPolygonMode(GL_FRONT_AND_BACK,GL_FILL); break;
  // show full screen
  case Qt::Key_F : showFullScreen(); break;
  // show windowed
  case Qt::Key_N : showNormal(); break;
   // update particle:
  case Qt::Key_U :
  {
    //Update grid
    Grid* grid=Grid::getGrid();
    grid->update(m_simTimeStep);

    //Update particles
    m_emitter->update(m_simTimeStep);

    //Increase elapsed time
    m_simElapsedTime+=m_simTimeStep;

    //Render new frame
    update();

    break;
  }
  default : break;
  }
  update();
}

void OpenGLWindow::timerEvent(QTimerEvent *_event )
{
  //Update grid
  Grid* grid=Grid::getGrid();
  grid->update(m_simTimeStep);

  //Update particles
  m_emitter->update(m_simTimeStep);

  //Increase elapsed time
  m_simElapsedTime+=m_simTimeStep;

  //Render new frame
  update();
}

