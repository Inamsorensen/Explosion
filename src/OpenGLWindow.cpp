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
  ngl::Vec3 from(0,0,1);
  ngl::Vec3 to(0,0,0);
  ngl::Vec3 up(0,1,0);
  m_camera.set(from,to,up);
  // set the shape using FOV 45 Aspect Ratio based on Width and Height
  // The final two are near and far clipping planes of 0.5 and 10
  m_camera.setShape(60,(float)720.0/576.0,0.5,150);

  //Shader setup
  //Setup Phong Gold shader - Used for burning particles
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

  //Setup Phong Silver shader - Used for soot
  shader->createShaderProgram("PhongSilver");
  shader->attachShader("PhongVertex",ngl::ShaderType::VERTEX);
  shader->attachShader("PhongFragment",ngl::ShaderType::FRAGMENT);
  shader->loadShaderSource("PhongVertex","shaders/Phong.vs");
  shader->loadShaderSource("PhongFragment","shaders/Phong.fs");
  shader->compileShader("PhongVertex");
  shader->compileShader("PhongFragment");
  shader->attachShaderToProgram("PhongSilver","PhongVertex");
  shader->attachShaderToProgram("PhongSilver","PhongFragment");
  shader->bindAttribute("PhongSilver",0,"inVert");
  shader->bindAttribute("PhongSilver",1,"inUV");
  shader->bindAttribute("PhongSilver",2,"inNormal");
  shader->linkProgramObject("PhongSilver");
  (*shader)["PhongSilver"]->use();
  shader->setShaderParam1i("Normalize",1);

  // the shader will use the currently active material and light0 so set them
  ngl::Material mSilver(ngl::STDMAT::SILVER);
  mSilver.loadToShader("material");
  //Attach light
  light.loadToShader("light");

  //Setup Phong Copper shader - Used for unignited particles
  shader->createShaderProgram("PhongCopper");
  shader->attachShader("PhongVertex",ngl::ShaderType::VERTEX);
  shader->attachShader("PhongFragment",ngl::ShaderType::FRAGMENT);
  shader->loadShaderSource("PhongVertex","shaders/Phong.vs");
  shader->loadShaderSource("PhongFragment","shaders/Phong.fs");
  shader->compileShader("PhongVertex");
  shader->compileShader("PhongFragment");
  shader->attachShaderToProgram("PhongCopper","PhongVertex");
  shader->attachShaderToProgram("PhongCopper","PhongFragment");
  shader->bindAttribute("PhongCopper",0,"inVert");
  shader->bindAttribute("PhongCopper",1,"inUV");
  shader->bindAttribute("PhongCopper",2,"inNormal");
  shader->linkProgramObject("PhongCopper");
  (*shader)["PhongCopper"]->use();
  shader->setShaderParam1i("Normalize",1);

  // the shader will use the currently active material and light0 so set them
  ngl::Material mCopper(ngl::STDMAT::COPPER);
  mCopper.loadToShader("material");
  //Attach light
  light.loadToShader("light");

  //Setup colour shader - Used to draw bounding box
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
  m_explosionController=ExplosionController::instance();

  //Setup VAO for bounding box
  buildVAO();

  //Send render variables to explosion controller
  m_explosionController->setRenderVariables(&m_camera, "Phong");


  // as re-size is not explicitly called we need to do this.
  glViewport(0,0,width(),height());

  //Start timer.
  startTimer(10);


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
  ModelMatrix_BBox.setPosition(m_explosionController->getGridPosition());
  float gridSize=m_explosionController->getGridSize();
  ModelMatrix_BBox.setScale(gridSize,gridSize,gridSize);


  //MVP calculated from multiplying ModelMatrix_BBox with ModelMatrix_Camera
  ngl::Mat4 ModelMatrix;
  ModelMatrix=ModelMatrix_BBox.getMatrix()*ModelMatrix_Camera;

  MVP=ModelMatrix*m_camera.getVPMatrix();
  shader->setShaderParamFromMat4("MVP",MVP);

  m_vao->bind();
  m_vao->draw();
  m_vao->unbind();


  //Draw particles
  m_explosionController->render(ModelMatrix_Camera);

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
  case Qt::Key_F : showFullScreen(); break;
  case Qt::Key_E : m_explosionController->toggleAlembicExport(); break;
  default : break;
  }
//  update();
}

void OpenGLWindow::timerEvent(QTimerEvent *_event )
{
  //Update explosion calculation
  m_explosionController->update();

  //Render new frame
  update();
}

