#include <QtGui/QGuiApplication>
#include <iostream>

#include <math.h>

#include "OpenGLWindow.h"


int main(int argc, char **argv)
{
  QGuiApplication app(argc, argv);
  // create an OpenGL format specifier
  QSurfaceFormat format;
  // set the number of samples for multisampling
  // will need to enable glEnable(GL_MULTISAMPLE); once we have a context
  format.setSamples(4);
  #if defined( DARWIN)
    // at present mac osx Mountain Lion only supports GL3.2
    // the new mavericks will have GL 4.x so can change
    format.setMajorVersion(3);
    format.setMinorVersion(2);
  #else
    // with luck we have the latest GL version so set to this
    format.setMajorVersion(4);
    format.setMinorVersion(3);
  #endif
  // now we are going to set to CoreProfile OpenGL so we can't use and old Immediate mode GL
  format.setProfile(QSurfaceFormat::CoreProfile);
  // now set the depth buffer to 24 bits
  format.setDepthBufferSize(24);
  // now we are going to create our scene window
  OpenGLWindow window;
  // and set the OpenGL format
  window.setFormat(format);
  // we can now query the version to see if it worked
  std::cout<<"Profile is "<<format.majorVersion()<<" "<<format.minorVersion()<<"\n";
  // set the window size
  window.resize(1024, 720);
  // and finally show
  window.show();

  return app.exec();

  return EXIT_SUCCESS;
}

//int main()
//{
//  //Test grid index
//  ngl::Vec3 gridpos=ngl::Vec3(-2.5,-2.5,-2.5);
////  m_gridPosition=ngl::Vec3(-5,-5,-5);
//  float gridSize=5.0;
//  int noCells=16;

//  float noise=0.2;
//  float vorticity=4.0;

//  Grid* grid=Grid::createGrid(gridpos,gridSize,noCells,noise, vorticity);

//  float explosionTemp=10000;
//  float addDiv=1.0;

//  ngl::Vec3 explosionPos=ngl::Vec3(0.0,-4.0,0.0);
//  float explosionRadius=10.0;

//  grid->setExplosion(explosionPos, explosionRadius, explosionTemp, addDiv);


//  return EXIT_SUCCESS;
//}
