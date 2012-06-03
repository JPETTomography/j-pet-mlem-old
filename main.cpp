#include<ios>
#include<iostream>

#include<stdio.h>

#include<GL/glew.h>
#include<GL/gl.h>
#include<GL/glut.h>


#include<glm/glm.hpp>
#include<glm/gtc/type_ptr.hpp>
#include<glm/gtc/matrix_transform.hpp>

#include"glmutils.h"
#include"glutils.h"
#include"Shader.h"
#include"Uniform.h"

#include"density_plot.h"
#include"geometry_plot.h"

glm::mat4 v_matrix;
glm::mat4 p_matrix;
glm::mat4 mvp_mat;


#include"pixel_grid.h"
#include"topet_simulator.h"


void ChangeSize(int w, int h) {


  
  glViewport(0, 0, w, h); 
  GLfloat aspect=(GLfloat)w/(GLfloat)h;
  std::cerr<<"seting projection\n";
  GLfloat height= 800.0f;
  p_matrix=glm::ortho(-height*aspect/2.0f,height*aspect/2.0f,
		      -height/2.0f,height/2.0f,99.0f,101.0f);
  std::cerr<<p_matrix<<std::endl;
  std::cerr<<"set projection\n";
  
}

Shader *shader;

UniformMatrix<4,4> mvp_matrix_uniform("mvp_matrix");


DensityPlot *density_plot;
GeometryPlot *geometry_plot;


#include"tof_event.h"
#include"tof_detector.h"
#include"phantom.h"
#include"tof_monte_carlo.h"

  TOPETSimulator<GLfloat>  simulator;

void SetupRC() {


  
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  //glEnable(GL_DEPTH_TEST);

  std::cerr<<"creating shader"<<std::endl;
  shader= new Shader("project.vp","coldtohot.fp"); 
  std::cerr<<"created shader"<<std::endl;
  shader->add_attribute( GLT_ATTRIBUTE_VERTEX, "v_position");
  shader->add_attribute( GLT_ATTRIBUTE_TEXTURE0, "v_texture0");
  shader->link();

  
  mvp_matrix_uniform.add_shader(
				shader->shader()
				);

  glUseProgram(shader->shader());
  glm::vec3 eye(-100.0,0.0,0.0);
  glm::vec3 center(0.0,0.0,0.0);
  glm::vec3 up(0.0,1.0,0.0);

  std::cerr<<"setting view matrix"<<std::endl;
  v_matrix=glm::lookAt(eye,center,up);
  std::cerr<<v_matrix<<std::endl;
  std::cerr<<"set view matrix"<<std::endl;
  
  
  GLuint width = 100;
  GLuint height= 100;
  GLuint texture_size=width*height;
  GLfloat *pBits=(GLfloat *)malloc(sizeof(GLfloat)*texture_size);
  for(int i=0;i<texture_size;i++) {
    pBits[i]=i*1.0/texture_size;
  }


  density_plot= new DensityPlot(100,100);
  density_plot->set_vertex(0,0,-250,-250,0,0);
  density_plot->set_vertex(1,0,-250,250,1,0);
  density_plot->set_vertex(2,0, 250,250,1,1);
  density_plot->set_vertex(3,0, 250,-250,0,1);



  simulator.init();
  //simulator.simulate_from_single_point(100000);
  simulator.simulate_from_phantom(10000000);
  density_plot->set_pixmap(simulator.emitted_density());
  
  geometry_plot=new GeometryPlot;
}


void keyboardHandler(unsigned char pressed,int x, int y) {
  switch (pressed) {
  case 'a':  
    density_plot->set_pixmap(simulator.acceptance_density());
    glutPostRedisplay();
    break;
  case 't':  
    density_plot->set_pixmap(simulator.tof_density());
    glutPostRedisplay();
    break;
  case 'e':  
    density_plot->set_pixmap(simulator.emitted_density());    
    glutPostRedisplay();
    break;
  case 'x':  
    exit(0);
    break;
  }

  
}

void RenderScene() {  

  std::cerr<<"renderimg scene"<<std::endl;
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
  glUseProgram(shader->shader());
  mvp_mat=p_matrix*v_matrix;
  mvp_matrix_uniform.load_4x4_matrix(shader->shader(),
				     glm::value_ptr(mvp_mat));
  density_plot->set_mvp_matrix(mvp_mat);
  
  density_plot->render();

  geometry_plot->set_mvp_matrix(mvp_mat);
  geometry_plot->renderStart();
  geometry_plot->renderZYEllipse(0,0,100,200,0.0);
  geometry_plot->renderZYEllipse(0,-100,50,70,M_PI/3.0);
  geometry_plot->renderZYEllipse(20,150,10,17,M_PI/4.0);
  geometry_plot->renderZYRectangle(-250,-360,250,-340);
  geometry_plot->renderZYRectangle(-250,340,250,360);


  geometry_plot->renderEnd();

  glutSwapBuffers();


}

 
int 
main(int argc, char *argv[]) {

  set_root_logger("tof");
  

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_STENCIL);
  glutInitWindowSize(600, 600);
  glutCreateWindow("ToF");
  glutReshapeFunc(ChangeSize);
  glutDisplayFunc(RenderScene);
  glutKeyboardFunc(keyboardHandler);

  GLenum err = glewInit();
  if (GLEW_OK != err) {
        fprintf(stderr, "GLEW Error: %s\n", glewGetErrorString(err));
  }


  SetupRC();
  glutMainLoop();

  std::cout<<"after main loop"<<std::endl;
  return 0;
  
}

