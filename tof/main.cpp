#include <ios>
#include <iostream>

#include <stdio.h>

#ifndef __APPLE__
  #include <GL/gl.h>
  #include <GL/glut.h>
#else
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "glmutils.h"
#include "glutils.h"
#include "Shader.h"
#include "Uniform.h"

#include "density_plot.h"
#include "geometry_plot.h"

glm::mat4 v_matrix;
glm::mat4 p_matrix;
glm::mat4 mvp_mat;
#include "pixel_grid.h"
#include "topet_simulator.h"
#include "detector_view.h"
#include "phantom_view.h"
#include "event_view.h"
bool event_mode = false;

void ChangeSize(int w, int h) {
  glViewport(0, 0, w, h);
  GLfloat aspect = (GLfloat)w/(GLfloat)h;
  std::cerr << "seting projection\n";
  GLfloat height = 800.0f;
  p_matrix = glm::ortho(-height*aspect/2.0f, height*aspect/2.0f, -height/2.0f, height/2.0f, 99.0f, 101.0f);
  std::cerr << p_matrix << std::endl;
  std::cerr << "set projection\n";
}

Shader *shader;

UniformMatrix<4, 4> mvp_matrix_uniform("mvp_matrix");
DensityPlot *density_plot;
GeometryPlot *geometry_plot;
DetectorView<GLfloat> *detector_view;
PhantomView *phantom_view;
EventView<GLfloat> *event_view;

#include "tof_event.h"
#include "tof_detector.h"
#include "phantom.h"
#include "tof_monte_carlo.h"

TOPETSimulator<GLfloat>  simulator;

void SetupRC() {
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  //glEnable(GL_DEPTH_TEST);

  std::cerr << "creating shader" << std::endl;
  shader = new Shader("project.vp", "coldtohot.fp");
  std::cerr << "created shader" << std::endl;
  shader->add_attribute( GLT_ATTRIBUTE_VERTEX, "v_position");
  shader->add_attribute( GLT_ATTRIBUTE_TEXTURE0, "v_texture0");
  shader->link();
  mvp_matrix_uniform.add_shader(
        shader->shader()
        );

  glUseProgram(shader->shader());
  glm::vec3 eye(-100.0, 0.0, 0.0);
  glm::vec3 center(0.0, 0.0, 0.0);
  glm::vec3 up(0.0, 1.0, 0.0);

  std::cerr << "setting view matrix" << std::endl;
  v_matrix = glm::lookAt(eye, center, up);
  std::cerr << v_matrix << std::endl;
  std::cerr << "set view matrix" << std::endl;
  GLuint width = 100;
  GLuint height = 100;
  GLuint texture_size = width*height;
  GLfloat *pBits = (GLfloat *)malloc(sizeof(GLfloat)*texture_size);
  for(int i = 0;i<texture_size;i++) {
    pBits[i] = i*1.0/texture_size;
  }
  density_plot = new DensityPlot(100, 100);
  density_plot->set_vertex(0, 0, -250, -250, 0, 0);
  density_plot->set_vertex(1, 0, -250, 250, 1, 0);
  density_plot->set_vertex(2, 0, 250, 250, 1, 1);
  density_plot->set_vertex(3, 0, 250, -250, 0, 1);
  simulator.init();
  //simulator.simulate_from_single_point(100000);
  simulator.simulate_from_phantom(10000000);
  density_plot->set_pixmap(simulator.emitted_density());

  geometry_plot = new GeometryPlot;

  detector_view = new DetectorView<GLfloat>(geometry_plot, simulator.detector());
  phantom_view = new PhantomView(geometry_plot, simulator.phantom());
  event_view = new EventView<GLfloat>(geometry_plot, simulator.detector(), simulator.tof_begin(), simulator.tof_end(), simulator.detected_begin()
              );
}
void keyboardHandler(unsigned char pressed, int x, int y) {
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
  case 'd':
    density_plot->set_pixmap(simulator.detected_density());
    glutPostRedisplay();
    break;
  case 'v':
    event_mode = !event_mode;
    glutPostRedisplay();
    break;
  case 'n':
    if(event_mode) {
      ++(*event_view);
      glutPostRedisplay();
    }
    break;
  case '2':
    if(event_mode) {
      event_view->togle_two_sets_mode();
      glutPostRedisplay();
    }
    break;
  case 'x':
  case 'q':
    exit(0);
    break;
  }
}

void RenderScene() {

  std::cerr << "renderimg scene" << std::endl;
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
  glUseProgram(shader->shader());
  mvp_mat = p_matrix*v_matrix;
  mvp_matrix_uniform.load_4x4_matrix(shader->shader(), glm::value_ptr(mvp_mat));
  density_plot->set_mvp_matrix(mvp_mat);

  density_plot->render();

  geometry_plot->set_mvp_matrix(mvp_mat);
  geometry_plot->renderStart();

  phantom_view->render();
  detector_view->render();
  if(event_mode)
    event_view->render();

  geometry_plot->renderEnd();
  glutSwapBuffers();
}
int
main(int argc, char *argv[]) {

  set_root_logger("tof");
  glutInit(&argc, argv);
#ifndef __APPLE__
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_STENCIL);
#else
  glutInitDisplayString("double rgba depth>=16 stencil profile=32");
#endif
  glutInitWindowSize(600, 600);
  glutCreateWindow("ToF");
  glutReshapeFunc(ChangeSize);
  glutDisplayFunc(RenderScene);
  glutKeyboardFunc(keyboardHandler);
  SetupRC();
  glutMainLoop();

  std::cout << "after main loop" << std::endl;
  return 0;
}

