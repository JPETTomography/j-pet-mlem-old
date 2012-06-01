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



glm::mat4 v_matrix;
glm::mat4 p_matrix;
glm::mat4 mvp_mat;

void ChangeSize(int w, int h) {


  
  glViewport(0, 0, w, h); 
  std::cerr<<"seting projection\n";
  p_matrix=glm::ortho(-250.0f,250.0f,-400.0f,400.0f,-200.0f,200.0f);
  std::cerr<<p_matrix<<std::endl;
  std::cerr<<"set projection\n";
  
}

Shader *shader;

UniformMatrix<4,4> mvp_matrix_uniform("mvp_matrix");




void SetupRC() {


  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glEnable(GL_DEPTH_TEST);

  std::cerr<<"creating shader"<<std::endl;
  shader= new Shader("project.vp","simple.fp"); 
  std::cerr<<"created shader"<<std::endl;
  shader->add_attribute( GLT_ATTRIBUTE_VERTEX, "v_position");
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
  
}


void RenderScene() {  

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
  glUseProgram(shader->shader());
  mvp_mat=p_matrix*v_matrix;
  mvp_matrix_uniform.load_4x4_matrix(shader->shader(),
				     glm::value_ptr(mvp_mat));
  

  glBegin(GL_QUADS);
  glVertex3f(0.0f, -100.0f, -100.0f);
  glVertex3f(0.0f, -100.0f,  100.0f);
  glVertex3f(0.0f,  100.0f,  100.0f);
  glVertex3f(0.0f,  100.0f, -100.0f);
  glEnd();

  glutSwapBuffers();


}


int 
main(int argc, char *argv[]) {

   set_root_logger("tof");
  
  

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_STENCIL);
  glutInitWindowSize(600, 600);
  glutCreateWindow("Cube");
  glutReshapeFunc(ChangeSize);
  glutDisplayFunc(RenderScene);

  GLenum err = glewInit();
  if (GLEW_OK != err) {
        fprintf(stderr, "GLEW Error: %s\n", glewGetErrorString(err));
  }


  SetupRC();
  glutMainLoop();

  std::cout<<"after main loop"<<std::endl;
  return 0;
  
}

