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



GLuint textureID;


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
UniformScalar<GLuint> texture_uniform("texture0");



void SetupRC() {


  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glEnable(GL_DEPTH_TEST);

  std::cerr<<"creating shader"<<std::endl;
  shader= new Shader("project.vp","simple.fp"); 
  std::cerr<<"created shader"<<std::endl;
  shader->add_attribute( GLT_ATTRIBUTE_VERTEX, "v_position");
  shader->add_attribute( GLT_ATTRIBUTE_TEXTURE0, "v_texture0");
  shader->link();

  
  mvp_matrix_uniform.add_shader(
				shader->shader()
				);
  texture_uniform.add_shader(shader->shader());

  glUseProgram(shader->shader());
  glm::vec3 eye(-100.0,0.0,0.0);
  glm::vec3 center(0.0,0.0,0.0);
  glm::vec3 up(0.0,1.0,0.0);

  std::cerr<<"setting view matrix"<<std::endl;
  v_matrix=glm::lookAt(eye,center,up);
  std::cerr<<v_matrix<<std::endl;
  std::cerr<<"set view matrix"<<std::endl;
  
  
  GLuint width=100;
  GLuint height=100;
  GLuint texture_size=width*height*3;
  GLbyte *pBits=(GLbyte *)malloc(sizeof(unsigned char)*texture_size);
  for(int i=0;i<texture_size;i+=3) {
    pBits[i]=128;
    pBits[i+1]=0;
    pBits[i+2]=32;
  }
  
  glGenTextures(1, &textureID);
  glBindTexture(GL_TEXTURE_2D, textureID);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0,
  GL_RGB, GL_UNSIGNED_BYTE, pBits);

}


void RenderScene() {  

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
  glUseProgram(shader->shader());
  mvp_mat=p_matrix*v_matrix;
  mvp_matrix_uniform.load_4x4_matrix(shader->shader(),
				     glm::value_ptr(mvp_mat));
  
  texture_uniform.load_concrete(shader->shader(),0);

  glBegin(GL_QUADS);
  glVertexAttrib2f(GLT_ATTRIBUTE_TEXTURE0,0.0f,0.0f);
  glVertex3f(0.0f, -100.0f, -100.0f);
  glVertexAttrib2f(GLT_ATTRIBUTE_TEXTURE0,0.0f,1.0f);
  glVertex3f(0.0f, -100.0f,  100.0f);
  glVertexAttrib2f(GLT_ATTRIBUTE_TEXTURE0,1.0f,1.0f);
  glVertex3f(0.0f,  100.0f,  100.0f);
  glVertexAttrib2f(GLT_ATTRIBUTE_TEXTURE0,1.0f,0.0f);
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

