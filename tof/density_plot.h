#ifndef __DENSITY_PLOT_H__
#define __DENSITY_PLOT_H__


#include<GL/glew.h>
#include<GL/gl.h>

#include"Uniform.h"


class DensityPlot  {
public:
  DensityPlot(GLuint  width, GLuint height):
    width_(width),height_(height),
    mvp_matrix_uniform_("mvp_matrix"),
    texture_uniform("texture0"),
    start_color_uniform("start_color"),
    end_color_uniform("end_color"),
    vertices_(4*5,0.0f) {

    std::cerr<<"creating density plot"<<std::endl;
  std::cerr<<"creating shader"<<std::endl;
  shader_= new Shader("project.vp","coldtohot.fp"); 
  std::cerr<<"created shader"<<std::endl;
  shader_->add_attribute( GLT_ATTRIBUTE_VERTEX, "v_position");
  shader_->add_attribute( GLT_ATTRIBUTE_TEXTURE0, "v_texture0");
  shader_->link();
  texture_uniform.add_shader(shader_->shader()); 
  mvp_matrix_uniform_.add_shader(
				shader_->shader()
				);
  glGenTextures(1, &map_texture_id_);
  glBindTexture(GL_TEXTURE_2D, map_texture_id_);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    
  //  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

  std::cerr<<"created density plot"<<std::endl;
  };

  void set_mvp_matrix(glm::mat4 mvp) {
    mvp_=mvp;
  }

  void set_pixmap(const GLfloat *p_bits) {
    glBindTexture(GL_TEXTURE_2D, map_texture_id_);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    std::cerr<<"p bits "<<p_bits[0]<<std::endl;
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, width_, height_, 0,
		 GL_RED, GL_FLOAT, p_bits);
  }

  

  void set_vertex(GLuint i,GLfloat x,GLfloat y,GLfloat z,GLfloat s,GLfloat t) {
    vertices_[5*i]=x;
    vertices_[5*i+1]=y;
    vertices_[5*i+2]=z;
    vertices_[5*i+3]=s;
    vertices_[5*i+4]=t;    
  }

  void render() {
    GLint current_shader;
    glGetIntegerv(GL_CURRENT_PROGRAM,&current_shader);

    glUseProgram(shader_->shader());
    texture_uniform.load_concrete(shader_->shader(),0);
    mvp_matrix_uniform_.load_4x4_matrix(shader_->shader(),
				     glm::value_ptr(mvp_));
    glBegin(GL_QUADS);
    for(int i=0;i<4;i++) {
      glVertexAttrib2f(GLT_ATTRIBUTE_TEXTURE0,vertices_[5*i+3],vertices_[5*i+4]);
      glVertex3f(vertices_[5*i],vertices_[5*i+1],vertices_[5*i+2] );

    }
    glEnd();
    glUseProgram(current_shader);
  }

 private:
  GLuint width_;
  GLuint height_;
  GLuint map_texture_id_;
  Shader *shader_;
  UniformMatrix<4,4> mvp_matrix_uniform_;
  UniformScalar<GLuint> texture_uniform;
  UniformVector<GLfloat,4> start_color_uniform;
  UniformVector<GLfloat,4> end_color_uniform;
  glm::mat4 mvp_;

  std::vector<GLfloat> vertices_;
  
};

#endif
