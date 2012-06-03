#ifndef __GEOMETRY_PLOT_H__
#define __GEOMETRY_PLOT_H__


class GeometryPlot  {
public:
  GeometryPlot():
    mvp_matrix_uniform_("mvp_matrix"),
    start_color_uniform("start_color"),
    end_color_uniform("end_color") {

    std::cerr<<"creating density plot"<<std::endl;
    std::cerr<<"creating shader"<<std::endl;
    shader_= new Shader("project.vp","simple.fp"); 
    std::cerr<<"created shader"<<std::endl;

    shader_->add_attribute( GLT_ATTRIBUTE_VERTEX, "v_position");
    shader_->link();
    mvp_matrix_uniform_.add_shader(
				   shader_->shader()
				   );
    std::cerr<<"created geometry plot"<<std::endl;
  };


  void set_mv_matrix(glm::mat4 mv) {
    mv_=mv;
  }

  void set_mvp_matrix(glm::mat4 mvp) {
    mvp_=mvp;
  }
  

  void renderStart() {
    std::cerr<<"renderStart"<<std::endl;

    std::cerr<<"get current program"<<std::endl;
    glGetIntegerv(GL_CURRENT_PROGRAM,&current_shader_);

    std::cerr<<"use shader"<<std::endl;
    glUseProgram(shader_->shader());
    
    std::cerr<<"load mvp uniform"<<std::endl;
    mvp_matrix_uniform_.load_4x4_matrix(shader_->shader(),
					glm::value_ptr(mvp_));
  }

  void drawCircle(GLfloat r) {
  const int n = 120;
    glBegin(GL_LINE_STRIP) ;
    for(int i=0;i<n;i++ ){      
      GLfloat z=r*cos((2.0*M_PI/n)*i);
      GLfloat y=r*sin((2.0*M_PI/n)*i);
      glVertex3f(0.0,y,z);
    }
    glEnd();
  }

  void drawDisk(GLfloat r) {
  const int n = 120;
    glBegin(GL_POLYGON) ;
    for(int i=0;i<n;i++ ){      
      GLfloat z=r*cos((2.0*M_PI/n)*i);
      GLfloat y=r*sin((2.0*M_PI/n)*i);
      glVertex3f(0.0,y,z);
    }
    glEnd();
  }

  void renderZYCircle(GLfloat zc, GLfloat yc, GLfloat r,bool filled = false) {
    std::cerr<<"rendering circle"<<std::endl;

    glm::mat4  M=mvp_;
    std::cerr<<M<<std::endl;
    M=glm::translate(M,glm::vec3(0,yc,zc));
    M=glm::scale(M,glm::vec3(0,r,r));

    
    std::cerr<<M<<std::endl;
     mvp_matrix_uniform_.load_4x4_matrix(shader_->shader(),
					glm::value_ptr(M));
     drawCircle(1);
     if(filled) 
       drawDisk(1);
    mvp_matrix_uniform_.load_4x4_matrix(shader_->shader(),
					glm::value_ptr(mvp_));
  }

  void renderZYLine(GLfloat z1, GLfloat y1,GLfloat z2, GLfloat y2) {
    glBegin(GL_LINE_STRIP);
    glVertex3f(0.0f,y1,z1);
    glVertex3f(0.0f,y2,z2);
    glEnd();
  }

  void renderZYPoint(GLfloat z1, GLfloat y1) {
    glBegin(GL_POINTS);
    glVertex3f(0.0f,y1,z1);
    glEnd();
  }

  void renderZYEllipse(GLfloat zc, GLfloat yc, GLfloat rz, GLfloat ry, 
		       GLfloat theta = 0,bool filled=false) {
    std::cerr<<"rendering circle"<<std::endl;
    GLfloat degree_theta=180.0*theta/M_PI;
    glm::mat4  M=mvp_;
    std::cerr<<M<<std::endl;
    M=glm::translate(M,glm::vec3(0,yc,zc));
    M=glm::rotate(M,degree_theta,glm::vec3(-1,0,0));
    M=glm::scale(M,glm::vec3(0,ry,rz));

    std::cerr<<M<<std::endl;
     mvp_matrix_uniform_.load_4x4_matrix(shader_->shader(),
					glm::value_ptr(M));
     drawCircle(1);
     if(filled)
       drawDisk(1);

    mvp_matrix_uniform_.load_4x4_matrix(shader_->shader(),
					glm::value_ptr(mvp_));
  }

  void renderZYRectangle(GLfloat ll_z, GLfloat ll_y, GLfloat ur_z, GLfloat ur_y, 
			 GLfloat theta =0,bool filled=false) {
    std::cerr<<"rendering rectangle"<<std::endl;
    GLfloat degree_theta=180.0*theta/M_PI;
    GLfloat a=ur_z-ll_z;
    GLfloat b=ur_y-ll_y;
    GLfloat zc=0.5f*(ur_z+ll_z);
    GLfloat yc=0.5f*(ur_y+ll_y);
    
    glm::mat4  M=mvp_;
    std::cerr<<M<<std::endl;
    M=glm::translate(M,glm::vec3(0,yc,zc));
    M=glm::rotate(M,degree_theta,glm::vec3(-1,0,0));
    M=glm::scale(M,glm::vec3(0,b,a));



    std::cerr<<M<<std::endl;
     mvp_matrix_uniform_.load_4x4_matrix(shader_->shader(),
					glm::value_ptr(M));
     
     glBegin(GL_LINE_LOOP);
     glVertex3f(0.0f,-0.5f,-0.5f);
     glVertex3f(0.0f, 0.5f,-0.5f);
     glVertex3f(0.0f, 0.5f, 0.5f);
     glVertex3f(0.0f,-0.5f, 0.5f);
     glEnd();
     if(filled) {
       glBegin(GL_POLYGON);
       glVertex3f(0.0f,-0.5f,-0.5f);
       glVertex3f(0.0f, 0.5f,-0.5f);
       glVertex3f(0.0f, 0.5f, 0.5f);
       glVertex3f(0.0f,-0.5f, 0.5f);
     glEnd();
     }
     mvp_matrix_uniform_.load_4x4_matrix(shader_->shader(),
					glm::value_ptr(mvp_));
  }

  void renderEnd() {
    glEnd();
    glUseProgram(current_shader_);
  }

private:
 
  Shader *shader_;
  GLint current_shader_;
  UniformMatrix<4,4> mvp_matrix_uniform_;
  
  UniformVector<GLfloat,4> start_color_uniform;
  UniformVector<GLfloat,4> end_color_uniform;
  glm::mat4 mvp_;
  glm::mat4 mv_;

  
};



#endif
