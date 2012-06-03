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


#include"tof_event.h"
#include"tof_detector.h"
#include"phantom.h"
#include"tof_monte_carlo.h"

void SetupRC() {


  
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glEnable(GL_DEPTH_TEST);

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


#if 1
  TOPETSimulator<GLfloat>  simulator;
  simulator.init();
  simulator.simulate_from_single_point(10000);
#else

  PixelGrid<GLfloat> grid(Point<GLfloat>(-250,-250),
			  Point<GLfloat>(250,250),width,height);

  typedef ToF_Detector_2D<GLfloat> detector_t;
  typedef ToF_Event_2D<GLfloat> event_t;
  ToF_Detector_2D<GLfloat> detector(350,500);
  detector.set_sigma(11,32);
  ToF_Monte_Carlo<GLfloat,ToF_Detector_2D<GLfloat> > mc(detector);
  mc.gen_seeds(5565665);
  const int  n=100000;
  std::vector<event_t> events(n);
  std::vector<event_t>::iterator it=events.begin();
  int count;
  count=mc.fill_with_events_from_single_point(it,0,0,n/2);
  it=events.begin()+count;
  count+=mc.fill_with_events_from_single_point(it,200,-150,n/2);
  std::vector<event_t> events_out(n);

  count=mc.add_noise_to_detected(events.begin(),events.begin()+count,events_out.begin());

  for(int i=0;i<count;i++) {
    grid.insert(events_out[i].z(),
		events_out[i].y(),
		1.0f
		);
  }

  GLfloat max=grid.max();
  grid.divide_by(max);

  density_plot->set_pixmap(grid.pixels());
 #endif
  density_plot->set_pixmap(simulator.emitted_density());
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

  GLenum err = glewInit();
  if (GLEW_OK != err) {
        fprintf(stderr, "GLEW Error: %s\n", glewGetErrorString(err));
  }


  SetupRC();
  glutMainLoop();

  std::cout<<"after main loop"<<std::endl;
  return 0;
  
}

