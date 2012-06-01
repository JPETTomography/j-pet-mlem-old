#version 330

out vec4 f_color;
in  vec2 f_texture0;

uniform sampler2D texture0;



void main() {

     vec4 start_color=vec4(0,0,1,1);
     vec4 end_color  =vec4(1,0,0,1);	

       f_color=vec4(1,0,0,1);
       float f=texture(texture0,f_texture0).r;
       f_color=mix(start_color,end_color,f);
}