#version 150

out vec4 f_color;
in  vec2 f_texture0;

uniform sampler2D texture0;

uniform     vec4 start_color;
uniform     vec4 end_color;	


void main() {


       f_color=vec4(1,0,0,1);
       float f=texture(texture0,f_texture0).r;
       f_color=mix(start_color,end_color,f);
}