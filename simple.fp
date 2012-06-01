#version 330

out vec4 f_color;
in  vec2 f_texture0;

uniform sampler2D texture0;


void main() {
       f_color=vec4(1,0,0,1);
       f_color=texture(texture0,f_texture0);
}