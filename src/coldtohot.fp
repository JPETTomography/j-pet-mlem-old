#version 330

out vec4 f_color;
in  vec2 f_texture0;

uniform sampler2D texture0;



void main() {

     vec4 blue=vec4(0,0,1,1);
     vec4 cyan=vec4(0,1,1,1);
     vec4 green=vec4(0,1,0,1);
     vec4 yellow=vec4(1,1,0,1);
     vec4 red=vec4(1,0,0,1);	


       float f=texture(texture0,f_texture0).r;
       if(f<0.25)	
       	 	f_color=mix(blue,cyan,4*f);
       else if(f<0.5)
       	    	f_color=mix(cyan,green,4*(f-0.25));
       else if(f<0.75)
       	    	f_color=mix(green,yellow,4*(f-0.50));
       else 
       	    	f_color=mix(yellow,red,4*(f-0.75));
}