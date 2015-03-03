uniform sampler2D R,T;
varying vec2 V;
void main(){gl_FragColor=texture2D(T,gl_FragCoord.xy/4096.,-20.);}
