uniform sampler2D _N,_T;
varying vec2 V;
void main(){gl_FragColor=texture2D(_T,gl_FragCoord.xy/4096.,-20.);}
