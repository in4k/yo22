uniform float P,t;
uniform vec2 S;
uniform sampler2D R,T;
varying vec2 V;
float h(vec2 p){return texture2D(T,p/4096.,-20.).w*1000.;}
void main(){
  vec2 uv=V;//uv.x*=S.x/S.y;
  vec3 O=vec3(0,1000.,0.),D=normalize(vec3(uv,-2.));
  float l=0.;
  for(int i=0;i<64;++i){
    vec3 p=O+D*l;
    float d=p.y-h(p.xz);
    if(d<.01)break;
    l+=d*.5;
  }
  //gl_FragColor=l/2000.;//vec4((O+D*l).y/1000.);
  vec3 color=(O+D*l)/2000.;
  gl_FragColor = mix(vec4(color,1.), vec4(step(V.x*.5+.5,P)), step(V.y,-.96));
}
