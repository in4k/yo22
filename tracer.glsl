uniform float _t, _p;
uniform vec2 _r;
uniform sampler2D _N,_T,_F;
float h(vec2 p){return texture2D(_T,p/4096.,-20.).w;}
vec3 n(vec2 p){
  vec2 e=vec2(1.,0.);
  vec3 dx=vec3(2.,h(p+e.xy)-h(p-e.xy),0.);
  vec3 dz=vec3(0.,h(p+e.yx)-h(p-e.yx),2.);
  return normalize(cross(dz,dx));
}
void main(){
  vec2 uv=V;uv.x*=_r.x/_r.y;
  //gl_FragColor = vec4(sin(_t));//texture2D(_N,V,-20.);return;
  //gl_FragColor = vec4(_p/1024.);return;
  vec3 O=vec3(0,1000.,0.),D=normalize(vec3(uv,-2.));
  float l=0.;
  for(int i=0;i<40;++i){
    vec3 p=O+D*l;
    float d=p.y-h(p.xz);
    if(d<1.)break;
    l+=d*1.5;
  }
  if(l>2000.){gl_FragColor=vec4(0.);return;}
  //gl_FragColor=l/2000.;//vec4((O+D*l).y/1000.);
  vec3 color=max(0.,dot(n(O+D*l), normalize(vec3(1.))));
  //vec3 color=((O+D*l)-980.)/20.;
  gl_FragColor=vec4(_p);
}
