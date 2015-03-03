uniform float _t, _p;
uniform vec2 _r;
uniform sampler2D _N,_T,_F;
varying vec2 V;
float t=_t*.06;
float h(vec2 p){
  //return 400. + 100.*(sin(p.y*.01)+sin(p.x*.01));
  return texture2D(_T,p/4096.,-20.).w;
}
vec3 n(vec2 p){
  vec2 e=vec2(1.,0.);
  vec3 dx=vec3(2.*e.x,h(p+e.xy)-h(p-e.xy),0.);
  vec3 dz=vec3(0.,    h(p+e.yx)-h(p-e.yx),2.*e.x);
  return normalize(cross(dz,dx));
}
void main(){
  vec2 uv=gl_FragCoord.xy/_r-.5;uv.x*=_r.x/_r.y;
  vec3 O=vec3(sin(t)*1000.,150.,cos(t)*1000.),D=normalize(vec3(uv,-2.));
  O.y+=h(O.xz);
  float l=0.,lp=l;
  for(int i=0;i<128;++i){
    vec3 p=O+D*l;
    float d=p.y-h(p.xz);
    if(l>2000.)break;
    if(d<.01*l) {
      for (int j=0;j<9;++j){
        float lm=(l+lp)*.5;
        vec3 p=O+D*lm;
        float d=p.y-h(p.xz);
        if(d<.01*lm)l=lm;else lp=lm;
      }
      break;
    }
    lp=l;
    l+=max(20.,d*.1);
  }
  if(l>2000.){gl_FragColor=vec4(0.);return;}
  vec3 p=O+D*l;
  if (p.y-h(p.xz)>.01*l){gl_FragColor=vec4(1.,0.,0.,1.);return;}
  if (h(p.xz)>p.y){gl_FragColor=vec4(1.,0.,1.,1.);return;}
  //gl_FragColor=l/2000.;//vec4((O+D*l).y/1000.);
  vec3 n = n((O+D*l).xz);
  //vec3 color=vec3(-n.y, n.y, 0.);
  //vec3 color=max(0.,dot(n, normalize(vec3(1.))))+vec3(.2);
  vec3 color=max(0.,dot(n, normalize(vec3(cos(t*10.),1.,sin(t*20.)))))+vec3(.2);
  //vec3 color = (O+D*l) / vec3(3000.,1000.,3000.);
  gl_FragColor = mix(vec4(color,1.), vec4(step(V.x*.5+.5,_p)), step(V.y,-.96));
}
