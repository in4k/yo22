uniform float _t, _p;
uniform vec2 _r;
uniform sampler2D _N,_T,_F;
varying vec2 V;
float t=_t*.001;
float h(vec2 p){
  vec4 c=texture2D(_T,p/4096.,-20.);
  return c.z+c.w*10.;
  //return c.z;
}
float h2(vec2 p){
  vec4 c=texture2D(_T,p/4096.,-20.);
  //return c.z+c.w;
  return c.z;
}
vec3 n(vec2 p){
  vec2 e=vec2(.5,0.);
  vec3 dx=vec3(2.*e.x,h(p+e.xy)-h(p-e.xy),0.);
  vec3 dz=vec3(0.,h(p+e.yx)-h(p-e.yx),2.*e.x);
  return normalize(cross(dz,dx));
}
vec3 albedo(vec3 p){
#if 0
  float m10 = mod(floor(p.x/10.)+floor(p.z/10.),2.);
  float m100 = mod(floor(p.x/100.)+floor(p.z/100.),2.);
  float m1000 = mod(floor(p.x/1000.)+floor(p.z/1000.),2.);
  return vec3(m10,m100,m1000)*(.3+.7*mod(floor(p.x)+floor(p.z),2.));
#else
  return texture2D(_T,p.xz/4096.,-20.).xyw*10.;
#endif
}
void main(){
  if (gl_FragCoord.y<16.){gl_FragColor=step(V.x*.5+.5,_p);return;}
  vec2 uv=gl_FragCoord.xy/_r-.5;uv.x*=_r.x/_r.y;
  vec3 O=vec3(sin(t)*1000.,50.,cos(t)*1000.),D=normalize(vec3(uv,-2.));
  //for(int i=0;i<16.
  O.y+=h2(O.xz);
  //vec3 O=vec3(sin(t)*1000.,600.,cos(t)*1000.),D=normalize(vec3(uv,-2.));
  float l=0.,lp=l;
  for(int i=0;i<128;++i){
    vec3 p=O+D*l;
    float d=p.y-h(p.xz);
    if(d<.001*l) {
      for (int j=0;j<9;++j){
        float lm=(l+lp)*.5;
        vec3 p=O+D*lm;
        float d=p.y-h(p.xz);
        if(d<.001*lm)l=lm;else lp=lm;
      }
      break;
    }
    lp=l;
    //l+=max(20.,d*.1);
    l+=d;
    if(l>3000.)break;
  }
  if(l>3000.){gl_FragColor=vec4(0.);return;}
  vec3 p=O+D*l;
  if (p.y-h(p.xz)>.01*l){gl_FragColor=vec4(1.,0.,0.,1.);return;}
  if (h(p.xz)>p.y){gl_FragColor=vec4(1.,0.,1.,1.);return;}
  //gl_FragColor=l/2000.;//vec4((O+D*l).y/1000.);
  vec3 n = n((O+D*l).xz);
  //vec3 color=vec3(-n.y, n.y, 0.);
  vec3 color=albedo(p)*(max(0.,dot(n, normalize(vec3(1., 1., -.5))))+vec3(.05));
  //vec3 color=max(0.,dot(n, normalize(vec3(cos(t*10.),1.,sin(t*20.)))))+vec3(.2);
  //vec3 color = (O+D*l) / vec3(3000.,1000.,3000.);
  gl_FragColor = vec4(pow(color,vec3(1./2.2)),1.);
}
