uniform float _t, _p, _pt;
uniform vec2 _r;
uniform sampler2D _N,_T,_P,_F;
varying vec2 V;

vec4 noise(vec2 p){return texture2D(_N,(p+.5)/1024.,-20.);}
vec4 terrain(vec2 p_m){return texture2D(_T,p_m/4096.,-20.);}
vec4 plan(vec2 p_m){return texture2D(_P,(floor(p_m/32.)+.5)/128.,-20.);}

float t=_t;
float h(vec2 p){
  vec4 c=texture2D(_T,p/4096.,-20.);
  return c.z;
  //return c.z;
}
float h2(vec2 p){
  vec4 c=terrain(p);
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
#elif 0
  float XZ = mod(floor(p.x/1.)+floor(p.z/1.),2.);
  float Y=mod(floor(p.y),2.);
  return vec3(.2+.5*Y+.3*XZ);
#else
  vec4 T=terrain(p.xz);
  return vec3(.2) + vec3(T.x, T.y, T.w);
#endif
}
vec4 trace(vec2 uv){
  vec3 O=vec3(0.,150.,300.),D=normalize(vec3(uv,-2.));
  O.y+=h2(O.xz);
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
  if(l>3000.){return vec4(0.);}
  vec3 p=O+D*l;
  if (p.y-h(p.xz)>.01*l){return vec4(1.,0.,0.,1.);}
  if (h(p.xz)>p.y){return vec4(1.,0.,1.,1.);}
  vec3 n = n((O+D*l).xz);
  vec3 color=albedo(p)*(max(0.,dot(n, normalize(vec3(1., 1., -.5))))+vec3(.05));
  return vec4(pow(color,vec3(1./2.2)),1.);
}

void main(){
  if (gl_FragCoord.y<16.){gl_FragColor=vec4(step(V.x*.5+.5,_p));return;}
#if 0
  vec2 uv=floor(gl_FragCoord.xy)/_r;
  vec2 pt=floor(uv*4096.);
  vec2 pp=floor(uv*128.);
  vec4 T=texture2D(_T,(pt+vec2(.5))/4096.,-20.);
#if 1
  gl_FragColor=T*vec4(1.,1.,.005,1.);
#else
  vec4 P=texture2D(_P,(pp+vec2(.5))/128.,-20.);
  float dH = P.w - T.z;
  float dh = T.z - P.z;
  gl_FragColor=vec4(-dH, -dh, (P.w-P.z)/100., 0.);
#endif
#else
  if (_pt < .0)
    gl_FragColor = trace((gl_FragCoord.xy/_r - vec2(.5)) * vec2(_r.x/_r.y,1.));
  else
    gl_FragColor = texture2D(_F,gl_FragCoord.xy/_r,-20.)*.2/(1.+256.*_pt);
#endif
}
