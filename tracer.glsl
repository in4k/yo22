uniform float _t,_p;
uniform vec2 _r;
uniform sampler2D _N,_T,_P,_F;
varying vec2 V;

vec4 noise(vec2 p){return texture2D(_N,(p+.5)/1024.,-20.);}
vec4 terrain(vec2 p_m){return texture2D(_T,p_m/4096.,-20.);}
vec4 plan(vec2 p_m){return texture2D(_P,(floor(p_m/32.)+.5)/128.,-20.);}
vec2 fc=floor(gl_FragCoord);
int rand_state_=fc.y+fc.x*720.+int(_t)*141437;//+(noise(vec2(_t,0.)).x*256.+noise(vec2(0.,_t)).z)*256.;
vec4 rand(){rand_state_=mod(rand_state_+1,1024*1024);return noise(vec2(float(rand_state_),floor(float(rand_state_)/1024.)));}

#define STEPS 128
#define HIT_EPS .001
#define FAR 3000.
#define BOUNCES 3

vec3 sundir = normalize(vec3(.1,.06,.07));

float h(vec2 p){
  vec4 c=terrain(p);
  return c.z;}
float h2(vec2 p){
  vec4 c=terrain(p);
  return c.z;
}

vec3 terrain_normal(vec2 p){
  vec2 e=vec2(.5,0.);
  vec3 dx=vec3(2.*e.x,h(p+e.xy)-h(p-e.xy),0.);
  vec3 dz=vec3(0.,h(p+e.yx)-h(p-e.yx),2.*e.x);
  return normalize(cross(dz,dx));
}

vec3 terrain_albedo(vec3 p){
#if 1
  //return step(80.,terrain(p.xz).z);
  vec4 P=plan(p.xz);
  //return P.g;//(P.w-P.z)/200.;
  float d=P.w-P.z;
  return vec3(.2, .2 + 1.-min(1.,floor(d/5.)/2.), .2);
  //return 1.-step(8., P.w-P.z);
  //return step(100.,P.w);
#elif 1
  return vec3(.2);
#elif 1
  float m10 = mod(floor(p.x/10.)+floor(p.z/10.),2.);
  float m100 = mod(floor(p.x/100.)+floor(p.z/100.),2.);
  float m1000 = mod(floor(p.x/1000.)+floor(p.z/1000.),2.);
  return vec3(m10,m100,m1000)*(.3+.7*mod(floor(p.x)+floor(p.z),2.));
#elif 1
  float XZ = mod(floor(p.x/1.)+floor(p.z/1.),2.);
  float Y=mod(floor(p.y),2.);
  return vec3(.2+.5*Y+.3*XZ);
#else
  return texture2D(_T,p.xz/4096.,-20.).xyw*10.;
#endif
}
//vec3 terrain_brdf(float s,vec3 p,vec3 n,vec3 inc,vec3 dir){}
vec3 terrain_bounce(float s,vec3 p,vec3 n,vec3 inc){
  vec3 r=noise(vec2(s)).ywx*2.-vec3(1.);
  return normalize(r*sign(dot(r,n)));
}

#define DEF_TRACE(NAME_,FUNC_,STEPS_,HIT_EPS_) \
float NAME_(vec3 O, vec3 D, float l, float lM){\
  float lp=l;\
  for(int i=0;i<STEPS_;++i){\
    vec3 p=O+D*l;\
    float dhit=HIT_EPS*l;\
    float d=FUNC_(p,D,dhit);\
    if(d<dhit){\
      for (int j=0;j<9;++j){\
        float lm=(l+lp)*.5;\
        p=O+D*lm;\
        d=FUNC_(p,D,dhit);\
        if(d<dhit)l=lm;else lp=lm;\
      }\
      break;\
    }\
    lp=l;\
    l+=d;\
    if(l>lM)return lM;\
  }\
  return l;\
}

float maxv(vec3 v){return max(v.x,max(v.y,v.z));}
float dbox(vec3 p,vec3 sz){
  return maxv(abs(p)-sz);
}

float dist_b0(vec3 p){
  vec3 ofs = vec3(0., 250., -500.);
  return length(p-ofs)-30.;
}

float dist_b1(vec3 p){
  vec3 ofs = vec3(0., 230., -500.);
  return length(p-ofs-vec3(85.,0.,0.))-24.;
}
float dist_b2(vec3 p){
  vec3 ofs = vec3(0., 230., -500.);
  return length(p-ofs+vec3(85.,0.,0.))-20.;
}

float dist_terrain(vec3 p,vec3 D,float dhit){
  float H = h(p.xz);
  return p.y-H;
}

float dist_city(vec3 p, vec3 D, float dhit){
  float SZ = 32.;
  vec2 cell = floor(p.xz/SZ);
  vec2 cellc = cell*SZ + SZ/2.;
  vec3 cellcenter = vec3(cellc.x, h(cellc), cellc.y);
  vec4 crand = noise(cell);
  crand.y = mix(crand.y, -100., step(.5, crand.w));
  return //100.;
    dbox(p-cellcenter, vec3((.1+crand.x*.35)*SZ*.5, crand.y*50., (.05+crand.z*.4)*SZ*.5));
}

mat3 geometry_material(vec3 p){
  float terrain = dist_terrain(p,vec3(0.),0.);
  float bld = dist_city(p,vec3(0.),0.);
  float b0 = dist_b0(p), b1 = dist_b1(p), b2 = dist_b2(p);
  float mm = min(min(min(min(terrain,b0),b1),b2),bld);

  mat3 m = mat3(0.,0.,0.,0.,0.,0.,0.,0.,0.);
  if (mm == bld) {
    vec2 cell = floor(p.xz/32.);
    vec4 crand = noise(cell);
    float ff = mod(p.y/4.,1.);
    vec3 c = (vec3(ff*.0)+1.*noise(floor(crand.yx*1024.)).wyx);
    m[0] = c * .9;
    m[1] = c * 1.2;
  } else if (mm == b2) {
    m[0] = vec3(1.);
    m[1] = vec3(100.,0.,0.);
  } else if (mm == b1) {
    m[0] = vec3(1.);
    m[1] = vec3(0.,100.,0.);
  } else if (mm == b0) {
    m[0] = vec3(1.);
    m[1] = vec3(1000.);
  } else {
    m[0] = .6*terrain_albedo(p);
  }

  return m;
}

float geometry_world(vec3 p,vec3 D,float dhit){
  float terrain = dist_terrain(p,D,dhit);
  float bld = dist_city(p,D,dhit);
  //return min(min(min(min(terrain,dist_b0(p)),dist_b1(p)),dist_b2(p)),bld);
  return min(terrain,bld);
}

float geometry_world_(vec3 p){return geometry_world(p,vec3(0.),.0);}

vec3 geometry_normal(vec3 p){
  vec2 e=vec2(.01,.0);
  return normalize(vec3(
    geometry_world_(p+e.xyy)-geometry_world_(p-e.xyy),
    geometry_world_(p+e.yxy)-geometry_world_(p-e.yxy),
    geometry_world_(p+e.yyx)-geometry_world_(p-e.yyx)
    ));
}

//struct hit_t {vec3 p,i,n;float l,s;};
//struct material_t {vec3 d,s,e;}

vec3 geometry_bounce(vec3 p,vec3 n,vec3 inc){
  vec3 r=rand().zyx*2.-vec3(1.);
  return normalize(r*sign(dot(r,n)));
}

vec3 brdf(vec3 p, vec3 i, vec3 n, vec3 d){
  mat3 m = geometry_material(p);
  return m[0] * max(0.,dot(n,d));
}

DEF_TRACE(trace_geometry,geometry_world,64,.01)

vec3 air(vec3 O, vec3 D) {
  return 30. * vec3(1.) * step(.99,dot(D,sundir)) + vec3(.1);
}

void main(){
  vec2 res=vec2(1280.,720.);
  vec2 uv=gl_FragCoord.xy/res-vec2(.5);uv.x*=res.x/res.y;
  vec3 O=vec3(0.,150.,300.),D=normalize(vec3(uv,-2.));
  //vec3 O=vec3(sin(t*.01)*100.,50.,cos(t*.01)*300.),D=normalize(vec3(uv,-2.));
  //vec3 O=vec3(sin(t)*1000.,50.,cos(t)*1000.),D=normalize(vec3(uv,-2.));
  O.y+=h2(O.xz);

  // FIXME replace with lens/sampler model
  O+=noise(vec2(_t)).xyz*.3;

  vec3 color = vec3(0.), cmask=vec3(1.);
  for (int i=0;i<BOUNCES;++i){
    if (dot(cmask,cmask) < .001) break;
    float lg = trace_geometry(O, D, 0., FAR);

    if (lg < FAR) {
      //if (p.y-h(p.xz)>HIT_EPS*lt){gl_FragColor=vec4(1.,0.,0.,1.);return;}
      //if (h(p.xz)>p.y){gl_FragColor=vec4(1.,0.,1.,1.);return;}
      vec3 p=O+D*lg;
      vec3 n = geometry_normal(p);
      mat3 m = geometry_material(p);
      O = p + n * 2. * HIT_EPS;
      vec3 c = m[1];
      // importance
      if (trace_geometry(O, sundir, 0., FAR) >= FAR) c += brdf(O, D, n, sundir) * air(O, sundir);
      vec3 nD = geometry_bounce(O, n, D);
      color += cmask * c;
      cmask *= brdf(O, D, n, nD);
      D = nD;
    } else {
      color += cmask * air(O, D);
    }
  }
  //color = max(vec3(0.), color);
  gl_FragColor = vec4(pow(color,vec3(1./2.2)),1.);
}
