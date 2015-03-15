uniform float _t,_p;
uniform vec2 _r;
uniform sampler2D _N,_T,_P,_F;
varying vec2 V;

vec4 noise(vec2 p){return texture2D(_N,(p+.5)/1024.,-20.);}
vec4 terrain(vec2 p_m){return texture2D(_T,(p_m+vec2(.5))/4096.,-20.);}
vec4 plan(vec2 p_m){return texture2D(_P,(floor(p_m/32.)+vec2(.5))/128.,-20.);}
vec2 fc=floor(gl_FragCoord.xy);
//int rand_state_=fc.y+fc.x*720.+int(_t)*141437;//+(noise(vec2(_t,0.)).x*256.+noise(vec2(0.,_t)).z)*256.;
int rand_state_=int(fc.y+fc.x*720.)+int(_t)*1023;//+(noise(vec2(_t,0.)).x*256.+noise(vec2(0.,_t)).z)*256.;
vec4 rand(){rand_state_=int(mod(float(rand_state_+1),1024.*1024.));return noise(vec2(float(rand_state_),floor(float(rand_state_)/1024.)));}

#define STEPS 128
#define HIT_EPS .001
#define FAR 3000.
#define BOUNCES 3
#define SKY FAR
#define GRIDSIZE 32.

vec3 sundir = normalize(vec3(.1,.06,.07));

float h(vec2 p){
  vec4 c=terrain(p);
  return c.z;}
float h2(vec2 p){
  vec4 c=terrain(p);
  return c.z;
}

vec3 terrain_normal(vec2 p){
  vec2 e=vec2(.5,.0);
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

float dist_terrain(vec3 p,vec3 D) {
  float H = h(p.xz);
  return p.y-H;
}

float dist_block(vec3 p, vec4 P) {
  vec2 cell = floor(p.xz/GRIDSIZE);
  vec2 cellc = cell*GRIDSIZE + GRIDSIZE/2.;
  vec3 cellcenter = vec3(cellc.x, h(cellc), cellc.y);
  vec4 crand = noise(cell);
  crand.y = mix(crand.y, -100., step(.5, crand.w));
  return 100.;
 //   dbox(p-cellcenter, vec3((.1+crand.x*.35)*GRIDSIZE*.5, crand.y*(P.y-P.z), (.05+crand.z*.4)*GRIDSIZE*.5));
}

vec3 block_normal(vec3 p, vec4 P){
  vec2 e=vec2(.01,.0);
  return normalize(vec3(
    dist_block(p+e.xyy,P)-dist_block(p-e.xyy,P),
    dist_block(p+e.yxy,P)-dist_block(p-e.yxy,P),
    dist_block(p+e.yyx,P)-dist_block(p-e.yyx,P)
    ));
}

//struct material_t {vec3 d,s,e;}

struct hit_t {
  int mid;
  float l;
  vec3 p,i,n;
  float _gc,_mc;
};

hit_t hit_block(hit_t h, vec4 P) {
  h.mid = 2;
  h.n = block_normal(h.p, P);
  return h;
}

hit_t hit_terrain(hit_t h) {
  h.mid = 1;
  h.n = terrain_normal(h.p.xz);
  return h;
}

#define ASSERT(cond,r,g,b) if(!(cond)){c=vec4(r,g,b,1.);break;}

float minpos(float a, float b) {
  return (a<0.)?b:((b<0.)?a:min(a,b));
}

float xp(float o, float d, vec2 p) {
  if (abs(d) < 1e-6) return 1e6;
  return ((d>0.) ? (p.y-o) : (p.x-o)) / d + HIT_EPS*2.;
}

#if 0
void main() {
  vec2 res=vec2(1280.,720.);
  vec2 uv=gl_FragCoord.xy/res-vec2(.5);uv.x*=res.x/res.y;
  vec3 O=vec3(16.,150.,300.),D=normalize(vec3(uv,-2.));
  O.y+=h2(O.xz);
  // FIXME replace with lens/sampler model
  //O+=noise(vec2(_t)).xyz*vec3(.3,.3,20.);
  vec3 c = vec3(0.);
#endif
hit_t trace_grid(vec3 O, vec3 D, float Lmax) {
  hit_t h;
  h.mid = -1;
  h.l = 0.;
  h.i = D;
  h._gc = h._mc = 0.;
  float dl = 0.;
  vec4 P=vec4(1e6);
  for (int i = 0; i < 256; ++i) {
    if (h.l > Lmax) break;
    h.p = O + D * h.l;
    if (dl < HIT_EPS || h.p.y > P.y) {
      h._gc += 1.;
      h.l += dl;
      h.p = O + D * h.l;
      vec2 cp = floor(h.p.xz/GRIDSIZE)*GRIDSIZE;
      float dx = xp(h.p.x,D.x,vec2(cp.x,cp.x+GRIDSIZE));
      float dz = xp(h.p.z,D.z,vec2(cp.y,cp.y+GRIDSIZE));
      //ASSERT(D.z<0.,0.,1.,0.)
      //ASSERT(dz > 0.,0.,1.,1.)
      //ASSERT(dx > 0.,0.,.5,1.)
      //ASSERT(h.p.z > cp.y,1.,1.,0.)
      //ASSERT(h.p.z < cp.y+GRIDSIZE,1.,0.,1.)
      P = plan(h.p.xz);
      //P(#,B,h,H)
      float dy = (h.p.y>P.y) ? xp(O.y,D.y,vec2(P.y,SKY)) :
       ((h.p.y>P.w) ? xp(O.y,D.y,vec2(P.w,P.y)) : xp(O.y,D.y,vec2(P.z,P.w)));
      //ASSERT(dy > 0.,.5,.5,1.)
      dl = min(dx,min(dy,dz));
      //if (dl < 0.) {c=vec3(1.,float(i),h.l);break;}
    } else {
      h._mc += 1.;
      float d = dist_block(h.p, P);
      if (d < HIT_EPS) return hit_block(h, P);
      if (h.p.y < P.w)
      {
        d = min(d, dist_terrain(h.p, D));
        if (d < HIT_EPS) return hit_terrain(h);
      }
      dl -= d;
      h.l += d;
    }
  }
  return h;
//  c = h.l/5000.;
// gl_FragColor=vec4(c,1.);
}
#if 1

struct sinfo_t {
  vec3 e,a;
};

sinfo_t solid_brdf(hit_t h, vec3 v) {
  float df = max(0.,dot(h.n,v));
  sinfo_t s;s.e=s.a=vec3(0.);
  switch (h.mid) {
  case 1:
    s.e = vec3(0.);
    s.a = vec3(.5);
    break;
  case 2:
    s.e = vec3(.5);
    s.a = vec3(1.);
    break;
  }
  s.a *= df;
  return s;
}

vec3 solid_bounce(hit_t h){
  vec3 r=rand().zyx*2.-vec3(1.);
  return normalize(r*sign(dot(r,h.n)));
}

//vec3 brdf(vec3 p, vec3 i, vec3 n, vec3 d){
//  mat3 m = geometry_material(p);
//  return m[0] * max(0.,dot(n,d));
//}

//DEF_TRACE(trace_geometry,geometry_world,64,.01)

vec3 air(vec3 O, vec3 D) {
  return 30. * vec3(1.) * step(.99,dot(D,sundir)) + vec3(.1);
}

#define MAKE_Q(T) T quantize(T a,T b,int n,T v){return floor(float(n)*(v-a)/(b-a));}
MAKE_Q(vec3)
MAKE_Q(float)

void main(){
  vec2 res=vec2(1280.,720.);
  vec2 uv=gl_FragCoord.xy/res-vec2(.5);uv.x*=res.x/res.y;
  vec3 O=vec3(16.,150.,300.),D=normalize(vec3(uv,-2.));
  //vec3 O=vec3(sin(t*.01)*100.,50.,cos(t*.01)*300.),D=normalize(vec3(uv,-2.));
  //vec3 O=vec3(sin(t)*1000.,50.,cos(t)*1000.),D=normalize(vec3(uv,-2.));
  O.y+=h2(O.xz);

  // FIXME replace with lens/sampler model
  //O+=noise(vec2(_t)).xyz*vec3(.3,.3,20.);

  vec3 color = vec3(0.), transm=vec3(1.);
#if 1
  hit_t t = trace_grid(O, D, FAR);
  //gl_FragColor = vec4(quantize(vec3(0.),vec3(256.),16,vec3(h._gc, h._mc, h._mc+h._gc)), 1.);return;
  //gl_FragColor = vec4(vec3(quantize(0.,256.,8,h._gc+h._mc)), 1.);return;
  {
    float dh=t.p.y-h(t.p.xz);
    if (dh>HIT_EPS){gl_FragColor=vec4(1.,0.,0.,1.)*100.;return;}
    if (dh<0.){gl_FragColor=vec4(1.,0.,1.,1.)*100.;return;}
  }

  if (t.mid > 0) {
    color = t.n*4.;
  } else {
    //color = h.l/FAR;//air(O, D);
    //color = air(O, D);
    color = vec3(1000.,0.,0.);
  }
#else
  for (int i=0;i<1;++i){
    if (dot(transm,transm) < .001) break;
    hit_t h = trace_grid(O, D, FAR);

    if (h.mid > 0) {
      //if (p.y-h(p.xz)>HIT_EPS*lt){gl_FragColor=vec4(1.,0.,0.,1.);return;}
      //if (h(p.xz)>p.y){gl_FragColor=vec4(1.,0.,1.,1.);return;}
      //mat3 m = geometry_material(p);
      //vec3 c = m[1];
      vec3 c = vec3(0.);
      O = h.p + h.n * 2. * HIT_EPS;
      // importance
      if (trace_grid(O, sundir, 100.).l >= 100.) c += solid_brdf(h, sundir).a * air(O, sundir);
      vec3 nD = solid_bounce(h);
      color += transm * c;
      transm *= solid_brdf(h, nD).a;
      D = nD;
    } else {
      color += transm * air(O, D);
    }
  }
  //color = max(vec3(0.), color);
#endif
  gl_FragColor = vec4(pow(color,vec3(1./2.2)),1.);
}
#endif
