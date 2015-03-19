uniform float _t,_p;
uniform vec2 _r;
uniform vec3 _s,_cp;
uniform mat3 _cm;
uniform sampler2D _N,_T,_P,_F;
varying vec2 V;

vec4 noise(vec2 p){return texture2D(_N,(p+.5)/1024.,-20.);}
vec4 terrain(vec2 p_m){return texture2D(_T,(p_m+vec2(.5,.5))/4096.,-20.);} // FIXME why -.5???
vec4 plan(vec2 p_m){return texture2D(_P,(floor(p_m/32.)+vec2(.5))/128.,-20.);}
vec2 fc=floor(gl_FragCoord.xy);
//int rand_state_=fc.y+fc.x*720.+int(_t)*141437;//+(noise(vec2(_t,0.)).x*256.+noise(vec2(0.,_t)).z)*256.;
//int rand_state_=int(fc.y*1280.+fc.x)+int(_t)*1023;//+(noise(vec2(_t,0.)).x*256.+noise(vec2(0.,_t)).z)*256.;
int rand_state_=int(fc.y*1280.+fc.x)+int(noise(vec2(_t,0.)).x*256.+noise(vec2(0.,_t)).z)*256;
vec4 rand(){rand_state_=int(mod(float(rand_state_+1),1024.*1024.));return noise(vec2(float(rand_state_),floor(float(rand_state_)/1024.)));}
vec4 n4(vec2 v){return texture2D(_N,(v+vec2(.5))/256.,-20.);}
float n(vec2 v){return n4(v).w;}
float vn(vec2 v,vec2 m) {
  vec2 e=vec2(1.,0.),V=floor(v);v=fract(v);v*=v*(3.-2.*v);
  return mix(mix(n(mod(v+e.yy,m)),n(mod(v+e.xy,m)),v.x),mix(n(mod(v+e.yx,m)),n(mod(v+e.xx,m)),v.x),v.y);
}
float fbm(vec2 v,float s) {
  float r=0.,k=.5;
  for(int i=0;i<4;++i,k*=.5,s*=2.)r+=k*vn(v*s,vec2(s));
  return r;
}

#define STEPS 256
#define EPS .01
#define FAR 5000.
#define BOUNCES 3
#define SKY FAR
#define GRIDSIZE 32.

vec3 sundir = normalize(_s);//vec3(.1,.036,.037));

float H(vec2 p){
  vec4 c=terrain(p);
  return c.z;}
float h2(vec2 p){
  vec4 c=terrain(p);
  return c.z;
}

vec3 terrain_normal(vec2 p){
  vec2 e=vec2(1.,.0);
  vec3 dx=vec3(2.*e.x,H(p+e.xy)-H(p-e.xy),0.);
  vec3 dz=vec3(0.,H(p+e.yx)-H(p-e.yx),2.*e.x);
  //return normalize(cross(dz,dx));//+1.5*noise(p*134.241));
  return normalize(cross(dz,dx)+3.5*(noise(p*7.31).yxz-vec3(.5)));
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

float maxv(vec3 v){return max(v.x,max(v.y,v.z));}
float dbox(vec3 p,vec3 sz){
  return maxv(abs(p)-sz);
}

float dist_terrain(vec3 p,vec3 D) {
  float H = H(p.xz);
  return (p.y - H);
}


#define MAKE_BINSEARCH(NAME_,DISTF_,STEPS_,ARG_) \
float NAME_(vec3 O, vec3 D, float e, float lo, float li){\
  for(int i = 0; i < STEPS_; ++i) {\
    float lm = (lo + li) * .5;\
    if (DISTF_(O + D * lm, ARG_) < e) li = lm; else lo = lm;\
  }\
  return lo;\
}
MAKE_BINSEARCH(bin_terrain, dist_terrain, 8, D)

struct cell_t {
  vec2 i;
  vec3 c,m,M,BB;
  vec4 S;
  float B,H,h,p,Bh,f;
};

cell_t cell_identify(vec3 p) {
  cell_t c;
  c.i = floor(p.xz/GRIDSIZE);
  c.m = c.i.xyy*GRIDSIZE;
  vec4 P = plan(c.m.xz); c.p = P.x;
  c.B = P.y; c.H = P.w; c.h = P.z; c.Bh = c.B - c.H;
  c.m.y = c.h;
  c.M = c.m + vec3(GRIDSIZE, c.B, GRIDSIZE);
  c.c = c.m + vec3(GRIDSIZE, c.H-c.h, GRIDSIZE) * .5;
  c.S = noise(c.i+vec2(23.,17.)).wyzx;
  c.BB = vec3(GRIDSIZE*.5-4., c.Bh, GRIDSIZE*.5-4.);
  c.f = 9.+.2*c.S.y;
  return c;
}

float bd1(vec3 p, cell_t c) {
vec2 e=vec2(1.,0.);
  float d = 1e6;
  vec3 s = c.BB * vec3(.2+.3*c.S.z,1.,.2+.3*c.S.y);
  for (int i = 0; i < 4; ++i)
  {
    vec3 pp = p-(c.BB-s)*(2.*vec3(c.S.w,0.,c.S.z)-vec3(1.));
    d = min(d, dbox(pp, s));
    s += vec3(.2+c.S.x, 0., .2+c.S.y) * .2 * c.BB;
    s.y *= .7+.1*c.S.z;
  }
  return d;
}

float dist_block(vec3 p, cell_t c) {
  p -= c.c; // p is now in cell-space
  float d = (c.Bh < 10.) ? 1e6 : bd1(p, c);
  return max(d, dbox(p, c.BB));
}

vec3 block_normal(vec3 p, cell_t c){
  vec2 e=vec2(.01,.0);
  return normalize(vec3(
    dist_block(p+e.xyy,c)-dist_block(p-e.xyy,c),
    dist_block(p+e.yxy,c)-dist_block(p-e.yxy,c),
    dist_block(p+e.yyx,c)-dist_block(p-e.yyx,c)
    ));
}

struct hit_t {
  int mid;
  float l;
  vec3 p, i, n;
  cell_t c;
  float _gc, _mc; // DEBUG
};

hit_t hit_block(hit_t h) {
  h.mid = 2;
  h.n = block_normal(h.p, h.c);
  return h;
}

hit_t hit_terrain(hit_t h) {
  h.mid = (h.c.Bh > 10.) ? 3 : 1;
  h.n = terrain_normal(h.p.xz);
  return h;
}

vec4 assertcolor_=vec4(0.);
#define ASSERT(cond,r,g,b) if(!(cond)){assertcolor_=vec4(r,g,b,1.);}

float minpos(float a, float b) {
  return (a<0.)?b:((b<0.)?a:min(a,b));
}

float xp(float o, float d, float m, float M) {
  if (abs(d) < 1e-6) return 1e6;
  return ((d>0.) ? (M-o) : (m-o)) / d;
}

hit_t trace_grid(vec3 O, vec3 D, float Lmax) {
  hit_t h;
  h.mid = -1;
  h.l = 0.;
  h.i = D;
  h._gc = h._mc = 0.; // DEBUG
  float dl = 0.;
  float pl = h.l;
  for (int i = 0; i < STEPS; ++i) {
    if (h.l > Lmax) break;
    float de = EPS * h.l * 2e-2;
    h.p = O + D * h.l;
    if (dl < 0. || h.p.y > h.c.B) {
      h._gc += 1.;
      h.l += max(0.,dl) + EPS;
      h.p = O + D * h.l;
      h.c = cell_identify(h.p);
      // TODO skipsize dependent on D.y
      //float dx = xp(h.p.x, D.x, h.c.m.x - GRIDSIZE, h.c.m.x + GRIDSIZE);
      //float dz = xp(h.p.z, D.z, h.c.m.z - GRIDSIZE, h.c.m.z + GRIDSIZE);
      float dx = xp(h.p.x, D.x, h.c.m.x, h.c.m.x + GRIDSIZE);
      float dz = xp(h.p.z, D.z, h.c.m.z, h.c.m.z + GRIDSIZE);
      //ASSERT(dz >= 0.,0.,1.,1.)
      //ASSERT(dx >= 0.,0.,.5,1.)
      //ASSERT(h.p.z >= cp.y,1.,1.,0.)
      //ASSERT(h.p.z <= cp.y+GRIDSIZE,1.,0.,1.)
      float dy =
        (h.p.y > h.c.B) ?
          xp(h.p.y, D.y, h.c.B, SKY) :
        ((h.p.y > h.c.H) ?// && (P.y-P.w) > 0.) ?
          xp(h.p.y, D.y, h.c.H, h.c.B) :
          1e6);//xp(h.p.y, D.y, h.c.h, h.c.H));
      //ASSERT(dy >= 0.,.5,.5,1.)
      dl = min(dx,min(dy,dz));
      //ASSERT(dl >= 0.,.5,.5,.5)
      dl += EPS;
      //if (dl < 0.) {c=vec3(1.,float(i),h.l);break;}
    } else
    //h._mc -= 1.;}
    {
      h._mc += 1.;
      float d = dist_block(h.p, h.c);
      if (d < de) return hit_block(h);
      if (h.p.y < h.c.H)
      {
        d = min(d, dist_terrain(h.p, D));
        if (d < de) {
          h.l = bin_terrain(O, D, de, pl, h.l);
          h.p = O + D * h.l;
          return hit_terrain(h);
        }
      }
      d = min(d,dl+EPS);
      dl -= d;
      pl = h.l;
      h.l += d;
    }
  }
  return h;
}

struct mat_t {
  vec3 e, Cd, Cs;
  float s;
};

vec3 mg(vec2 p) {
  return vec3(.2,.6,.23)*.3 + .1 * (noise(p*.2).xxx-vec3(.5));
}

mat_t material(hit_t h) {
  mat_t m;
  m.e = vec3(1., 0., 1.);
  m.Cd = vec3(0.);
  m.Cs = vec3(0.);
  m.s = 0.;
  vec3 cp = abs(abs(h.p - h.c.c) - vec3(GRIDSIZE,0.,GRIDSIZE)*.5);

  if(h.mid == 1) {
    float f = step(h.p.y+h.c.S.w*50.+20.*noise(h.p.xz).x,180.);
    m.e = vec3(0.);
    m.Cd = mg(h.p.xz) + vec3(.1,.1,.02)*h.c.S.xyz*f;
    m.Cd = mix(m.Cd, vec3(.11, .07, .03), .4*f*step(min(cp.x,cp.z), 1.));
    //if (h.p.y < H(h.p.xz)) m.e = vec3(1.,0.,0.);
  } else {
  if (h.mid == 2) {
    m.e = vec3(0.);
    vec3 wn=floor(h.p);
    if (mod(wn,vec3(3.)) == vec3(0.)) {
      m.e = 10. * (vec3(.6) + .2 * (noise(wn.yz).wzy+noise(wn.xx).xwz)) * step(1.1, noise(wn.xy).x + noise(wn.zz).z);
    }
    m.Cd = vec3(.2) + .2*h.c.S.x+.1*h.c.S.ywz;
  } else {
  if (h.mid == 3) {
    m.e = vec3(0.);
    m.Cd = mix(mg(h.p.xz), vec3(.1), step(min(cp.x,cp.z), 2.));
    m.Cs = vec3(.3);
    //m.s = 0.; // well well, this is some nasty bug we have here, nvidia
  }}}
  return m;
}

vec3 brdf(hit_t h, mat_t m, vec3 wi) {
  vec3 wo=-h.i,wh=normalize(wi+wo);
  float ci=max(0.,dot(h.n,wi)),co=max(0.,dot(h.n,wo)),ch=dot(wi,wh);
  //return m.Cd * ci;
  return // diffuse
    .3875077 * m.Cd * (vec3(1.) - m.Cs) * (1. - pow(1.-.5*ci,5.)) * (1. - pow(1.-.5*co,5.))
    //+ // specular (borken)
    //.0015831 * (m.s + 4.) * pow(dot(h.n, wh), m.s) * (m.Cs + pow(1. - ch, 5.) * (vec3(1.) - m.Cs)) / max(1e-1, ch * max(ci, co))
    ;
}

vec3 ref(hit_t h){
  vec3 r=rand().zyx*2.-vec3(1.);
  return normalize(r*sign(dot(r,h.n)));
}

vec3 air(vec3 O, vec3 D) {
  return 30. * vec3(1.) * smoothstep(.999,.9999,dot(D,sundir))
    + pow(
    vec3(128., 218., 235.)/255.
    , vec3(2.2));

}

// DEBUG
#define MAKE_Q(T) T quantize(T a,T b,int n,T v){return floor(float(n)*(v-a)/(b-a));}
MAKE_Q(vec3)
MAKE_Q(float)

#define CHECK_ASSERT if(assertcolor_.a!=0.){gl_FragColor=assertcolor_*1e3;return;}

void main(){
  //vec2 res=vec2(1280.,720.);
  vec2 res=vec2(1920.,1080.);
  vec2 uv=gl_FragCoord.xy/res-vec2(.5);uv.x*=res.x/res.y;
  vec3 O=_cp,D=normalize(vec3(uv,2.))*_cm;
  O.y = max(O.y, h2(O.xz)+10.);

  // FIXME replace with lens/sampler model
  O+=noise(vec2(_t)*gl_FragCoord.xy).xyz*vec3(.3,.3,.3);

  vec3 color = vec3(0.), transm=vec3(1.);
#if 0
  hit_t t = trace_grid(O, D, FAR);
  CHECK_ASSERT
  gl_FragColor = vec4(quantize(vec3(0.),vec3(256.),8,vec3(t._gc, t._mc, t._mc+t._gc)), 1.);return;
  //gl_FragColor = vec4(vec3(quantize(0.,256.,8,t._gc)), 1.);return;
  //gl_FragColor = vec4(vec3(quantize(0.,256.,8,t._gc+t._mc)), 1.);return;
  //{
  //  float dh=t.p.y-H(t.p.xz);
  //  if (dh>EPS){gl_FragColor=vec4(1.,0.,0.,1.)*100.;return;}
  //  if (dh<0.){gl_FragColor=vec4(1.,0.,1.,1.)*100.;return;}
  //}

  if (t.mid > 0) {
    color = t.n*4.;
  } else {
    //color = h.l/FAR;//air(O, D);
    //color = air(O, D);
    color = vec3(1000.,0.,0.);
  }
#else
  float Lmax = FAR;
  for (int i=0;i<BOUNCES;++i){
    if (dot(transm,transm) < .001) break;
    hit_t h = trace_grid(O, D, Lmax);
  CHECK_ASSERT

    if (h.mid > 0) {
      mat_t m = material(h);
      vec3 c = m.e;
      O = h.p + h.n * EPS;
      // importance
      if (trace_grid(O, sundir, 100.).l >= 100.) c += brdf(h, m, sundir) * air(O, sundir);
      D = ref(h);
      color += transm * c;
      transm *= brdf(h, m, D);
      Lmax *= .6;
    } else {
      color += transm * air(O, D);
    }
  }
  //color = max(vec3(0.), color);
#endif
  gl_FragColor = vec4(color,1.);
}
