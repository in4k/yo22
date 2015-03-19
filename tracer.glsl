const float EPS = .01, FAR = 5000., SKY = FAR, GS = 32.;

vec3 sundir = normalize(_s);//vec3(.1,.036,.037));

float H(vec2 p){vec4 c=TT(p);return c.z;}

vec3 TN(vec2 p){
  vec2 e=vec2(1.,.0);
  vec3 dx=vec3(2.*e.x,H(p+e.xy)-H(p-e.xy),0.);
  vec3 dz=vec3(0.,H(p+e.yx)-H(p-e.yx),2.*e.x);
  //return normalize(cross(dz,dx));//+1.5*n4(p*134.241));
  return normalize(cross(dz,dx)+3.5*(n4(p*7.31).yxz-vec3(.5)));
}

float maxv(vec3 v){return max(v.x,max(v.y,v.z));}
float dbox(vec3 p,vec3 sz){
  return maxv(abs(p)-sz);
}

float TD(vec3 p,vec3 D) {
  float H = H(p.xz);
  return (p.y - H);
}

//#define MAKE_BINSEARCH(NAME_,DISTF_,STEPS_,ARG_) \
//float NAME_(vec3 O, vec3 D, float e, float lo, float li){\
//  for(int i = 0; i < STEPS_; ++i) {\
//    float lm = (lo + li) * .5;\
//    if (DISTF_(O + D * lm, ARG_) < e) li = lm; else lo = lm;\
//  }\
//  return lo;\
//}
//MAKE_BINSEARCH(TB, TD, 8, D)
float TB(vec3 O, vec3 D, float e, float lo, float li) {
  for(int i = 0; i < 8; ++i) {
    float lm = (lo + li) * .5;
    if (TD(O + D * lm, D) < e) li = lm; else lo = lm;
  }
  return lo;
}

struct cell_t {
  vec2 i;
  vec3 c,m,M,BB;
  vec4 S;
  float B,H,h,p,Bh,f;
};

cell_t CI(vec3 p) {
  cell_t c;
  c.i = floor(p.xz/GS);
  c.m = c.i.xyy*GS;
  vec4 P = PP(c.m.xz); c.p = P.x;
  c.B = P.y; c.H = P.w; c.h = P.z; c.Bh = c.B - c.H;
  c.m.y = c.h;
  c.M = c.m + vec3(GS, c.B, GS);
  c.c = c.m + vec3(GS, c.H-c.h, GS) * .5;
  c.S = n4(c.i+vec2(23.,17.)).wyzx;
  c.BB = vec3(GS*.5-4., c.Bh, GS*.5-4.);
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

float BD(vec3 p, cell_t c) {
  p -= c.c; // p is now in cell-space
  float d = (c.Bh < 10.) ? 1e6 : bd1(p, c);
  return max(d, dbox(p, c.BB));
}

vec3 block_normal(vec3 p, cell_t c){
  vec2 e=vec2(.01,.0);
  return normalize(vec3(
    BD(p+e.xyy,c)-BD(p-e.xyy,c),
    BD(p+e.yxy,c)-BD(p-e.yxy,c),
    BD(p+e.yyx,c)-BD(p-e.yyx,c)
    ));
}

struct hit_t {
  int mid;
  float l;
  vec3 p, i, n;
  cell_t c;
  float _gc, _mc; // DEBUG
};

hit_t BH(hit_t h) {
  h.mid = 2;
  h.n = block_normal(h.p, h.c);
  return h;
}

hit_t TH(hit_t h) {
  h.mid = (h.c.Bh > 10.) ? 3 : 1;
  h.n = TN(h.p.xz);
  return h;
}

//vec4 assertcolor_=vec4(0.);
//#define ASSERT(cond,r,g,b) if(!(cond)){assertcolor_=vec4(r,g,b,1.);}

//float minpos(float a, float b) {
//  return (a<0.)?b:((b<0.)?a:min(a,b));
//}

float xp(float o, float d, float m, float M) {
  if (abs(d) < 1e-6) return 1e6;
  return ((d>0.) ? (M-o) : (m-o)) / d;
}

hit_t RT(vec3 O, vec3 D, float Lmax) {
  hit_t h;
  h.mid = -1;
  h.l = 0.;
  h.i = D;
  h._gc = h._mc = 0.; // DEBUG
  float dl = 0.;
  float pl = h.l;
  for (int i = 0; i < 256; ++i) {
    if (h.l > Lmax) break;
    float de = EPS * h.l * 2e-2;
    h.p = O + D * h.l;
    if (dl < 0. || h.p.y > h.c.B) {
      h._gc += 1.;
      h.l += max(0.,dl) + EPS;
      h.p = O + D * h.l;
      h.c = CI(h.p);
      // TODO skipsize dependent on D.y
      //float dx = xp(h.p.x, D.x, h.c.m.x - GS, h.c.m.x + GS);
      //float dz = xp(h.p.z, D.z, h.c.m.z - GS, h.c.m.z + GS);
      float dx = xp(h.p.x, D.x, h.c.m.x, h.c.m.x + GS);
      float dz = xp(h.p.z, D.z, h.c.m.z, h.c.m.z + GS);
      //ASSERT(dz >= 0.,0.,1.,1.)
      //ASSERT(dx >= 0.,0.,.5,1.)
      //ASSERT(h.p.z >= cp.y,1.,1.,0.)
      //ASSERT(h.p.z <= cp.y+GS,1.,0.,1.)
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
      float d = BD(h.p, h.c);
      if (d < de) return BH(h);
      if (h.p.y < h.c.H)
      {
        d = min(d, TD(h.p, D));
        if (d < de) {
          h.l = TB(O, D, de, pl, h.l);
          h.p = O + D * h.l;
          return TH(h);
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
  return vec3(.2,.6,.23)*.3 + .1 * (n4(p*.2).xxx-vec3(.5));
}

mat_t M(hit_t h) {
  mat_t m;
  m.e = vec3(1., 0., 1.);
  m.Cd = vec3(0.);
  m.Cs = vec3(0.);
  m.s = 0.;
  vec3 cp = abs(abs(h.p - h.c.c) - vec3(GS,0.,GS)*.5);

  if(h.mid == 1) {
    float f = step(h.p.y+h.c.S.w*50.+20.*n4(h.p.xz).x,180.);
    m.e = vec3(0.);
    m.Cd = mg(h.p.xz) + vec3(.1,.1,.02)*h.c.S.xyz*f;
    m.Cd = mix(m.Cd, vec3(.11, .07, .03), .4*f*step(min(cp.x,cp.z), 1.));
    //if (h.p.y < H(h.p.xz)) m.e = vec3(1.,0.,0.);
  } else {
  if (h.mid == 2) {
    m.e = vec3(0.);
    vec3 wn=floor(h.p);
    if (mod(wn,vec3(3.)) == vec3(0.)) {
      m.e = 10. * (vec3(.6) + .2 * (n4(wn.yz).wzy+n4(wn.xx).xwz)) * step(1.1, n4(wn.xy).x + n4(wn.zz).z);
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
//#define MAKE_Q(T) T quantize(T a,T b,int n,T v){return floor(float(n)*(v-a)/(b-a));}
//MAKE_Q(vec3)
//MAKE_Q(float)
//#define CHECK_ASSERT if(assertcolor_.a!=0.){gl_FragColor=assertcolor_*1e3;return;}

void main(){
  //vec2 res=vec2(1280.,720.);
  vec2 res=vec2(1920.,1080.);
  vec2 uv=gl_FragCoord.xy/res-vec2(.5);uv.x*=res.x/res.y;
  vec3 O=_cp,D=normalize(vec3(uv,2.))*_cm;
  //O.y = max(O.y, h2(O.xz)+10.);

  // FIXME replace with lens/sampler model
  O+=n4(vec2(_t)*gl_FragCoord.xy).xyz*vec3(.3,.3,.3);

  vec3 color = vec3(0.), ct=vec3(1.);
//#if 0
  //hit_t t = trace_grid(O, D, FAR);
  //CHECK_ASSERT
  //gl_FragColor = vec4(quantize(vec3(0.),vec3(256.),8,vec3(t._gc, t._mc, t._mc+t._gc)), 1.);return;
  //gl_FragColor = vec4(vec3(quantize(0.,256.,8,t._gc)), 1.);return;
  //gl_FragColor = vec4(vec3(quantize(0.,256.,8,t._gc+t._mc)), 1.);return;
  //{
  //  float dh=t.p.y-H(t.p.xz);
  //  if (dh>EPS){gl_FragColor=vec4(1.,0.,0.,1.)*100.;return;}
  //  if (dh<0.){gl_FragColor=vec4(1.,0.,1.,1.)*100.;return;}
  //}

//  if (t.mid > 0) {
//    color = t.n*4.;
//  } else {
    //color = h.l/FAR;//air(O, D);
    //color = air(O, D);
//    color = vec3(1000.,0.,0.);
//  }
//#else
  float Lmax = FAR;
  for (int i = 0; i < 3; ++i){
    if (dot(ct,ct) < .001) break;
    hit_t h = RT(O, D, Lmax);
  //CHECK_ASSERT

    if (h.mid > 0) {
      mat_t m = M(h);
      vec3 c = m.e;
      O = h.p + h.n * EPS;
      // importance
      if (RT(O, sundir, 100.).l >= 100.) c += brdf(h, m, sundir) * air(O, sundir);
      D = ref(h);
      color += ct * c;
      ct *= brdf(h, m, D);
      Lmax *= .6;
    } else {
      color += ct * air(O, D);
    }
  }
  //color = max(vec3(0.), color);
//#endif
  gl_FragColor = vec4(color,1.);
}
