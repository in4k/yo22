const float EPS = .01, FAR = 5000., SKY = FAR, GS = 32.;

float Mv(vec3 v){return max(v.x,max(v.y,v.z));}
float bx(vec3 p,vec3 sz){
  return Mv(abs(p)-sz);
}


float H(vec2 p){vec4 c=TT(p);return c.z
  // grass
  +.5*n4(p).z*step(.5,n4(p*.17).w)
  // forests
  +4.*n4(p*.3).z*step(.6,n4(p*.02).w)
;}
float TD(vec3 p,vec3 D) {
  return p.y - H(p.xz);
}
vec3 TN(vec2 p){
  vec2 e=vec2(.01,.0);
  vec3 dx=vec3(2.*e.x,H(p+e.xy)-H(p-e.xy),0.);
  vec3 dz=vec3(0.,H(p+e.yx)-H(p-e.yx),2.*e.x);
  //return normalize(cross(dz,dx));//+1.5*n4(p*134.241));
  return normalize(cross(dz,dx)+.001*(n4(p*7.31).yxz-vec3(.5)));
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
  vec3 c, m, M, BB;
  vec4 S;
  float B, H, h, p, Bh, f;
};
cell_t c;

int hm;
float hl;
vec3 hp, hi, hn;

// colors: emission, diffuse, specular
vec3 Ce, Cd, Cs;

void CI(vec3 p) {
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
}

//float bd1(vec3 p) {
//  vec2 e=vec2(1.,0.);
//  float d = 1e6;
//  vec3 s = c.BB * vec3(.4+.3*c.S.z,1.,.4+.3*c.S.y);
//  for (int i = 0; i < 4; ++i)
//  {
//    vec3 pp = p-(c.BB-s)*(2.*vec3(c.S.w,0.,c.S.z)-vec3(1.));
//    d = min(d, bx(pp, s));
//    s += vec3(.2+c.S.x, 0., .2+c.S.y) * .2 * c.BB;
//    s.y *= .7+.1*c.S.z;
//  }
//  return d;
//}

float BD(vec3 p) {
  p -= c.c; // p is now in cell-space
  if (c.Bh < 10.) return 1e6;
  //float d = (c.Bh < 10.) ? 1e6 : bd1(p);
  vec2 e=vec2(1.,0.);
  float d = 1e6;
  vec3 s = c.BB * vec3(.4+.3*c.S.z,1.,.4+.3*c.S.y);
  for (int i = 0; i < 4; ++i)
  {
    vec3 pp = p-(c.BB-s)*(2.*vec3(c.S.w,0.,c.S.z)-vec3(1.));
    d = min(d, bx(pp, s));
    s += vec3(.2+c.S.x, 0., .2+c.S.y) * .2 * c.BB;
    s.y *= .7+.1*c.S.z;
  }
  return max(d, bx(p, c.BB));
}

vec3 BN(vec3 p){
  vec2 e=vec2(.01,.0);
  return normalize(vec3(
    BD(p+e.xyy)-BD(p-e.xyy),
    BD(p+e.yxy)-BD(p-e.yxy),
    BD(p+e.yyx)-BD(p-e.yyx)
    ));
}

// grass color
vec3 mg(vec2 p) {
  return vec3(.2,.6,.23)*.3 + .2 * (n4(p*.2).xxx-vec3(.5));
}

void BH() {
  hm = 1;
  vec3 wn = floor(hp);
  if (mod(wn,vec3(3.)) == vec3(0.))
    Ce = 10. * (vec3(.6) + .2 * (n4(wn.yz).wzy+n4(wn.xx).xwz)) * step(1.1, n4(wn.xy).x + n4(wn.zz).z);
  else
    Ce = vec3(0.);
  Cd = vec3(.2) + .2*c.S.x+.1*c.S.ywz;
  //Cs = vec3(.0);
  hn = BN(hp);
}

void TH() {
  hm = 1;
  vec3 cp = abs(abs(hp - c.c) - vec3(GS,0.,GS)*.5);
  if(c.Bh < 10.) {
    float f = step(hp.y+c.S.w*50.+20.*n4(hp.xz).x,180.);
    //Ce = vec3(0.);
    Cd = mg(hp.xz) + vec3(.1,.1,.02)*c.S.xyz*f;
    Cd = mix(Cd, vec3(.11, .07, .03), .4*f*step(min(cp.x,cp.z), 1.));
  } else {
    //Ce = vec3(0.);
    Cd = mix(mg(hp.xz), vec3(.1), step(min(cp.x,cp.z), 2.));
    //Cs = vec3(.3);
    //Ce = vec3(1., 1., 0.);
    Ce =
      (vec3(1.,1.,.4) * step(.92,n4(cp.xz).z)
      +vec3(1.,0.,0.) * step(.92,n4(cp.xz).y))
      * step(min(cp.x,cp.z), 2.) * 3.;

    //m.s = 0.; // well well, this is some nasty bug we have here, nvidia
  }
  hn = TN(hp.xz);
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

float RT(vec3 O, vec3 D, float L) {
  hm = -1;
  hl = 0.;
  hi = D;
  float dl = 0., pl = hl, de, dx, dy, dz, d;
  for (int i = 0; i < 256; ++i) {
    if (hl > L) break;
    float de = EPS * hl * 2e-2;
    hp = O + D * hl;
    if (dl < 0. || hp.y > c.B) {
      hl += max(0.,dl) + EPS;
      hp = O + D * hl;
      CI(hp);
      // TODO skipsize dependent on D.y
      //float dx = xp(hp.x, D.x, c.m.x - GS, hc.m.x + GS);
      //float dz = xp(hp.z, D.z, c.m.z - GS, hc.m.z + GS);
      dx = xp(hp.x, D.x, c.m.x, c.m.x + GS);
      dz = xp(hp.z, D.z, c.m.z, c.m.z + GS);
      //ASSERT(dz >= 0.,0.,1.,1.)
      //ASSERT(dx >= 0.,0.,.5,1.)
      //ASSERT(hp.z >= cp.y,1.,1.,0.)
      //ASSERT(hp.z <= cp.y+GS,1.,0.,1.)
      dy =
        (hp.y > c.B) ?
          xp(hp.y, D.y, c.B, SKY) :
        ((hp.y > c.H) ?// && (P.y-P.w) > 0.) ?
          xp(hp.y, D.y, c.H, c.B) :
          1e6);//xp(hp.y, D.y, hc.h, hc.H));
      //ASSERT(dy >= 0.,.5,.5,1.)
      dl = min(dx,min(dy,dz));
      //ASSERT(dl >= 0.,.5,.5,.5)
      dl += EPS;
      //if (dl < 0.) {c=vec3(1.,float(i),hl);break;}
    } else
    //h_mc -= 1.;}
    {
      d = BD(hp);
      if (d < de) {BH(); break;}
      if (hp.y < c.H)
      {
        d = min(d, TD(hp, D));
        if (d < de) {
          hl = TB(O, D, de, pl, hl);
          hp = O + D * hl;
          TH(); break;
        }
      }
      d = min(d,dl+EPS);
      dl -= d;
      pl = hl;
      hl += d;
    }
  }
  return hl;
}

vec3 brdf(vec3 wi) {
  vec3 wo=-hi,wh=normalize(wi+wo);
  float ci=max(0.,dot(hn,wi)),co=max(0.,dot(hn,wo)),ch=dot(wi,wh);
  //return mCd * ci;
  return // diffuse
    .3875077 * Cd * (vec3(1.) - Cs) * (1. - pow(1.-.5*ci,5.)) * (1. - pow(1.-.5*co,5.))
    //+ // specular (borken)
    //.0015831 * (m.s + 4.) * pow(dot(hn, wh), m.s) * (mCs + pow(1. - ch, 5.) * (vec3(1.) - mCs)) / max(1e-1, ch * max(ci, co))
    ;
}

vec3 ref(){
  //vec3 r=rand().zyx*2.-vec3(1.);
  vec3 r = n4(fc+floor(hp.zy*313.)).zyx + n4(floor(hp.zx*171.)-fc).wxy - vec3(1.);
  return normalize(r*sign(dot(r,hn)));
}

const float I=10.,g=.76,g2=g*g,R0=6360e3,R1=6420e3;
const vec3 C=vec3(0.,-R0,0.),Km=vec3(21e-6),Kr=vec3(5.8,13.5,33.1)*1e-6;

float ER(vec3 o, vec3 d, float r) {
  o -= C;
  float
    b = dot(o, d),
    c = dot(o, o) - r * r,
    e = b * b - c,
    p, q;
  if (e < 0.) return -1.;
  e = sqrt(e);
  p = -b - e,
  q = -b + e;
  return (p >= 0.) ? p : q;
}

vec2 SD(vec3 O) {
  float h = R0 - length(O - C);
  return vec2(exp(h / 120.), exp(h / 8000.));
}

vec2 Id(vec3 O, vec3 D, float L) {
  float dx, bdx = L / 16., l = 0.;
  vec3 p = O;
  vec2 mr = vec2(0.);
  for (int i = 0; i < 16; ++i) {
    dx = bdx * (1. + n4(p.xz * 2e3).w);
    p = O + D * (l + dx);
    mr += SD(p) * dx;
    l += bdx;
  }
  return mr;
}

vec3 S(vec3 O, vec3 D, float L, vec3 cc) {
	float c = max(0., dot(D, _s)),
  Pr = .0596831 * (1. + c * c),
  Pm = .1193662 * (1. - g2) * (1. + c * c) / ((2. + g2) * pow(1. + g2 - 2.*g*c, 1.5)),
  dm = 0.,
  dr = 0.,
  ll = 0.,
  bdx = L / 32., sL;

  vec3 im = vec3(0.), ir = vec3(0.), p = O, Is;
  vec2 mr;
  for (int i = 0; i < 32; ++i) {
    float dx = bdx * (1. + n4(p.xz * 2e3).w);
    p = O + D * (ll + dx);
    ll += bdx;
    mr = SD(p) * dx;
    dm += mr.x;
    dr += mr.y;

    // if shadowed by terrain
    //float Lt = TT(p, _s);
    //if (Lt < Far) continue;

    sL = ER(p,_s,R1);

    vec2 dd = Id(p, _s, sL);

    Is = exp(-(
        Km * (dd.x + dm) +
        Kr * (dd.y + dr)
      ));
    im += mr.x * Is;
    ir += mr.y * Is;
  }

  return max(vec3(0.),I * (
      Km * Pm * im +
      Kr * Pr * ir)
    + cc * exp(-(Km * dm + Kr * dr)));
}

// DEBUG
//#define MAKE_Q(T) T quantize(T a,T b,int n,T v){return floor(float(n)*(v-a)/(b-a));}
//MAKE_Q(vec3)
//MAKE_Q(float)
//#define CHECK_ASSERT if(assertcolor_.a!=0.){gl_FragColor=assertcolor_*1e3;return;}

void main(){
  vec2 res=vec2(1280.,720.),
  //vec2 res=vec2(1920.,1080.);
  uv=V;uv.x*=res.x/res.y;
  vec3 C = vec3(0.), ct = vec3(1.),O=_cp+n4(vec2(_t)*fc).xyz*.3,D=normalize(vec3(uv,2.))*_cm,pO,loc;
  //O.y = max(O.y, h2(O.xz)+10.);

  // FIXME replace with lens/sampler model
  //O;

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

//  if (t.m > 0) {
//    color = t.n*4.;
//  } else {
    //color = hl/FAR;//air(O, D);
    //color = air(O, D);
//    color = vec3(1000.,0.,0.);
//  }
//#else
  float f = FAR;
  pO = O;
  for (int i = 0; i < 3; ++i){
    if (dot(ct,ct) < .001) break;
    Ce = Cd = Cs = vec3(0.);
    RT(O, D, f);
  //CHECK_ASSERT

    if (hm > 0) {
      loc = Ce;
      O = hp + hn * EPS;
      // importance
      vec3 cd = Cd, cs = Cs, chi = hi, chn = hn, chp = hp;
      cell_t cc = c;
      float chl = hl,
      shl = RT(O, _s, 100.);
      Cd = cd; Cs = cs; hn = chn; hp = chp; hi = chi; c = cc; hl = chl;
      if (shl > 100.) loc += brdf(_s) * S(O, _s, ER(O, _s, R1), vec3(I));
      C += ct * S(pO, D, hl, loc);
      D = ref();
      ct *= brdf(D);
      f *= .6;
    } else {
      C += ct * S(O, D, ER(O, D, R1), vec3(0.));
      break;
    }
  }
//#endif
  gl_FragColor = vec4(C,1.);
}
