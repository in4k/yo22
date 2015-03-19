static const char *fragment_shaders[Prog_COUNT] = {
"uniform float _t,_p;\n"
"uniform vec2 _r;\n"
"uniform sampler2D _N,_T,_P,_F;\n"
"varying vec2 V;\n"
"\n"
"vec4 n4(vec2 v){return texture2D(_N,(v+vec2(.5))/256.,-20.);}\n"
"float n(vec2 v){return n4(v).w;}\n"
"float vn(vec2 v,vec2 m){\n"
"vec2 e=vec2(1.,0.),V=floor(v);v=fract(v);v*=v*(3.-2.*v);\n"
"return mix(mix(n(mod(V+e.yy,m)),n(mod(V+e.xy,m)),v.x),mix(n(mod(V+e.yx,m)),n(mod(V+e.xx,m)),v.x),v.y);\n"
"}\n"
"float fbm(vec2 v,float s){\n"
"float r=0.,k=.5;\n"
"for(int i=0;i<12;++i,k*=.5,s*=2.)r+=k*vn(v*s,vec2(s));\n"
"return r;\n"
"}\n"
"void main(){\n"
"vec2 uv=floor(gl_FragCoord.xy)/4096.;\n"
"float height=fbm(uv,5.);\n"
"height=pow(height,3.)*1200.;\n"
"gl_FragColor=vec4(0.,0.,height,0.);\n"
"}\n"
, 
"uniform float _t,_p;\n"
"uniform sampler2D _N,_T;\n"
"varying vec2 V;\n"
"vec4 n4(vec2 v){return texture2D(_N,(v+vec2(.5))/256.,-20.);}\n"
"float n(vec2 v){return n4(v).w;}\n"
"float vn(vec2 v,vec2 m){\n"
"vec2 e=vec2(1.,0.),V=floor(v);v=fract(v);v*=v*(3.-2.*v);\n"
"return mix(mix(n(mod(v+e.yy,m)),n(mod(v+e.xy,m)),v.x),mix(n(mod(v+e.yx,m)),n(mod(v+e.xx,m)),v.x),v.y);\n"
"}\n"
"float fbm(vec2 v,float s){\n"
"float r=0.,k=.5;\n"
"for(int i=0;i<12;++i,k*=.5,s*=2.)r+=k*vn(v*s,vec2(s));\n"
"return r;\n"
"}\n"
"float\n"
"Kr=0.27,\n"
"Ke=0.005,\n"
"Ks=0.8;\n"
"vec2 uv=floor(gl_FragCoord.xy);\n"
"float transfer(float H,float w,vec2 ofs){\n"
"vec4 c=texture2D(_T,(uv+ofs+vec2(.5))/4096.,-20.);\n"
"float h=c.w+c.z,d=H-h;\n"
"if(d>0.)return-min(w,d);\n"
"else return min(c.w,-d);\n"
"}\n"
"void main(){\n"
"vec4 c=texture2D(_T,(uv+vec2(.5))/4096.,-20.);\n"
"\n"
"\n"
"vec4 c0=c;\n"
"float h=c.w+c.z,dw=0.;\n"
"vec2 e=vec2(1.,0.);\n"
"dw+=transfer(h,c.w,e.xy);\n"
"dw+=transfer(h,c.w,e.yx);\n"
"dw+=transfer(h,c.w,-e.xy);\n"
"dw+=transfer(h,c.w,-e.yx);\n"
"#if 0\n"
"float kk=1./(4.*sqrt(2.));\n"
"dw+=transfer(h,c.w,vec2(1.,1.))*kk;\n"
"dw+=transfer(h,c.w,vec2(1.,-1.))*kk;\n"
"dw+=transfer(h,c.w,vec2(-1.,1.))*kk;\n"
"dw+=transfer(h,c.w,-vec2(1.,1.))*kk;\n"
"c.w+=dw/4.;\n"
"#else\n"
"c.w+=dw/4.;\n"
"#endif\n"
"\n"
"\n"
"dw=c.w*Ke;\n"
"c.z+=dw*Ks;\n"
"c.w-=dw;\n"
"\n"
"float sa=sin(_t*373.4111),ca=cos(_t*373.4111);\n"
"mat2 m=mat2(sa,ca,ca,-sa);\n"
"dw=Kr*n4((uv+vec2(_t*123.24,317.1141*_t))*m).x*(0==mod(floor(_t/3.),12.)?.1:.0);\n"
"\n"
"c.z-=dw*Ks;\n"
"c.w+=dw;\n"
"\n"
"gl_FragColor=vec4(\n"
"1000.*(c.w-c0.w),\n"
"-1000.*(c.z-c0.z),\n"
"c.z,\n"
"c.w);\n"
"}\n"
, 
"uniform float _t,_p;\n"
"uniform vec2 _r;\n"
"uniform sampler2D _N,_T,_P,_F;\n"
"varying vec2 V;\n"
"\n"
"vec4 n4(vec2 p){return texture2D(_N,(p+vec2(.5))/1024.,-20.);}\n"
"\n"
"void main(){\n"
"vec2 fc=floor(gl_FragCoord.xy);\n"
"vec2 p=fc*32.+vec2(.5);\n"
"float h=1e6,H=0.;\n"
"for(int x=0;x<32;++x){\n"
"for(int y=0;y<32;++y){\n"
"vec2 P=p+vec2(float(x),float(y));\n"
"vec4 s=texture2D(_T,P/4096.,-20.);\n"
"H=max(H,s.z);h=min(h,s.z);\n"
"}\n"
"}\n"
"float d=H-h;\n"
"vec4 r=n4(fc);\n"
"float pop=pow(smoothstep(25.,0.,d)*smoothstep(230.,110.,H),14.);\n"
"pop*=.2+.8*r.z;\n"
"pop+=step(0.,pop)*.08;\n"
"gl_FragColor=vec4(r.y,max(0.,pop*90.)+H,h,H);\n"
"\n"
"}\n"
, 
"uniform float _t,_p;\n"
"uniform vec2 _r;\n"
"uniform vec3 _s,_cp;\n"
"uniform mat3 _cm;\n"
"uniform sampler2D _N,_T,_P,_F;\n"
"varying vec2 V;\n"
"\n"
"vec4 noise(vec2 p){return texture2D(_N,(p+.5)/1024.,-20.);}\n"
"vec4 terrain(vec2 p_m){return texture2D(_T,(p_m+vec2(.5,.5))/4096.,-20.);}\n"
"vec4 plan(vec2 p_m){return texture2D(_P,(floor(p_m/32.)+vec2(.5))/128.,-20.);}\n"
"vec2 fc=floor(gl_FragCoord.xy);\n"
"\n"
"\n"
"int rand_state_=int(fc.y*1280.+fc.x)+int(noise(vec2(_t,0.)).x*256.+noise(vec2(0.,_t)).z)*256;\n"
"vec4 rand(){rand_state_=int(mod(float(rand_state_+1),1024.*1024.));return noise(vec2(float(rand_state_),floor(float(rand_state_)/1024.)));}\n"
"vec4 n4(vec2 v){return texture2D(_N,(v+vec2(.5))/256.,-20.);}\n"
"float n(vec2 v){return n4(v).w;}\n"
"float vn(vec2 v,vec2 m){\n"
"vec2 e=vec2(1.,0.),V=floor(v);v=fract(v);v*=v*(3.-2.*v);\n"
"return mix(mix(n(mod(v+e.yy,m)),n(mod(v+e.xy,m)),v.x),mix(n(mod(v+e.yx,m)),n(mod(v+e.xx,m)),v.x),v.y);\n"
"}\n"
"float fbm(vec2 v,float s){\n"
"float r=0.,k=.5;\n"
"for(int i=0;i<4;++i,k*=.5,s*=2.)r+=k*vn(v*s,vec2(s));\n"
"return r;\n"
"}\n"
"\n"
"#define STEPS 256\n"
"#define EPS.01\n"
"#define FAR 5000.\n"
"#define BOUNCES 3\n"
"#define SKY FAR\n"
"#define GRIDSIZE 32.\n"
"\n"
"vec3 sundir=normalize(_s);\n"
"\n"
"float H(vec2 p){\n"
"vec4 c=terrain(p);\n"
"return c.z;}\n"
"float h2(vec2 p){\n"
"vec4 c=terrain(p);\n"
"return c.z;\n"
"}\n"
"\n"
"vec3 terrain_normal(vec2 p){\n"
"vec2 e=vec2(1.,.0);\n"
"vec3 dx=vec3(2.*e.x,H(p+e.xy)-H(p-e.xy),0.);\n"
"vec3 dz=vec3(0.,H(p+e.yx)-H(p-e.yx),2.*e.x);\n"
"\n"
"return normalize(cross(dz,dx)+3.5*(noise(p*7.31).yxz-vec3(.5)));\n"
"}\n"
"\n"
"vec3 terrain_albedo(vec3 p){\n"
"#if 1\n"
"\n"
"vec4 P=plan(p.xz);\n"
"\n"
"float d=P.w-P.z;\n"
"return vec3(.2,.2+1.-min(1.,floor(d/5.)/2.),.2);\n"
"\n"
"\n"
"#elif 1\n"
"return vec3(.2);\n"
"#elif 1\n"
"float m10=mod(floor(p.x/10.)+floor(p.z/10.),2.);\n"
"float m100=mod(floor(p.x/100.)+floor(p.z/100.),2.);\n"
"float m1000=mod(floor(p.x/1000.)+floor(p.z/1000.),2.);\n"
"return vec3(m10,m100,m1000)*(.3+.7*mod(floor(p.x)+floor(p.z),2.));\n"
"#elif 1\n"
"float XZ=mod(floor(p.x/1.)+floor(p.z/1.),2.);\n"
"float Y=mod(floor(p.y),2.);\n"
"return vec3(.2+.5*Y+.3*XZ);\n"
"#else\n"
"return texture2D(_T,p.xz/4096.,-20.).xyw*10.;\n"
"#endif\n"
"}\n"
"\n"
"float maxv(vec3 v){return max(v.x,max(v.y,v.z));}\n"
"float dbox(vec3 p,vec3 sz){\n"
"return maxv(abs(p)-sz);\n"
"}\n"
"\n"
"float dist_terrain(vec3 p,vec3 D){\n"
"float H=H(p.xz);\n"
"return(p.y-H);\n"
"}\n"
"\n"
"\n"
"#define MAKE_BINSEARCH(NAME_,DISTF_,STEPS_,ARG_)\\\n"
"float NAME_(vec3 O,vec3 D,float e,float lo,float li){\\\n"
"for(int i=0;i < STEPS_;++i){\\\n"
"float lm=(lo+li)*.5;\\\n"
"if(DISTF_(O+D*lm,ARG_)< e)li=lm;else lo=lm;\\\n"
"}\\\n"
"return lo;\\\n"
"}\n"
"MAKE_BINSEARCH(bin_terrain,dist_terrain,8,D)\n"
"\n"
"struct cell_t{\n"
"vec2 i;\n"
"vec3 c,m,M,BB;\n"
"vec4 S;\n"
"float B,H,h,p,Bh,f;\n"
"};\n"
"\n"
"cell_t cell_identify(vec3 p){\n"
"cell_t c;\n"
"c.i=floor(p.xz/GRIDSIZE);\n"
"c.m=c.i.xyy*GRIDSIZE;\n"
"vec4 P=plan(c.m.xz);c.p=P.x;\n"
"c.B=P.y;c.H=P.w;c.h=P.z;c.Bh=c.B-c.H;\n"
"c.m.y=c.h;\n"
"c.M=c.m+vec3(GRIDSIZE,c.B,GRIDSIZE);\n"
"c.c=c.m+vec3(GRIDSIZE,c.H-c.h,GRIDSIZE)*.5;\n"
"c.S=noise(c.i+vec2(23.,17.)).wyzx;\n"
"c.BB=vec3(GRIDSIZE*.5-4.,c.Bh,GRIDSIZE*.5-4.);\n"
"c.f=9.+.2*c.S.y;\n"
"return c;\n"
"}\n"
"\n"
"float bd1(vec3 p,cell_t c){\n"
"vec2 e=vec2(1.,0.);\n"
"float d=1e6;\n"
"vec3 s=c.BB*vec3(.2+.3*c.S.z,1.,.2+.3*c.S.y);\n"
"for(int i=0;i < 4;++i)\n"
"{\n"
"vec3 pp=p-(c.BB-s)*(2.*vec3(c.S.w,0.,c.S.z)-vec3(1.));\n"
"d=min(d,dbox(pp,s));\n"
"s+=vec3(.2+c.S.x,0.,.2+c.S.y)*.2*c.BB;\n"
"s.y*=.7+.1*c.S.z;\n"
"}\n"
"return d;\n"
"}\n"
"\n"
"float dist_block(vec3 p,cell_t c){\n"
"p-=c.c;\n"
"float d=(c.Bh < 10.)? 1e6 : bd1(p,c);\n"
"return max(d,dbox(p,c.BB));\n"
"}\n"
"\n"
"vec3 block_normal(vec3 p,cell_t c){\n"
"vec2 e=vec2(.01,.0);\n"
"return normalize(vec3(\n"
"dist_block(p+e.xyy,c)-dist_block(p-e.xyy,c),\n"
"dist_block(p+e.yxy,c)-dist_block(p-e.yxy,c),\n"
"dist_block(p+e.yyx,c)-dist_block(p-e.yyx,c)\n"
"));\n"
"}\n"
"\n"
"struct hit_t{\n"
"int mid;\n"
"float l;\n"
"vec3 p,i,n;\n"
"cell_t c;\n"
"float _gc,_mc;\n"
"};\n"
"\n"
"hit_t hit_block(hit_t h){\n"
"h.mid=2;\n"
"h.n=block_normal(h.p,h.c);\n"
"return h;\n"
"}\n"
"\n"
"hit_t hit_terrain(hit_t h){\n"
"h.mid=(h.c.Bh > 10.)? 3 : 1;\n"
"h.n=terrain_normal(h.p.xz);\n"
"return h;\n"
"}\n"
"\n"
"vec4 assertcolor_=vec4(0.);\n"
"#define ASSERT(cond,r,g,b)if(!(cond)){assertcolor_=vec4(r,g,b,1.);}\n"
"\n"
"float minpos(float a,float b){\n"
"return(a<0.)?b:((b<0.)?a:min(a,b));\n"
"}\n"
"\n"
"float xp(float o,float d,float m,float M){\n"
"if(abs(d)< 1e-6)return 1e6;\n"
"return((d>0.)?(M-o):(m-o))/ d;\n"
"}\n"
"\n"
"hit_t trace_grid(vec3 O,vec3 D,float Lmax){\n"
"hit_t h;\n"
"h.mid=-1;\n"
"h.l=0.;\n"
"h.i=D;\n"
"h._gc=h._mc=0.;\n"
"float dl=0.;\n"
"float pl=h.l;\n"
"for(int i=0;i < STEPS;++i){\n"
"if(h.l > Lmax)break;\n"
"float de=EPS*h.l*2e-2;\n"
"h.p=O+D*h.l;\n"
"if(dl < 0.|| h.p.y > h.c.B){\n"
"h._gc+=1.;\n"
"h.l+=max(0.,dl)+EPS;\n"
"h.p=O+D*h.l;\n"
"h.c=cell_identify(h.p);\n"
"\n"
"\n"
"\n"
"float dx=xp(h.p.x,D.x,h.c.m.x,h.c.m.x+GRIDSIZE);\n"
"float dz=xp(h.p.z,D.z,h.c.m.z,h.c.m.z+GRIDSIZE);\n"
"\n"
"\n"
"\n"
"\n"
"float dy=\n"
"(h.p.y > h.c.B)?\n"
"xp(h.p.y,D.y,h.c.B,SKY):\n"
"((h.p.y > h.c.H)?\n"
"xp(h.p.y,D.y,h.c.H,h.c.B):\n"
"1e6);\n"
"\n"
"dl=min(dx,min(dy,dz));\n"
"\n"
"dl+=EPS;\n"
"\n"
"}else\n"
"\n"
"{\n"
"h._mc+=1.;\n"
"float d=dist_block(h.p,h.c);\n"
"if(d < de)return hit_block(h);\n"
"if(h.p.y < h.c.H)\n"
"{\n"
"d=min(d,dist_terrain(h.p,D));\n"
"if(d < de){\n"
"h.l=bin_terrain(O,D,de,pl,h.l);\n"
"h.p=O+D*h.l;\n"
"return hit_terrain(h);\n"
"}\n"
"}\n"
"d=min(d,dl+EPS);\n"
"dl-=d;\n"
"pl=h.l;\n"
"h.l+=d;\n"
"}\n"
"}\n"
"return h;\n"
"}\n"
"\n"
"struct mat_t{\n"
"vec3 e,Cd,Cs;\n"
"float s;\n"
"};\n"
"\n"
"vec3 mg(vec2 p){\n"
"return vec3(.2,.6,.23)*.3+.1*(noise(p*.2).xxx-vec3(.5));\n"
"}\n"
"\n"
"mat_t material(hit_t h){\n"
"mat_t m;\n"
"m.e=vec3(1.,0.,1.);\n"
"m.Cd=vec3(0.);\n"
"m.Cs=vec3(0.);\n"
"m.s=0.;\n"
"vec3 cp=abs(abs(h.p-h.c.c)-vec3(GRIDSIZE,0.,GRIDSIZE)*.5);\n"
"\n"
"if(h.mid==1){\n"
"float f=step(h.p.y+h.c.S.w*50.+20.*noise(h.p.xz).x,180.);\n"
"m.e=vec3(0.);\n"
"m.Cd=mg(h.p.xz)+vec3(.1,.1,.02)*h.c.S.xyz*f;\n"
"m.Cd=mix(m.Cd,vec3(.11,.07,.03),.4*f*step(min(cp.x,cp.z),1.));\n"
"\n"
"}else{\n"
"if(h.mid==2){\n"
"m.e=vec3(0.);\n"
"vec3 wn=floor(h.p);\n"
"if(mod(wn,vec3(3.))==vec3(0.)){\n"
"m.e=10.*(vec3(.6)+.2*(noise(wn.yz).wzy+noise(wn.xx).xwz))*step(1.1,noise(wn.xy).x+noise(wn.zz).z);\n"
"}\n"
"m.Cd=vec3(.2)+.2*h.c.S.x+.1*h.c.S.ywz;\n"
"}else{\n"
"if(h.mid==3){\n"
"m.e=vec3(0.);\n"
"m.Cd=mix(mg(h.p.xz),vec3(.1),step(min(cp.x,cp.z),2.));\n"
"m.Cs=vec3(.3);\n"
"\n"
"}}}\n"
"return m;\n"
"}\n"
"\n"
"vec3 brdf(hit_t h,mat_t m,vec3 wi){\n"
"vec3 wo=-h.i,wh=normalize(wi+wo);\n"
"float ci=max(0.,dot(h.n,wi)),co=max(0.,dot(h.n,wo)),ch=dot(wi,wh);\n"
"\n"
"return \n"
".3875077*m.Cd*(vec3(1.)-m.Cs)*(1.-pow(1.-.5*ci,5.))*(1.-pow(1.-.5*co,5.))\n"
"\n"
"\n"
";\n"
"}\n"
"\n"
"vec3 ref(hit_t h){\n"
"vec3 r=rand().zyx*2.-vec3(1.);\n"
"return normalize(r*sign(dot(r,h.n)));\n"
"}\n"
"\n"
"vec3 air(vec3 O,vec3 D){\n"
"return 30.*vec3(1.)*smoothstep(.999,.9999,dot(D,sundir))\n"
"+pow(\n"
"vec3(128.,218.,235.)/255.\n"
",vec3(2.2));\n"
"\n"
"}\n"
"\n"
"\n"
"#define MAKE_Q(T)T quantize(T a,T b,int n,T v){return floor(float(n)*(v-a)/(b-a));}\n"
"MAKE_Q(vec3)\n"
"MAKE_Q(float)\n"
"\n"
"#define CHECK_ASSERT if(assertcolor_.a!=0.){gl_FragColor=assertcolor_*1e3;return;}\n"
"\n"
"void main(){\n"
"\n"
"vec2 res=vec2(1920.,1080.);\n"
"vec2 uv=gl_FragCoord.xy/res-vec2(.5);uv.x*=res.x/res.y;\n"
"vec3 O=_cp,D=normalize(vec3(uv,2.))*_cm;\n"
"O.y=max(O.y,h2(O.xz)+10.);\n"
"\n"
"\n"
"O+=noise(vec2(_t)*gl_FragCoord.xy).xyz*vec3(.3,.3,.3);\n"
"\n"
"vec3 color=vec3(0.),transm=vec3(1.);\n"
"#if 0\n"
"hit_t t=trace_grid(O,D,FAR);\n"
"CHECK_ASSERT\n"
"gl_FragColor=vec4(quantize(vec3(0.),vec3(256.),8,vec3(t._gc,t._mc,t._mc+t._gc)),1.);return;\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"\n"
"if(t.mid > 0){\n"
"color=t.n*4.;\n"
"}else{\n"
"\n"
"\n"
"color=vec3(1000.,0.,0.);\n"
"}\n"
"#else\n"
"float Lmax=FAR;\n"
"for(int i=0;i<BOUNCES;++i){\n"
"if(dot(transm,transm)<.001)break;\n"
"hit_t h=trace_grid(O,D,Lmax);\n"
"CHECK_ASSERT\n"
"\n"
"if(h.mid > 0){\n"
"mat_t m=material(h);\n"
"vec3 c=m.e;\n"
"O=h.p+h.n*EPS;\n"
"\n"
"if(trace_grid(O,sundir,100.).l >=100.)c+=brdf(h,m,sundir)*air(O,sundir);\n"
"D=ref(h);\n"
"color+=transm*c;\n"
"transm*=brdf(h,m,D);\n"
"Lmax*=.6;\n"
"}else{\n"
"color+=transm*air(O,D);\n"
"}\n"
"}\n"
"\n"
"#endif\n"
"gl_FragColor=vec4(color,1.);\n"
"}\n"
, 
"uniform float _t,_p,_pt;\n"
"uniform vec2 _r;\n"
"uniform vec3 _s,_cp;\n"
"uniform mat3 _cm;\n"
"uniform sampler2D _N,_T,_P,_F;\n"
"varying vec2 V;\n"
"\n"
"vec4 noise(vec2 p){return texture2D(_N,(p+.5)/1024.,-20.);}\n"
"vec4 terrain(vec2 p_m){return texture2D(_T,p_m/4096.,-20.);}\n"
"vec4 plan(vec2 p_m){return texture2D(_P,(floor(p_m/32.)+.5)/128.,-20.);}\n"
"\n"
"float t=_t;\n"
"float h(vec2 p){\n"
"vec4 c=texture2D(_T,p/4096.,-20.);\n"
"return c.z;\n"
"\n"
"}\n"
"float h2(vec2 p){\n"
"vec4 c=terrain(p);\n"
"\n"
"return c.z;\n"
"}\n"
"vec3 n(vec2 p){\n"
"vec2 e=vec2(.5,0.);\n"
"vec3 dx=vec3(2.*e.x,h(p+e.xy)-h(p-e.xy),0.);\n"
"vec3 dz=vec3(0.,h(p+e.yx)-h(p-e.yx),2.*e.x);\n"
"return normalize(cross(dz,dx));\n"
"}\n"
"vec3 albedo(vec3 p){\n"
"#if 0\n"
"float m10=mod(floor(p.x/10.)+floor(p.z/10.),2.);\n"
"float m100=mod(floor(p.x/100.)+floor(p.z/100.),2.);\n"
"float m1000=mod(floor(p.x/1000.)+floor(p.z/1000.),2.);\n"
"return vec3(m10,m100,m1000)*(.3+.7*mod(floor(p.x)+floor(p.z),2.));\n"
"#elif 0\n"
"float XZ=mod(floor(p.x/1.)+floor(p.z/1.),2.);\n"
"float Y=mod(floor(p.y),2.);\n"
"return vec3(.2+.5*Y+.3*XZ);\n"
"#else\n"
"vec4 T=terrain(p.xz);\n"
"return vec3(.2)+vec3(T.x,T.y,T.w);\n"
"#endif\n"
"}\n"
"vec4 trace(vec2 uv){\n"
"vec3 O=_cp,D=normalize(vec3(uv,2.))*_cm;\n"
"O.y=max(O.y,h2(O.xz)+10.);\n"
"float l=0.,lp=l;\n"
"for(int i=0;i<128;++i){\n"
"vec3 p=O+D*l;\n"
"float d=p.y-h(p.xz);\n"
"if(d<.001*l){\n"
"for(int j=0;j<9;++j){\n"
"float lm=(l+lp)*.5;\n"
"vec3 p=O+D*lm;\n"
"float d=p.y-h(p.xz);\n"
"if(d<.001*lm)l=lm;else lp=lm;\n"
"}\n"
"break;\n"
"}\n"
"lp=l;\n"
"\n"
"l+=d;\n"
"if(l>3000.)break;\n"
"}\n"
"if(l>3000.){return vec4(0.);}\n"
"vec3 p=O+D*l;\n"
"if(p.y-h(p.xz)>.01*l){return vec4(1.,0.,0.,1.);}\n"
"if(h(p.xz)>p.y){return vec4(1.,0.,1.,1.);}\n"
"vec3 n=n((O+D*l).xz);\n"
"vec3 color=albedo(p)*(max(0.,dot(n,normalize(_s)))+vec3(.05));\n"
"return vec4(pow(color,vec3(1./2.2)),1.);\n"
"}\n"
"\n"
"void main(){\n"
"vec2 fc=floor(gl_FragCoord.xy);\n"
"vec2 uv=fc/_r;\n"
"if(fc.y<16.){gl_FragColor=vec4(step(V.x*.5+.5,_p));return;}\n"
"#if 0\n"
"vec2 pt=floor(uv*4096.);\n"
"vec2 pp=floor(uv*128.);\n"
"vec4 T=texture2D(_T,(pt+vec2(.5))/4096.,-20.);\n"
"#if 0\n"
"gl_FragColor=T*vec4(1.,1.,.005,1.);\n"
"#else\n"
"vec4 P=texture2D(_P,(pp+vec2(.5))/128.,-20.);\n"
"float dH=P.w-T.z;\n"
"float dh=T.z-P.z;\n"
"float dHB=P.y-P.w;\n"
"\n"
"gl_FragColor=vec4(-dH*1e7,-dh*1e7,dHB,0.);\n"
"#endif\n"
"#else\n"
"if(_pt <.0)\n"
"gl_FragColor=trace((uv-vec2(.5))*vec2(_r.x/_r.y,1.));\n"
"else \n"
"gl_FragColor=pow(texture2D(_F,uv,-20.)*.6/(1.+63.*_pt),vec4(1./2.2));\n"
"#endif\n"
"}\n"
, 
};
