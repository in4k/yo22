uniform float _t,_p;
uniform vec2 _r;
uniform vec3 _s,_cp;
uniform mat3 _cm;
uniform sampler2D _N,_T,_P,_F;
varying vec2 V;
vec2 fc=floor(gl_FragCoord.xy);
vec4 n4(vec2 v){return texture2D(_N,(v+vec2(.5))/256.,-20.);}
float n(vec2 v){return n4(v).w;}
float vn(vec2 v,vec2 m) {
  vec2 e=vec2(1.,0.),V=floor(v);v=fract(v);v*=v*(3.-2.*v);
  return mix(mix(n(mod(V+e.yy,m)),n(mod(V+e.xy,m)),v.x),mix(n(mod(V+e.yx,m)),n(mod(V+e.xx,m)),v.x),v.y);
}
float fbm(vec2 v,float s) {
  float r=0.,k=.5;
  for(int i=0;i<12;++i,k*=.5,s*=2.)r+=k*vn(v*s,vec2(s));
  return r;
}
int rs=int(fc.y*1280.+fc.x)+int(n4(vec2(_t,0.)).x*256.+n4(vec2(0.,_t)).z)*256;
vec4 rand(){rs=int(mod(float(rs+1),1024.*1024.));return n4(vec2(float(rs),floor(float(rs)/1024.)));}
vec4 TT(vec2 p){return texture2D(_T,(p+vec2(.5))/4096.,-20.);}
vec4 PP(vec2 p){return texture2D(_P,(floor(p/32.)+vec2(.5))/128.,-20.);}