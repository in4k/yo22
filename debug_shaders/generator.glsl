uniform float _t, _p;
uniform vec2 _r;
uniform sampler2D _N,_T,_P,_F;
varying vec2 V;

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
void main() {
  vec2 uv=floor(gl_FragCoord.xy)/4096.;
  float height = fbm(uv,5.);
  height = pow(height, 3.) * 1200.; //+ 20. * step(.7,n4(uv*16.).x);
  gl_FragColor = vec4(0., 0., height, 0.);
}
