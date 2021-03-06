uniform float _t, _p;
uniform vec2 _r;
uniform sampler2D _N,_T,_P,_F;
varying vec2 V;

vec4 n4(vec2 p){return texture2D(_N,(p+vec2(.5))/1024.,-20.);}

void main(){
  vec2 fc=floor(gl_FragCoord.xy);
  vec2 p = fc*32.+vec2(.5);
  float h=1e6,H=0.;
  for(int x=0;x<32;++x){
    for(int y=0;y<32;++y){
      vec2 P = p + vec2(float(x), float(y));
      vec4 s = texture2D(_T,P/4096.,-20.);
      H=max(H,s.z);h=min(h,s.z);
    }
  }
  float d = H - h;
  vec4 r = n4(fc);
  float pop = pow(smoothstep(25.,0.,d)*smoothstep(230.,110.,H),14.);
  pop *= .2+.8*r.z;
  pop += step(0., pop) * .08;
  gl_FragColor=vec4(r.y,max(0.,pop*90.)+H,h,H);
  //gl_FragColor=vec4(0.,H,h,H);
}
