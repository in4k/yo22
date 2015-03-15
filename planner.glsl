uniform float _t, _p;
uniform vec2 _r;
uniform sampler2D _N,_T,_P,_F;
varying vec2 V;

void main(){
  vec2 p=floor(gl_FragCoord.xy)*32.+vec2(.5);
  float h=1e6,H=0.;
  for(int x=0;x<32;++x){
    for(int y=0;y<32;++y){
      vec2 P = p + vec2(float(x), float(y));
      vec4 s = texture2D(_T,P/4096.,-20.);
      H=max(H,s.z);h=min(h,s.z);
    }
  }
  //gl_FragColor=vec4(0.,step(20.,h-H)*20. + h,h,H);
  gl_FragColor=vec4(0.,H,h,H);
}
