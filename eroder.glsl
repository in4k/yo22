uniform float _p;
uniform sampler2D _N,_T;
varying vec2 V;
float
Kr = 0.1,
Ke = 0.6,
Ks = 1.0;
float transfer(float H,float w,vec2 ofs){
  vec4 c=texture2D(_T,(gl_FragCoord.xy+ofs)/4096.,-20.);
  float h=c.w+c.z,d=H-h;
  if(d>0.)return -min(w,d);
  else return min(c.w,-d);
}
void main(){
  vec4 c=texture2D(_T,(gl_FragCoord.xy)/4096.,-20.);
  // 1. transfer
  vec4 c0=c;
  float h=c.w+c.z,dw=0.;
  vec2 e=vec2(1.,0.);
  dw+=transfer(h,c.w,e.xy);
  dw+=transfer(h,c.w,e.yx);
  dw+=transfer(h,c.w,-e.xy);
  dw+=transfer(h,c.w,-e.yx);
  dw+=transfer(h,c.w,vec2(1.,1.));
  dw+=transfer(h,c.w,vec2(1.,1.));
  dw+=transfer(h,c.w,-vec2(1.,1.));
  dw+=transfer(h,c.w,-vec2(1.,1.));
  c.w+=dw/8.;
  // 2. evaporate
  dw=c.w*Ke;
  c.z+=dw*Ks;
  c.w-=dw;
  // 3. rain
  dw=Kr;
  c.z-=dw*Ks;
  c.w+=dw;
  // 4. PROFIT
  gl_FragColor=vec4(
    //(mod(_p*20.,2.)*2.-1.)
    1000.*(c.w-c0.w),
    //c.w-Kr,
    -1000.*(c.z-c0.z),
    c.z,
    c.w);
}
