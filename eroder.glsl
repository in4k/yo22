float
Kr = .27,
Ke = .005,
Ks = .8;
vec2 uv=floor(gl_FragCoord.xy);
float tf(float H,float w,vec2 o){
  vec4 c=TT(uv+o);
  float h=c.w+c.z,d=H-h;
  return mix(-min(w,d),min(c.w,-d),step(d,.0));
}
void main() {
  vec4 c=TT(uv);
  //gl_FragColor=c;return;
  // 1. tf
  vec4 c0=c;
  float h=c.w+c.z,dw=0.;
  vec3 e=vec3(1.,0.,-1.);
  float kk = .25;
  dw+=kk*tf(h,c.w,e.xy);
  dw+=kk*tf(h,c.w,e.yx);
  dw+=kk*tf(h,c.w,-e.xy);
  dw+=kk*tf(h,c.w,-e.yx);
  //kk*=1./sqrt(2.);
  //dw+=kk*tf(h,c.w,e.xx);
  //dw+=kk*tf(h,c.w,e.xz);
  //dw+=kk*tf(h,c.w,e.zx);
  //dw+=kk*tf(h,c.w,e.zz);
  c.w+=dw;
  // 2. evaporate
  //dw=min(c.w,Ke);
  dw=c.w*Ke;
  c.z+=dw*Ks;
  c.w-=dw;
  // 3. rain
  float sa=sin(_t*373.4111),ca=cos(_t*373.4111);
  mat2 m=mat2(sa,ca,ca,-sa);
  dw=Kr*n4((uv+vec2(_t*123.24,317.1141*_t))*m).x*(0==mod(floor(_t/3.),12.)?.1:.0);
  //dw=Kr*n4((uv+vec2(_t*123.24,317.1141*_t))*m)*(0==mod(floor(_t/3.),12.)?1.:0.);
  c.z-=dw*Ks;
  c.w+=dw;
  // 4. PROFIT
  gl_FragColor=vec4(
    1000.*(c.w-c0.w),
    -1000.*(c.z-c0.z),
    c.z,
    c.w);
}
