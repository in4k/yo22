void main() {
  vec2 uv=fc/4096.;
  float h = fbm(uv,5.);
  h = pow(h, 3.) * 1200.; //+ 20. * step(.7,n4(uv*16.).x);
  gl_FragColor = vec4(0., 0., h, 0.);
}
