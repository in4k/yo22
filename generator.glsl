void main() {
  vec2 uv=floor(gl_FragCoord.xy)/4096.;
  float height = fbm(uv,5.);
  height = pow(height, 3.) * 1200.; //+ 20. * step(.7,n4(uv*16.).x);
  gl_FragColor = vec4(0., 0., height, 0.);
}
