void main(){
  vec2 uv=fc/_r;
  //gl_FragColor = mix(
  //  vec4(mix(0.,step(uv.x,_p),step(uv.y,.025))), // loading
  //  pow(texture2D(_F,uv,-20.)/64.,vec4(1./2.2)), // final
  //  step(1.,_p)
  //  );
  gl_FragColor = pow(texture2D(_F,uv,-20.)*.5/(1.+63.*_pt),vec4(1./2.2));
}
