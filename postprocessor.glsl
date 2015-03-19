void main(){
  vec2 uv=fc/_r;
  if (_p < 1.)
    gl_FragColor = vec4(mix(0.,step(uv.x,_p),step(uv.y,.025))); // loading
  else
    gl_FragColor = pow(texture2D(_F,uv,-20.)*.7/32.,vec4(1./2.2));
  //  pow(texture2D(_F,uv,-20.)/64.,vec4(1./2.2)), // final
  //  step(1.,_p)
  //  );
  gl_FragColor = pow(texture2D(_F,uv,-20.)*.7/(1.+31.*_pt),vec4(1./2.2));
}
