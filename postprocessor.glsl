uniform float P,t;
uniform sampler2D R;
varying vec2 V;
void main() {
  gl_FragColor = vec4(step(V.y,-.9) * step(V.x*.5+.5,P) + texture2D(R,V*(sin(t)+1.)).x);
}
