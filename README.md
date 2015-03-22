# yo22
Terrain and city generator in less than 4096 bytes

Made for 4k procedural graphics compo at [NVScene 2015](http://nv.scene.org/), where it took 3rd place.

Pouet link: http://www.pouet.net/prod.php?which=65263

Demozoo link: http://demozoo.org/productions/135216/

Requires:
 - OpenGL 2.1
 - Fast GPU; takes ~10 seconds to generate multisample images, or about 100..1000 ms for individual samples on Titan-class cards

Contains:
- compute pipeline using fragment shaders, with textures bound to framebuffers as in/out
- multioctave value noise terrain generator and hydraulic erosion (erosion is turned off for the party version build due to size limitations)
- mixed heightmap and signed distance field raymarcher
- atmosphere scattering (rayleigh + mie)
- monte-carlo global illumination (path tracing)

Can be used to produce images like this:
![](http://yolp.omgwtf.ru/img/yo22_stage_7.jpg)

or this:
![](http://yolp.omgwtf.ru/img/yo22_stage_4.jpg)

Supported platforms:
 - Windows (target compo platform): viewer only
 - Mac OS X (used for development):  shader livecoding tool
 - Linux (used for development): shader livecoding tool and wsad+mouselook

I have no plans to continue development or support this, but may occasionally play with it and tweak things out of boredom.

Released under WTFPL.
