#ifdef __linux__
#define PLATFORM_POSIX
#define PLATFORM_LINUX
#define PLATFORM_X11
#elif defined(__MACH__) && defined(__APPLE__)
#define PLATFORM_POSIX
#define PLATFORM_OSX
#elif defined(_WIN32)
#define PLATFORM_WINDOWS
#else
#error Not ported
#endif

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <errno.h>

#ifdef PLATFORM_POSIX
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

#ifdef PLATFORM_OSX
#include <sys/event.h>
#include <GLUT/GLUT.h>
#include <OpenGL/gl3.h>
#endif

#ifdef PLATFORM_LINUX
#include <sys/inotify.h>
#endif

#ifdef PLATFORM_X11
#define GL_GLEXT_PROTOTYPES 1
#include <X11/Xlib.h>
#include <GL/glx.h>
#include <GL/gl.h>
#include <GL/glext.h>
#endif

#ifdef PLATFORM_WINDOWS
#define WIN32_LEAN_AND_MEAN 1
#define NOMINMAX
#include <Windows.h>
#include <GL/gl.h>
#include "glext.h"
#endif

#define FUNCLIST \
  FUNCLIST_DO(PFNGLACTIVETEXTUREPROC, ActiveTexture) \
  FUNCLIST_DO(PFNGLCREATESHADERPROC, CreateShader) \
  FUNCLIST_DO(PFNGLSHADERSOURCEPROC, ShaderSource) \
  FUNCLIST_DO(PFNGLCOMPILESHADERPROC, CompileShader) \
  FUNCLIST_DO(PFNGLGENFRAMEBUFFERSPROC, GenFramebuffers) \
  FUNCLIST_DO(PFNGLBINDFRAMEBUFFERPROC, BindFramebuffer) \
  FUNCLIST_DO(PFNGLFRAMEBUFFERTEXTURE2DPROC, FramebufferTexture2D) \
  FUNCLIST_DO(PFNGLUSEPROGRAMPROC, UseProgram) \
  FUNCLIST_DO(PFNGLUNIFORM1FPROC, Uniform1f) \
  FUNCLIST_DO(PFNGLUNIFORM2FPROC, Uniform2f) \
  FUNCLIST_DO(PFNGLUNIFORM1IPROC, Uniform1i) \
  FUNCLIST_DO(PFNGLCREATEPROGRAMPROC, CreateProgram) \
  FUNCLIST_DO(PFNGLATTACHSHADERPROC, AttachShader) \
  FUNCLIST_DO(PFNGLBINDATTRIBLOCATIONPROC, BindAttribLocation) \
  FUNCLIST_DO(PFNGLLINKPROGRAMPROC, LinkProgram) \
  FUNCLIST_DO(PFNGLGETUNIFORMLOCATIONPROC, GetUniformLocation)
#ifndef DEBUG
#define FUNCLIST_DBG
#else
#define FUNCLIST_DBG \
  FUNCLIST_DO(PFNGLDELETESHADERPROC, DeleteShader) \
  FUNCLIST_DO(PFNGLGETPROGRAMINFOLOGPROC, GetProgramInfoLog) \
  FUNCLIST_DO(PFNGLDELETEPROGRAMPROC, DeleteProgram) \
  FUNCLIST_DO(PFNGLGETSHADERIVPROC, GetShaderiv) \
  FUNCLIST_DO(PFNGLGETSHADERINFOLOGPROC, GetShaderInfoLog) \
  FUNCLIST_DO(PFNGLGETPROGRAMIVPROC, GetProgramiv) \
  FUNCLIST_DO(PFNGLCHECKFRAMEBUFFERSTATUSPROC, CheckFramebufferStatus) \
  FUNCLIST_DO(PFNGLENABLEVERTEXATTRIBARRAYPROC, EnableVertexAttribArray) \
  FUNCLIST_DO(PFNGLVERTEXATTRIBPOINTERPROC, VertexAttribPointer)
#endif

#ifdef PLATFORM_POSIX
static void report_n_abort(const char *file, int line, const char *message) {
  fprintf(stderr, "error @ %s:%d : %s\n", file, line, message);
  exit(-1);
}
#endif

#ifdef PLATFORM_WINDOWS
#define FUNCLIST_DO(T,N) "gl" #N "\0"
static const char *gl_names =
FUNCLIST FUNCLIST_DBG
;
#undef FUNCLIST_DO

static struct {
#define FUNCLIST_DO(T,N) T N;
FUNCLIST FUNCLIST_DBG
#undef FUNCLIST_DO
} gl;

#if DEBUG
static void report_n_abort(const char *file, int line, const char *message) {
  char buf[256];
  _snprintf_s(buf, sizeof(buf), sizeof(buf)-1, "error @ %s:%d : %s\n", file, line, message);
  MessageBoxA(0, message, buf, MB_ICONSTOP);
  ExitProcess(-1);
}
#endif
#else // not windows
static struct {
#define FUNCLIST_DO(T,N) T N;
FUNCLIST FUNCLIST_DBG
#undef FUNCLIST_DO
} gl = {
#define FUNCLIST_DO(T,N) gl##N,
FUNCLIST FUNCLIST_DBG
#undef FUNCLIST_DO
};
#endif

#ifdef DEBUG
#define CHECK(cond, errmsg) if (!(cond)) report_n_abort(__FILE__, __LINE__, errmsg);
#else
#define CHECK(cond, errmsg)
#endif

#define WIDTH 1280
#define HEIGHT 720

#define NOISE_SIZE 1024
#define TERRAIN_SIZE 4096

#define TERRAIN_ITERATIONS 8/*256*/
#define TRACE_ITERATIONS 1024
#define TOTAL_ITERATIONS (TERRAIN_ITERATIONS + TRACE_ITERATIONS)



static unsigned int noise_buffer[NOISE_SIZE * NOISE_SIZE];

static unsigned int prng_state = 5323u;

static unsigned int prng__() {
  prng_state = (1103515245u * prng_state + 12345u);
  return prng_state >> 16;
}

static unsigned int prng() {
  unsigned int s0 = prng__();
  unsigned int s1 = prng__();
  unsigned int s2 = prng__();
  return (s0 | (s1 << 16)) ^ (s2 << 8);
}

enum {
  TexNoise8U,
  TexTerrain0, TexTerrain1,
  TexFrame,
  Tex_COUNT
};

static int tex_widths[Tex_COUNT] = {
  NOISE_SIZE, TERRAIN_SIZE, TERRAIN_SIZE, WIDTH
};
static int tex_heights[Tex_COUNT] = {
  NOISE_SIZE, TERRAIN_SIZE, TERRAIN_SIZE, HEIGHT
};
static void *tex_data[Tex_COUNT] = {
  noise_buffer, NULL, NULL, NULL
};

enum {
  ProgTerrainGenerate,
  ProgTerrainErode,
  ProgTrace,
  ProgPostprocess,
  Prog_COUNT
};

#ifndef TOOL
#include "shaders.h"
#else
static const char *fragment_shaders[Prog_COUNT] = {NULL, NULL, NULL, NULL};
#endif

static const char *vertex_shader_source[] = {"attribute vec4 v;varying vec2 V;void main(){gl_Position=v;V=v.xy;}"};
/*static char *common_shader_header = 0;*/

static GLuint textures[Tex_COUNT];
static GLuint framebuffer;
static GLuint vertex_shader;
static GLuint programs[Prog_COUNT];
static struct {
  GLint time, progress, progress_erosion, progress_trace, target_res, noise, terrain, frame;
} program_locs[Prog_COUNT];
static const char *uniform_names[] = {
  "_t", "_p", "_pe", "_pt", "_r", "_N", "_T", "_F", 0
};

static void bind_texture(int index, int sampler) {
  gl.ActiveTexture(GL_TEXTURE0 + sampler);
  glBindTexture(GL_TEXTURE_2D, textures[index]);
}

static void bind_framebuffer(int target_index) {
  gl.BindFramebuffer(GL_FRAMEBUFFER, framebuffer);
  gl.FramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textures[target_index], 0);
  CHECK(GL_FRAMEBUFFER_COMPLETE == gl.CheckFramebufferStatus(GL_FRAMEBUFFER), "framebuffer is not complete");
  glViewport(0, 0, tex_widths[target_index], tex_heights[target_index]);
}

static GLuint program;

static float u_time = 0, u_progress, u_progress_erosion, u_progress_trace;
static const int noise_sampler = 0, terrain_sampler = 1, frame_sampler = 2;
static int width, height;

static void use_program(int index) {
  program = programs[index];
  gl.UseProgram(program);
  gl.Uniform1f(program_locs[index].time, u_time*.001);
  gl.Uniform1f(program_locs[index].progress, u_progress);
  gl.Uniform1f(program_locs[index].progress_erosion, u_progress_erosion);
  gl.Uniform1f(program_locs[index].progress_trace, u_progress_trace);
  gl.Uniform2f(program_locs[index].target_res, (float)width, (float)height);
  gl.Uniform1i(program_locs[index].noise, noise_sampler);
  gl.Uniform1i(program_locs[index].terrain, terrain_sampler);
  gl.Uniform1i(program_locs[index].frame, frame_sampler);
}

static void compute() {
  glRects(-1, -1, 1, 1);
}

static GLuint create_and_compile_shader(int type, int n, const char **source) {
  GLuint shader = gl.CreateShader(type);
  CHECK(shader != 0, "glCreateShader");
  gl.ShaderSource(shader, n, source, NULL);
  gl.CompileShader(shader);

#ifdef DEBUG
  {
    GLint status;
    gl.GetShaderiv(shader, GL_COMPILE_STATUS, &status);
    if (status != GL_TRUE) {
      char errorbuf[4096];
      gl.GetShaderInfoLog(shader, sizeof(errorbuf), NULL, errorbuf);
      fprintf(stderr, "compile error:\n%s\n", errorbuf);
      gl.DeleteShader(shader);
      shader = 0;
    }
  }
#endif
  return shader;
}

static GLuint create_and_compile_program(int index) {
  int i;
#ifdef DEBUG
  GLint status;
  if (fragment_shaders[index] == 0) return 0;
#endif

  /*const char *sources[] = {common_shader_header, fragment_source};
  GLuint fragment = create_and_compile_shader(GL_FRAGMENT_SHADER, 2, sources);*/
  GLuint fragment = create_and_compile_shader(GL_FRAGMENT_SHADER, 1, fragment_shaders + index);
  GLuint program = gl.CreateProgram();
  gl.AttachShader(program, fragment);
  gl.AttachShader(program, vertex_shader);
  gl.BindAttribLocation(program, 0, "v");
  gl.LinkProgram(program);

#ifdef DEBUG
  gl.DeleteShader(fragment);
  gl.GetProgramiv(program, GL_LINK_STATUS, &status);
  if (status != GL_TRUE) {
    char errorbuf[4096];
    gl.GetProgramInfoLog(program, sizeof(errorbuf), NULL, errorbuf);
    fprintf(stderr, "link error:\n%s\n", errorbuf);
    gl.DeleteProgram(program);
    gl.DeleteShader(fragment);
    return 0;
  }
  fprintf(stderr, "Linked\n");

#endif
    for (i = 0; uniform_names[i] != 0; ++i) {
    ((GLint*)(program_locs + index))[i] = gl.GetUniformLocation(program, uniform_names[i]);
#ifdef DEBUG
    fprintf(stderr, "  %s loc = %d  ", uniform_names[i], ((GLint*)(program_locs + index))[i]);
#endif
  }
#ifdef DEBUG
  fprintf(stderr, "\n");
  gl.DeleteProgram(programs[index]);
#endif
  programs[index] = program;
  return program;
}

static void yo22_init() {
  int i;
  
  for (i = 0; i < NOISE_SIZE * NOISE_SIZE; ++i)
    noise_buffer[i] = prng();

  glGenTextures(Tex_COUNT, textures);
  gl.GenFramebuffers(1, &framebuffer);
  for (i = 0; i < Tex_COUNT; ++i) {
    glBindTexture(GL_TEXTURE_2D, textures[i]);
    glTexImage2D(GL_TEXTURE_2D, 0, (i==TexNoise8U)?GL_RGBA:GL_RGBA32F,
        tex_widths[i], tex_heights[i], 0, GL_RGBA, GL_UNSIGNED_BYTE, tex_data[i]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  }

  vertex_shader = create_and_compile_shader(GL_VERTEX_SHADER, 1, vertex_shader_source);
  CHECK(vertex_shader != 0, "vertex shader");
  for (i = 0; i < Prog_COUNT; ++i)
    create_and_compile_program(i);
}

/*
static float camera[16];
static void cross(const float *a, const float *b, float *out) {
  out[0] = a[1] * b[2] - a[2] * b[1];
  out[1] = a[2] * b[0] - a[0] * b[2];
  out[2] = a[0] * b[1] - a[1] * b[0];
}
static float dot(const float *a, const float *b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
static float sqrtf(float f){asm("fld %0;fsqrt;fstp %0;":"+m"(f));return f;}
static void normalize(float *v) {
  float k = 1. / sqrtf(dot(v,v));
  v[0] *= k; v[1] *= k; v[2] *= k;
}
static void camera_set(float px, float py, float pz, float ax, float ay, float az) {
}
static void camera_move(float fwd, float right, float up) {
}
static void camera_rotate_pitch(float fwd) {
}
static void camera_rotate_yaw(float fwd) {
}
*/

static void yo22_size(int w, int h) {
  width = w, height = h;
}

static int program_counter = 0;

static void yo22_paint() {
  u_time += 1;
  u_progress = (float)program_counter / TOTAL_ITERATIONS;
  u_progress_erosion = program_counter < TERRAIN_ITERATIONS ? (float)program_counter / TERRAIN_ITERATIONS : 1.f;
  u_progress_trace = program_counter < TRACE_ITERATIONS ? (float)(program_counter - TERRAIN_ITERATIONS) / (TRACE_ITERATIONS-TERRAIN_ITERATIONS) : 1.f;

  bind_texture(TexNoise8U, noise_sampler);
  if (program_counter == 0) {
    bind_framebuffer(TexTerrain0);
    use_program(ProgTerrainGenerate);
    compute();
  }

  if (program_counter < TERRAIN_ITERATIONS) {
    int j;for(j=0;j<8;++j){
    bind_framebuffer(TexTerrain1);
    use_program(ProgTerrainErode);
    bind_texture(TexTerrain0, terrain_sampler);
    compute();
    bind_framebuffer(TexTerrain0);
    bind_texture(TexTerrain1, terrain_sampler);
    compute();}
  }

  if (program_counter >= TERRAIN_ITERATIONS && program_counter < TOTAL_ITERATIONS) {
    bind_framebuffer(TexFrame);
    if (program_counter == TERRAIN_ITERATIONS) {
      glClearColor(0,0,0,0);
      glClear(GL_COLOR_BUFFER_BIT);
    }
    use_program(ProgTrace);
    bind_texture(TexTerrain0, terrain_sampler);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    compute();
    glDisable(GL_BLEND);
  }

  gl.BindFramebuffer(GL_FRAMEBUFFER, 0);
  glViewport(0, 0, width, height);
  use_program(ProgPostprocess);
  bind_texture(TexFrame, frame_sampler);
  bind_texture(TexTerrain0, terrain_sampler);
  compute();

  ++program_counter;
  if (program_counter >= TOTAL_ITERATIONS)
    program_counter = TOTAL_ITERATIONS;
  //else
  //  fprintf(stderr, "%d%c", program_counter, program_counter%16==15 ? '\n' : ' ');
}

#ifdef PLATFORM_POSIX
static char *read_file(const char *filename) {
  ssize_t rd;
  struct stat st;
  char *ret = NULL;
  int fd = open(filename, 0);
  if (fd == -1) {
    fprintf(stderr, "error: cannot open file \"%s\"\n", filename);
    goto error;
  }

  if (-1 == fstat(fd, &st)) {
    fprintf(stderr, "error: cannot read file \"%s\"\n", filename);
    goto error;
  }

  if (st.st_size < sizeof("void main(){gl_FragColor=vec4(0.);}")) {
    fprintf(stderr, "error: file \"%s\" has unreasonable size (%d bytes)\n",
      filename, (int)st.st_size);
    goto error;
  }

  ret = malloc(st.st_size + 1);

  rd = read(fd, ret, st.st_size);
  if (rd != st.st_size) {
    fprintf(stderr, "error: cannot read file \"%s\"\n", filename);
    goto error;
  }

  ret[st.st_size] = 0;
  close(fd);
  return ret;

error:
  if (ret != NULL) free(ret);
  if (fd != -1) close(fd);
  return NULL;
}
#endif

#ifdef PLATFORM_LINUX
struct file_program_t {
  const char *filename;
  int program_counter;
  int watch;
  int updated;
};
static struct file_program_t files[Prog_COUNT];
static int inotifyfd;
static void monitor_changes() {
  char buffer[sizeof(struct inotify_event) + NAME_MAX + 1];
  int i;

  for (i = 0; i < Prog_COUNT; ++i)
    if (files[i].watch == -1) {
      files[i].watch = inotify_add_watch(inotifyfd, files[i].filename,
        IN_MODIFY | IN_DELETE_SELF | IN_MOVE_SELF);
      files[i].updated = (files[i].watch == -1) ? 0 : 1;
    }

  for (;;) {
    const struct inotify_event *e = (const struct inotify_event*)buffer;
    ssize_t rd = read(inotifyfd, buffer, sizeof(buffer));
    if (rd == -1) {
      CHECK(errno == EAGAIN, "read inotify fd, errno != EAGAIN");
      break;
    }

    while (rd >= sizeof(struct inotify_event)) {
      int evt_size = sizeof(struct inotify_event) + e->len;
      struct file_program_t *p = 0;

      for (i = 0; i < Prog_COUNT; ++i)
        if (files[i].watch == e->wd)
          p = files + i;

      fprintf(stderr, "%d %d %p\n", i, e->wd, p);

      if (p) {
        p->updated = 1;

        if (e->mask & (IN_IGNORED | IN_DELETE_SELF | IN_MOVE_SELF)) {
          if (e->mask & IN_MOVE_SELF) inotify_rm_watch(inotifyfd, p->watch);
          p->watch = -1;
          p->updated = 0;
        }
      } // if (p)

      rd -= evt_size;
      e = (struct inotify_event*)(((char*)e) + evt_size);
    } // while (events)
  }

  for (i = 0; i < Prog_COUNT; ++i)
    if (files[i].updated != 0) {
      char *src = read_file(files[i].filename);
      files[i].updated = 0;
      if (!src) continue;
      if (fragment_shaders[i] != NULL)
        free((char*)fragment_shaders[i]);
      fragment_shaders[i] = src;
      fprintf(stderr, "loading program %d ... ", i);
      if (0 == create_and_compile_program(i)) continue;
      if (program_counter > files[i].program_counter)
        program_counter = files[i].program_counter;
    }
}
#endif // PLATFORM_LINUX

#ifdef PLATFORM_X11
static const int glxattribs[] = {
  GLX_X_RENDERABLE, True,
  GLX_DRAWABLE_TYPE, GLX_WINDOW_BIT,
  GLX_RENDER_TYPE, GLX_RGBA_BIT,
  GLX_CONFIG_CAVEAT, GLX_NONE,
  GLX_RED_SIZE, 8,
  GLX_GREEN_SIZE, 8,
  GLX_BLUE_SIZE, 8,
  GLX_ALPHA_SIZE, 8,
  GLX_DOUBLEBUFFER, True,
  0
};
int main(int argc, char *argv[]) {
  Display *display;
  Window window;
  GLXDrawable drawable;
  GLXContext context;
  int nglxconfigs = 0;
  GLXFBConfig *glxconfigs = NULL;
  XVisualInfo *vinfo = NULL;
  XSetWindowAttributes winattrs;
  Atom delete_message;

/*  if (argc != 5) {
    fprintf(stderr, "Usage: %s generator.glsl eroder.glsl tracer.glsl postprocessor.glsl\n", argv[0]);
    return 1;
  }*/
  files[ProgTerrainGenerate].filename = "generator.glsl";/*argv[1];*/
  files[ProgTerrainGenerate].program_counter = 0;
  files[ProgTerrainGenerate].watch = -1;
  files[ProgTerrainErode].filename = "eroder.glsl";/*argv[2];*/
  files[ProgTerrainErode].program_counter = 0;
  files[ProgTerrainErode].watch = -1;
  files[ProgTrace].filename = "tracer.glsl";/*argv[3];*/
  files[ProgTrace].program_counter = TERRAIN_ITERATIONS;
  files[ProgTrace].watch = -1;
  files[ProgPostprocess].filename = "postprocessor.glsl"/*argv[4]*/;
  files[ProgPostprocess].program_counter = TOTAL_ITERATIONS;
  files[ProgPostprocess].watch = -1;
  inotifyfd = inotify_init1(IN_NONBLOCK);
  CHECK(inotifyfd != -1, "inotify_init1");

  display = XOpenDisplay(NULL);
  CHECK(display, "XOpenDisplay");

  glxconfigs = glXChooseFBConfig(display, 0, glxattribs, &nglxconfigs);
  CHECK(glxconfigs && nglxconfigs, "glXChooseFBConfig");

  vinfo = glXGetVisualFromFBConfig(display, glxconfigs[0]);
  CHECK(vinfo, "glXGetVisualFromFBConfig");

  memset(&winattrs, 0, sizeof(winattrs));
  winattrs.event_mask =
    ExposureMask | VisibilityChangeMask | StructureNotifyMask |
    KeyPressMask | PointerMotionMask;
  winattrs.border_pixel = 0;
  winattrs.bit_gravity = StaticGravity;
  winattrs.colormap = XCreateColormap(display,
    RootWindow(display, vinfo->screen),
    vinfo->visual, AllocNone);
  winattrs.override_redirect = False;

  window = XCreateWindow(display, RootWindow(display, vinfo->screen),
    0, 0, 1280, 720,
    0, vinfo->depth, InputOutput, vinfo->visual,
    CWBorderPixel | CWBitGravity | CWEventMask | CWColormap,
    &winattrs);
  CHECK(window, "XCreateWindow");

  XStoreName(display, window, "shapa");

  delete_message = XInternAtom(display, "WM_DELETE_WINDOW", True);
  XSetWMProtocols(display, window, &delete_message, 1);

  XMapWindow(display, window);

  context = glXCreateNewContext(display, glxconfigs[0], GLX_RGBA_TYPE, 0, True);
  CHECK(context, "glXCreateNewContext");

  drawable = glXCreateWindow(display, glxconfigs[0], window, 0);
  CHECK(drawable, "glXCreateWindow");

  glXMakeContextCurrent(display, drawable, drawable, context);

  XSelectInput(display, window, StructureNotifyMask | KeyReleaseMask);

  yo22_init();

  for (;;) {
    while (XPending(display)) {
      XEvent e;
      XNextEvent(display, &e);
      switch (e.type) {
        case ConfigureNotify:
          yo22_size(e.xconfigure.width, e.xconfigure.height);
          break;

        case KeyRelease:
          if (XLookupKeysym(&e.xkey, 0) != XK_Escape)
            break;
        case ClientMessage:
        case DestroyNotify:
        case UnmapNotify:
          return 0;
      }
    }

    monitor_changes();

    yo22_paint();
    glXSwapBuffers(display, drawable);
  }

  glXMakeContextCurrent(display, 0, 0, 0);
  glXDestroyWindow(display, drawable);
  glXDestroyContext(display, context);
  XDestroyWindow(display, window);
  XCloseDisplay(display);

  return 0;
}
#endif // PLATFORM_X11

#ifdef PLATFORM_OSX
struct file_program_t {
  const char *filename;
  int fd;
  int program_counter;
  int updated;
};
static struct file_program_t files[Prog_COUNT];
static int kqfd;
static void monitor_changes() {
  int i;
  for (i = 0; i < Prog_COUNT; ++i) {
    struct kevent e;
    struct file_program_t *f = files + i;
    if (f->fd > 0) continue;
    
    f->fd = open(f->filename, 0);
    if (f->fd == -1) continue;
    e.ident = f->fd;
    e.filter = EVFILT_VNODE;
    e.flags = EV_ADD | EV_CLEAR;
    e.fflags = NOTE_DELETE | NOTE_WRITE | NOTE_EXTEND | NOTE_RENAME | NOTE_REVOKE;
    e.data = 0;
    e.udata = f;
    int result = kevent(kqfd, &e, 1, NULL, 0, NULL);
    if (result == -1) {
      close(f->fd);
      f->fd = -1;
      continue;
    }
    f->updated = 1;
    fprintf(stderr, "Listening for events on %s\n", f->filename);
  }
  
  for (;;) {
    struct kevent e;
    struct file_program_t *f;
    struct timespec ts = {0, 0};
    int result = kevent(kqfd, NULL, 0, &e, 1, &ts);
    CHECK(result != -1, "kevent");
    if (result == 0) break;
    
    f = e.udata;
    CHECK(f < files + Prog_COUNT, "kevent.udata");
    if (e.fflags & (NOTE_DELETE | NOTE_RENAME | NOTE_REVOKE)) {
      close(f->fd); /* will also EV_DELETE the kevent */
      f->fd = -1;
      continue;
    }
    
    f->updated = 1;
  }
  
  for (i = 0; i < Prog_COUNT; ++i)
    if (files[i].updated > 0) {
      char *src = read_file(files[i].filename);
      files[i].updated = 0;
      if (!src) continue;
      if (fragment_shaders[i] != NULL)
        free((char*)fragment_shaders[i]);
      fragment_shaders[i] = src;
      fprintf(stderr, "loading program %d ... ", i);
      if (0 == create_and_compile_program(i)) continue;
      if (program_counter > files[i].program_counter)
        program_counter = files[i].program_counter;
    }
}

static void glut_display() {
    monitor_changes();
    yo22_paint();
    glutSwapBuffers();
    glutPostRedisplay();
}

int main(int argc, char *argv[]) {
  glutInit(&argc, argv);
  
  files[ProgTerrainGenerate].filename = "generator.glsl";/*argv[1];*/
  files[ProgTerrainGenerate].program_counter = 0;
  files[ProgTerrainGenerate].updated = 1;
  files[ProgTerrainGenerate].fd = -1;
  files[ProgTerrainErode].filename = "eroder.glsl";/*argv[2];*/
  files[ProgTerrainErode].program_counter = 0;
  files[ProgTerrainErode].updated = 1;
  files[ProgTerrainErode].fd = -1;
  files[ProgTrace].filename = "tracer.glsl";/*argv[3];*/
  files[ProgTrace].program_counter = TERRAIN_ITERATIONS;
  files[ProgTrace].updated = 1;
  files[ProgTrace].fd = -1;
  files[ProgPostprocess].filename = "postprocessor.glsl"/*argv[4]*/;
  files[ProgPostprocess].program_counter = TOTAL_ITERATIONS;
  files[ProgPostprocess].updated = 1;
  files[ProgPostprocess].fd = -1;

  kqfd = kqueue();
  CHECK(kqfd != -1, "kqueue");
  
  glutCreateWindow("yo22");
  glutReshapeWindow(WIDTH, HEIGHT);
  yo22_init();
  glutDisplayFunc(glut_display);
  glutReshapeFunc(yo22_size);
  //glutKeyboardFunc(void (*func)(unsigned char key, int x, int y));
  glutMainLoop();
}
#endif

#ifdef PLATFORM_WINDOWS
static PIXELFORMATDESCRIPTOR pfd = {sizeof(pfd), 0, PFD_DOUBLEBUFFER | PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL,
  PFD_TYPE_RGBA, 24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
#ifdef __cplusplus
#error no
#endif
#ifdef DEBUG
int WINAPI WinMain(HINSTANCE a, HINSTANCE b, LPCWSTR c, int d) {
#else
void WinMainCRTStartup() {
#endif
  HWND w;
  HDC dc;
  RECT r;
  MSG m;
  const char *fnn = gl_names;
  PROC* fn = &gl;

  w = CreateWindowExA(0,"static",0,0x91000000,0,0,0,0,0,0,0,0);
  dc = GetDC(w);
  GetWindowRect(w, &r);
  SetPixelFormat(dc, ChoosePixelFormat(dc, &pfd), &pfd);
  wglMakeCurrent(dc, wglCreateContext(dc));

  for (; gl_names[0] != 0; ++fn) {
    *fn = wglGetProcAddress(gl_names);
     while (gl_names[0] != 0) ++gl_names;
     ++gl_names;
  }

  yo22_init();
  yo22_size(r.right - r.left, r.bottom - r.top);

  //ShowWindow(w, SW_SHOW);
  ShowCursor(FALSE);

  while (!GetAsyncKeyState(VK_ESCAPE)) {
    while (0 != PeekMessage(&m, w, 0, 0, PM_NOREMOVE)) {
      if (0 == GetMessage(&m, NULL, 0, 0)) goto exit;
      //TranslateMessage(&m);
      DispatchMessage(&m);
    }

    yo22_paint();
    SwapBuffers(dc);
  }

exit:
  ExitProcess(0);
}
#endif
