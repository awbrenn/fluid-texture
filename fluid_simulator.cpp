//------------------------------------------------
//
//  Program: fluid_simulator
//  Author:  Austin Brennan
//  Course:  Realtime Fluid Simulation (CPSC 8810)
//  School:  Clemson University
//
//-------------------------------------------------

//-------------------------------------------------
//
//  usage:
//
//  fluid_simulator is an interactive paint program
//  in which the user paints density, color, and
//  or divergence sources that flow using
//  computational fluid dynamics and react with 
//  obstructions in the space.
//
//  There are several paint modes.  Typing 'o' puts the
//  program in obstruction painting mode. When you
//  hold down the left mouse button and paint, you
//  will see a black obstruction painted.  This 
//  obstruction may be any shape.
//
//  Typing 's' puts the program in source painting 
//  mode.  Now painting with the left mouse button
//  down injects density into the simulation.
//  The flow it produces evolves as you
//  continue to paint.  The flow bounces off any
//  obstructions that have been painted or are
//  subsequently painted.
//
//  Typing 'b' puts the program in painting positive
//  divergence mode. Similiarly typing 'r' puts the
//  program in painting negative divergence mode.
//  Painting in this mode injects divergence into the
//  simulation.
//
//  Typing ',' or '.' increases or decreases the brush
//  size respectively.
//
//  Typing '=' and '-' brightens and darkens the display.
//
//  Pressing the spacebar starts and stops the flow 
//  evolution.
//
//
//-------------------------------------------------
#define GL_GLEXT_PROTOTYPES 1
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>
#include <cmath>
#include "CmdLineFind.h"
#include <stdio.h>
#include <unistd.h>

#include "cfd.h"

#ifdef __APPLE__
  #include <OpenGL/gl.h>   // OpenGL itself.
  #include <OpenGL/glu.h>  // GLU support library.
  #include <GLUT/glut.h> // GLUT support library.
#else
  #include <GL/gl.h>   // OpenGL itself.
  #include <GL/glu.h>  // GLU support library.
  #include <GL/glut.h> // GLUT support library.
  #include <omp.h>
#endif

#include <OpenImageIO/imageio.h>



using namespace std;
using namespace lux;
OIIO_NAMESPACE_USING

int iwidth, iheight;
unsigned int shader_program;
unsigned char* display_map;
float* density_source;
float* color_source;
float* obstruction_source;
float* divergance_source;
cfd *fluid;
int frame_count = 0;
string output_path;
bool capture_mode;

int paint_mode;
enum{ PAINT_OBSTRUCTION, PAINT_SOURCE, PAINT_DIVERGENCE_POSITIVE, PAINT_DIVERGENCE_NEGATIVE, PAINT_COLOR };

bool toggle_animation_on_off;

float scaling_factor;

int BRUSH_SIZE = 11;
float **obstruction_brush = NULL;
float **source_brush = NULL;

int xmouse_prev, ymouse_prev;

////////  OpenImageIO reader


void handleError(const char* error_message, int kill)
{
  fprintf(stderr, "Error: %s\n\n", error_message);

  if (kill == 1)
    exit(-1);
}


//----------------------------------------------------
//
//  Read and Write Images
//
//----------------------------------------------------


int readOIIOImage( const char* fname)
{
  int channels;
  ImageInput *in = ImageInput::create (fname);
  if (! in) { return -1; }
  ImageSpec spec;
  in->open (fname, spec);
  iwidth = spec.width; // note iwidth and iheight are set to to image size
  iheight = spec.height;
  channels = spec.nchannels;
  float* pixels = new float[iwidth*iheight*channels];
  color_source = new float[iwidth*iheight*channels]; // allocate appropriate space for image

  in->read_image (TypeDesc::FLOAT, pixels);
  long index = 0;
  for( int j=0;j<iheight;j++)
  {
    for( int i=0;i<iwidth;i++ )
    {
      for( int c=0;c<channels;c++ )
      {
        color_source[ (i + iwidth*(iheight - j - 1))*channels + c ] = pixels[index++];
      }
    }
  }
  delete pixels;

  in->close ();
  delete in;

  return 0;
}


void writeImage() {
  char buffer[256];

  if (sprintf(buffer, "%sfluid_simulator_%04d.jpg", output_path.c_str(), frame_count++) < 0) {
    handleError((const char *) "creating filename in writeImage() failed", 0);
    return;
  }
  const char *filename = buffer;
  const unsigned int channels = 3; // RGB
  float *write_pixels = new float[1024 * 1024 * channels];
  float *window_pixels = new float[1024 * 1024 * channels];
  ImageOutput *out = ImageOutput::create(filename);
  if (!out) {
    handleError((const char *) "creating output file in writeImage() failed", 0);
    return;
  }

  glReadPixels(0, 0, 1024, 1024, GL_RGB, GL_FLOAT, window_pixels);
  long index = 0;
  for (unsigned int j = 0; j < 1024; j++) {
    for (unsigned int i = 0; i < 1024; i++) {
      for (unsigned int c = 0; c < channels; c++) {
        write_pixels[(i + 1024 * (1024 - j - 1)) * channels + c] = window_pixels[index++]; //color[index++];
      }
    }
  }

  ImageSpec spec (1024, 1024, channels, TypeDesc::FLOAT);
  out->open (filename, spec);
  out->write_image (TypeDesc::FLOAT, write_pixels);
  out->close ();
  delete out;
  delete write_pixels;
  delete window_pixels;
}


//----------------------------------------------------
//
//  Initialize brushes and set number of cores
//
//----------------------------------------------------


void InitializeBrushes(int new_brush_size)
{
  // deallocate old brushes
  if (source_brush != NULL) {
    for (int i = 0; i < BRUSH_SIZE; ++i) {delete source_brush[i]; }
    delete source_brush;
  }
  if (obstruction_brush != NULL) {
    for (int i = 0; i < BRUSH_SIZE; ++i) {delete obstruction_brush[i]; }
    delete obstruction_brush;
  }

  // set BRUSH_SIZE to the new brush size. clamp min size to 3
  if (new_brush_size < 3)
    BRUSH_SIZE = 3;
  else
    BRUSH_SIZE = new_brush_size;

  // allocate memory for new brushes
  source_brush = new float*[BRUSH_SIZE];
  for (int i = 0; i < BRUSH_SIZE; ++i) {source_brush[i] = new float[BRUSH_SIZE]; }
  obstruction_brush = new float*[BRUSH_SIZE];
  for (int i = 0; i < BRUSH_SIZE; ++i) {obstruction_brush[i] = new float[BRUSH_SIZE]; }

  int brush_width = (BRUSH_SIZE-1)/2;
  for( int j=-brush_width;j<=brush_width;j++ )
  {
    int jj = j + brush_width;
    float jfactor =  (float(brush_width) - (float)fabs(j) )/float(brush_width);
    for( int i=-brush_width;i<=brush_width;i++ )
    {
      int ii = i + brush_width;
      float ifactor =  (float(brush_width) - (float)fabs(i) )/float(brush_width);
      float radius = (float) ((jfactor * jfactor + ifactor * ifactor) / 2.0);
      source_brush[ii][jj] = powf(radius,0.5);
      obstruction_brush[ii][jj] = (float)(1.0 - pow(radius, 1.0/4.0));
    }
  }
}

#ifdef __linux__
void setNbCores( int nb )
{
  omp_set_num_threads( nb );
}
#endif


//----------------------------------------------------
//
//  Painting and Display Code
//
//----------------------------------------------------


void ConvertToDisplay()
{
  float *color = fluid->getColorPointer();
#ifdef __linux__
#pragma omp parallel for
#endif
  for (int i = 0; i < iwidth*iheight*3; ++i) { display_map[i] = (unsigned char)(color[i] * 255.0f); }
}

void resetScaleFactor( float amount )
{
   scaling_factor *= amount;
}


void DabSomePaint( int x, int y ) {
  float divergence_source_magnitude = 250.0f;
  int brush_width = (BRUSH_SIZE - 1) / 2;
  int xstart = x - brush_width;
  int ystart = y - brush_width;
  if (xstart < 0) { xstart = 0; }
  if (ystart < 0) { ystart = 0; }

  int xend = x + brush_width;
  int yend = y + brush_width;
  if (xend >= iwidth) { xend = iwidth - 1; }
  if (yend >= iheight) { yend = iheight - 1; }

  if (paint_mode == PAINT_OBSTRUCTION) {
    for (int ix = xstart; ix <= xend; ix++) {
      for (int iy = ystart; iy <= yend; iy++) {
        int index = ix + iwidth * (iheight - iy - 1);
        obstruction_source[index] *= obstruction_brush[ix - xstart][iy - ystart];
      }
    }
    fluid->setObstructionSourceField(obstruction_source);
  }
  else if (paint_mode == PAINT_SOURCE) {
    for (int ix = xstart; ix <= xend; ix++) {
      for (int iy = ystart; iy <= yend; iy++) {
        int index = ix + iwidth * (iheight - iy - 1);
        color_source[3 * index] += source_brush[ix - xstart][iy - ystart];
        color_source[3 * index + 1] += source_brush[ix - xstart][iy - ystart];
        color_source[3 * index + 2] += source_brush[ix - xstart][iy - ystart];
        density_source[index] += source_brush[ix - xstart][iy - ystart];
      }
    }
    fluid->setColorSourceField(color_source);
    fluid->setDensitySourceField(density_source);
  }
  else if (paint_mode == PAINT_DIVERGENCE_POSITIVE ) {
    for (int ix = xstart; ix <= xend; ix++) {
      for (int iy = ystart; iy <= yend; iy++) {
        int index = ix + iwidth * (iheight - iy - 1);
        //color_source[3 * index + 2] += source_brush[ix - xstart][iy - ystart];
        divergance_source[index] += source_brush[ix - xstart][iy - ystart]*divergence_source_magnitude;
      }
    }
    fluid->setColorSourceField(color_source);
    fluid->setDivergenceSourceField(divergance_source);
  }
  else if ( paint_mode == PAINT_DIVERGENCE_NEGATIVE ) {
    for (int ix = xstart; ix <= xend; ix++) {
      for (int iy = ystart; iy <= yend; iy++) {
        int index = ix + iwidth * (iheight - iy - 1);
        //color_source[3 * index] += source_brush[ix - xstart][iy - ystart];
        divergance_source[index] += source_brush[ix - xstart][iy - ystart]*divergence_source_magnitude*(-1.0f);
      }
    }
    fluid->setColorSourceField(color_source);
    fluid->setDivergenceSourceField(divergance_source);
  }

  return;
}


//----------------------------------------------------
//
//  GL and GLUT callbacks
//
//----------------------------------------------------


void cbDisplay( void )
{
  glClear(GL_COLOR_BUFFER_BIT );
  glDrawPixels( iwidth, iheight, GL_RGB, GL_FLOAT, display_map );
  glutSwapBuffers();
}


void update()
{
  fluid->advect();
  fluid->sources();
}

// animate and display new result
void cbIdle()
{
  if (toggle_animation_on_off)
    update();
  ConvertToDisplay();
  if (capture_mode)
    writeImage();
  glutPostRedisplay(); 
}


void cbOnKeyboard( unsigned char key, int x, int y )
{
  switch (key) 
  {
    case '-': case '_':
      resetScaleFactor( 0.9 );
      break;

    case '+': case '=':
      resetScaleFactor( (float)(1.0/0.9) );
      break;

    case 'c':
      scaling_factor = 1.0;
      break;

    case ' ':
      toggle_animation_on_off = !toggle_animation_on_off;
      if (toggle_animation_on_off)
        cout << "Animation Toggled On" << endl;
      else
        cout << "Animation Toggled Off" << endl;
      break;

    case ',' : case '<':
      InitializeBrushes(BRUSH_SIZE-2);
      cout << "Setting Brush Size To " << BRUSH_SIZE << endl;
      break;

    case '.': case '>':
      InitializeBrushes(BRUSH_SIZE+2);
      cout << "Setting Brush Size To " << BRUSH_SIZE << endl;
      break;

    case 'o':
      paint_mode = PAINT_OBSTRUCTION;
      cout << "Paint Obstruction Mode" << endl;
      break;

    case 's':
      paint_mode = PAINT_SOURCE;
      cout << "Paint Source Density Mode" << endl;
      break;

    case 'b':
      paint_mode = PAINT_DIVERGENCE_POSITIVE;
      cout << "Paint Positive Divergence Mode" << endl;
      break;

    case 'r':
      paint_mode = PAINT_DIVERGENCE_NEGATIVE;
      cout << "Paint Negative Divergence Mode" << endl;
      break;

    case 'w':
      capture_mode = !capture_mode;
      if (capture_mode)
        cout << "Starting Capture..." << endl;
      else
        cout << "...Ending Capture" << endl;
      break;

    case 'q':
      cout << "Exiting Program" << endl;
      exit(0);

    default:
    break;
  }
}


void cbMouseDown( int button, int state, int x, int y )
{
  if( button != GLUT_LEFT_BUTTON ) { return; }
  if( state != GLUT_DOWN ) { return; }
  xmouse_prev = x;
  ymouse_prev = y;
  DabSomePaint( x, y );
}


void cbMouseMove( int x, int y )
{
  xmouse_prev = x;
  ymouse_prev = y;
  DabSomePaint( x, y ); 
}


//----------------------------------------------------
//
//  Printing Usage
//
//----------------------------------------------------


void PrintUsage()
{
  cout << "fluid_simulator keyboard choices\n";
  cout << "s        turns on painting source strength\n";
  cout << "o        turns on painting obstructions\n";
  cout << "b        turns on painting positive divergence\n";
  cout << "r        turns on painting negative divergence\n";
  cout << "+/-      increase/decrease brightness of display\n";
  cout << ",/.      increase/decrease brush size\n";
  cout << "c        clears changes to brightness\n";
  cout << "w        starts capture mode. file path can be set with -output_path flag\n";
  cout << "spacebar paused the simulation. pressing it again un-pauses the simulation\n";
  cout << "q        exits the program\n";
}


//----------------------------------------------------
//
// Main
//
//----------------------------------------------------

struct point {
    float x, y, z;
};


void set_texture() {

  glBindTexture(GL_TEXTURE_2D,1);
  glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,iwidth,iheight,0,GL_RGB,
               GL_UNSIGNED_BYTE,display_map);
  glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
  glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
}

void drawStuff() {
  int i;
  struct point tile0[4]={{0.0,0.0,0.0},{0.5,0.0,0.0},{0.5,0.5,0.0},{0.0,0.5,0.0}};
  struct point tile1[4]={{0.0,0.0,0.1},{0.75,0.0,0.1},{0.75,0.75,0.1},{0.0,0.75,0.1}};
  struct point tile2[4]={{-0.25,0.0,-0.1},{0.25,0.0,-0.1},{0.25,0.5,-0.1},{-0.25,0.5,-0.1}};
  struct point tile3[4]={{0.0,0.0,0.5},{0.75,0.0,0.5},{0.75,0.75,0.5},{0.0,0.75,0.5}};
  struct point tile4[4]={{0.5,0.0,0.6},{1.0,0.0,0.6},{1.0,1.0,0.6},{0.5,1.0,0.6}};
  struct point tile5[4]={{0.5,0.0,-0.2},{1.0,0.0,-0.2},{1.0,1.0,-0.2},{0.5,1.0,-0.2}};
  struct point tile6[4]={{0.2,0.0,-0.3},{0.95,0.0,-0.3},{0.95,0.75,-0.3},{0.2,0.75,-0.3}};
  struct point tile7[4]={{0.0,0.0,-0.4},{0.5,0.0,-0.4},{0.5,0.5,-0.4},{0.0,0.5,-0.4}};
  struct point tile8[4]={{0.1,0.0,-0.41},{0.6,0.0,-0.41},{0.6,0.5,-0.41},{0.1,0.5,-0.41}};
  float mytexcoords[4][2] = {{0.0,0.0},{1.0,0.0},{1.0,1.0},{0.0,1.0}};

  set_texture();
  glClearColor(0.0,0.0,0.0,0.0);
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

  glUseProgram((GLuint) shader_program);		// THIS IS IT!
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D,1);
  glEnable(GL_TEXTURE_2D);
  glBegin(GL_QUADS);
  glNormal3f(0.0,0.0,1.0);

  for(i=0;i<4;i++) {
    glTexCoord2fv(mytexcoords[i]);
    glVertex3f(tile8[i].x,tile8[i].y,tile8[i].z);
  }

  for(i=0;i<4;i++) {
    glTexCoord2fv(mytexcoords[i]);
    glVertex3f(tile7[i].x,tile7[i].y,tile7[i].z);
  }

  for(i=0;i<4;i++) {
    glTexCoord2fv(mytexcoords[i]);
    glVertex3f(tile6[i].x,tile6[i].y,tile6[i].z);
  }

  for(i=0;i<4;i++) {
    glTexCoord2fv(mytexcoords[i]);
    glVertex3f(tile5[i].x,tile5[i].y,tile5[i].z);
  }

  for(i=0;i<4;i++) {
    glTexCoord2fv(mytexcoords[i]);
    glVertex3f(tile2[i].x,tile2[i].y,tile2[i].z);
  }

  for(i=0;i<4;i++) {
    glTexCoord2fv(mytexcoords[i]);
    glVertex3f(tile0[i].x,tile0[i].y,tile0[i].z);
  }

  for(i=0;i<4;i++) {
    glTexCoord2fv(mytexcoords[i]);
    glVertex3f(tile1[i].x,tile1[i].y,tile1[i].z);
  }

  for(i=0;i<4;i++) {
    glTexCoord2fv(mytexcoords[i]);
    glVertex3f(tile3[i].x,tile3[i].y,tile3[i].z);
  }

  for(i=0;i<4;i++) {
    glTexCoord2fv(mytexcoords[i]);
    glVertex3f(tile4[i].x,tile4[i].y,tile4[i].z);
  }

  glEnd();
  glFlush();
}

void setupViewVolume()
{
  struct point eye, view, up;

// specify size and shape of view volume
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45.0,1.0,0.1,20.0);

// specify position for view volume
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  eye.x = 0.5; eye.y = 0.5; eye.z = 2.0;
  view.x = 0.5; view.y = 0.5; view.z = 0.0;
  up.x = 0.0; up.y = 1.0; up.z = 0.0;

  gluLookAt(eye.x,eye.y,eye.z,view.x,view.y,view.z,up.x,up.y,up.z);
}


char *read_shader_program(char *filename)
{
  FILE *fp;
  char *content = NULL;
  int fd, count;
  fd = open(filename,O_RDONLY);
  count = (int) lseek(fd,0,SEEK_END);
  close(fd);
  content = (char *)calloc(1,(size_t)(count+1));
  fp = fopen(filename,"r");
  count = fread(content,sizeof(char),count,fp);
  content[count] = '\0';
  fclose(fp);
  return content;
}


unsigned int setShaders()
{
  GLint vertCompiled, fragCompiled;
  char *vs, *fs;
  GLuint v, f, p;

  v = glCreateShader(GL_VERTEX_SHADER);
  f = glCreateShader(GL_FRAGMENT_SHADER);
  vs = read_shader_program((char *) "/home/awbrenn/Documents/workspace/fluid2D/midterm_show/sim_tex.vert");
  fs = read_shader_program((char *) "/home/awbrenn/Documents/workspace/fluid2D/midterm_show/sim_tex.frag");
  glShaderSource(v,1,(const char **)&vs,NULL);
  glShaderSource(f,1,(const char **)&fs,NULL);
  free(vs);
  free(fs);
  glCompileShader(v);
  glCompileShader(f);
  p = glCreateProgram();
  glAttachShader(p,f);
  glAttachShader(p,v);
  glLinkProgram(p);
  return(p);
}


void set_uniform_parameters(unsigned int p)
{
  int location;
  location = glGetUniformLocation(p,"mytexture");
  glUniform1i(location,0);
}

void lights()
{
  float light0_ambient[] = { 0.0, 0.0, 0.0, 0.0 };
  float light0_diffuse[] = { 1.0, 1.0, 1.0, 0.0 };
  float light0_specular[] = { 1.0, 1.0, 1.0, 0.0 };
  float light0_position[] = { M_SQRT2, 2.0, 2.0, 1.0 };
  float light0_direction[] = { -M_SQRT2, -2.0, -2.0, 1.0};

  glLightModelfv(GL_LIGHT_MODEL_AMBIENT,light0_ambient);
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,1);
  glLightfv(GL_LIGHT0,GL_AMBIENT,light0_ambient);
  glLightfv(GL_LIGHT0,GL_DIFFUSE,light0_diffuse);
  glLightfv(GL_LIGHT0,GL_SPECULAR,light0_specular);
  glLightf(GL_LIGHT0,GL_SPOT_EXPONENT,0.0);
  glLightf(GL_LIGHT0,GL_SPOT_CUTOFF,180.0);
  glLightf(GL_LIGHT0,GL_CONSTANT_ATTENUATION,1.0);
  glLightf(GL_LIGHT0,GL_LINEAR_ATTENUATION,0.0);
  glLightf(GL_LIGHT0,GL_QUADRATIC_ATTENUATION,0.0);
  glLightfv(GL_LIGHT0,GL_POSITION,light0_position);
  glLightfv(GL_LIGHT0,GL_SPOT_DIRECTION,light0_direction);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
}


void material()
{
  float mat_diffuse[] = {1.0,1.0,0.0,1.0};
  float mat_specular[] = {1.0,1.0,1.0,1.0};
  float mat_shininess[] = {2.0};

  glMaterialfv(GL_FRONT,GL_DIFFUSE,mat_diffuse);
  glMaterialfv(GL_FRONT,GL_SPECULAR,mat_specular);
  glMaterialfv(GL_FRONT,GL_SHININESS,mat_shininess);
}


int main(int argc, char** argv)
{
  CmdLineFind clf(argc, argv);

  iwidth = clf.find("-NX", 512, "Horizontal grid points");
  iheight = clf.find("-NY", iwidth, "Vertical grid points");

  int nloops = clf.find("-nloops", 3, "Number of loops over pressure.");
  int oploops = clf.find("-oploops", 1, "Number of orthogonal projection loops.");

  output_path = clf.find("-output_path", "output_images/", "Output path for writing image sequence");

#ifdef __linux__
  setNbCores(4);
#endif

//  string imagename = clf.find("-image", "dali1.jpeg", "Image to drive color");

  clf.usage("-h");
  clf.printFinds();
  PrintUsage();
  cout << "\n\nPROGRAM OUTPUT:\n";

  // initialize a few variables
  scaling_factor = 1.0;
  toggle_animation_on_off = true;
  capture_mode = true;

//  // if reading the image fails we need to allocate space for color_source
//  if (readOIIOImage(imagename.c_str()) != 0)
//    color_source = new float[iwidth*iheight*3]();
  iwidth = 128;
  iheight = 128;

  color_source = new float[iwidth*iheight*3]();

  density_source = new float[iwidth*iheight]();

  // create obstruction source and initialize it to 1.0
  obstruction_source = new float[iwidth*iheight];
  for(int i=0;i<iwidth*iheight;i++ ) { obstruction_source[i] = 1.0; }

  divergance_source = new float[iwidth*iheight*3]();

  display_map = new unsigned char[iwidth*iheight*3];

  // initialize fluid
  fluid = new cfd(iwidth, iheight, 1.0, (float)(1.0/24.0), nloops, oploops);
  fluid->setColorSourceField(color_source);
  update();
  ConvertToDisplay();

  InitializeBrushes(BRUSH_SIZE);

  paint_mode = PAINT_SOURCE;

  DabSomePaint(64, 64);

  paint_mode = PAINT_DIVERGENCE_NEGATIVE;
  DabSomePaint(60, 60);
  DabSomePaint(30, 30);
  DabSomePaint(70, 70);
  DabSomePaint(100, 100);
  DabSomePaint(64, 64);
  DabSomePaint(64, 64);
  DabSomePaint(64, 64);
  DabSomePaint(64, 64);

  // GLUT routines
  glutInit(&argc, argv);

  glutInitDisplayMode(GLUT_RGBA| GLUT_MULTISAMPLE);
  glutInitWindowPosition(700, 300);
  glutInitWindowSize(1024, 1024);

  // Open a window
  char title[] = "Fluid Simulator";
  glutCreateWindow(title);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_MULTISAMPLE_ARB);
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  set_texture();
  setupViewVolume();
  lights();
  material();
  shader_program = setShaders();
  set_uniform_parameters(shader_program);

  glutDisplayFunc(drawStuff);
  glutIdleFunc(&cbIdle);
  glutKeyboardFunc(cbOnKeyboard);

  cout << glGetString(GL_VERSION) << endl;
//  glutMouseFunc(&cbMouseDown);
//  glutMotionFunc(&cbMouseMove);

  glutMainLoop();
  return 1;
}