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




#include <cmath>
#include "CmdLineFind.h"

#include <GL/gl.h>   // OpenGL itself.
#include <GL/glu.h>  // GLU support library.
#include <GL/glut.h> // GLUT support library.
#include "cfd.h"

#include <iostream>
#include <OpenImageIO/imageio.h>
#include <omp.h>


using namespace std;
using namespace lux;
OIIO_NAMESPACE_USING

int iwidth, iheight;
float* display_map;
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


void writeImage()
{
  char buffer [256];
  if (sprintf(buffer, "%sfluid_simulator_%04d.jpg", output_path.c_str(), frame_count++) < 0) {
    handleError((const char *) "creating filename in writeImage() failed", 0);
    return;
  }
  const char *filename = buffer;
  const int xres = iwidth, yres = iheight;
  const int channels = 3; // RGB
  float* pixels = new float[xres*yres*channels];
  ImageOutput *out = ImageOutput::create (filename);
  if (! out) {
    handleError((const char *) "creating output file in writeImage() failed", 0);
    return;
  }

  long index = 0;
  float* current_fluid_color = fluid->getColorPointer();

  for( int j=0;j<iheight;j++)
  {
    for (int i = 0; i < iwidth; i++)
    {
      for (int c = 0; c < channels; c++)
      {
        pixels[ (i + iwidth*(iheight - j - 1))*channels + c ] = current_fluid_color[index++];
      }
    }
  }

  ImageSpec spec (xres, yres, channels, TypeDesc::FLOAT);
  out->open (filename, spec);
  out->write_image (TypeDesc::FLOAT, pixels);
  out->close ();
  delete out;

  delete pixels;
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
      source_brush[ii][jj] = pow(radius,0.5);
      obstruction_brush[ii][jj] = (float)(1.0 - pow(radius, 1.0/4.0));
    }
  }
}


void setNbCores( int nb )
{
  omp_set_num_threads( nb );
}


//----------------------------------------------------
//
//  Painting and Display Code
//
//----------------------------------------------------


void ConvertToDisplay()
{
  float *color = fluid->getColorPointer();
  for( int j=0;j<iheight;j++ )
  {
#pragma omp parallel for
    for(int i=0;i<iwidth;i++ )
    {
      int index = i + iwidth*j;
      float r,g,b;
      r = color[index*3];
      g = color[index*3+1];
      b = color[index*3+2];
      display_map[3*index  ] = r * scaling_factor;
      display_map[3*index+1] = g * scaling_factor;
      display_map[3*index+2] = b * scaling_factor;
    }
  }
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
        color_source[3 * index + 2] += source_brush[ix - xstart][iy - ystart];
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
        color_source[3 * index] += source_brush[ix - xstart][iy - ystart];
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


int main(int argc, char** argv)
{
  CmdLineFind clf(argc, argv);

  iwidth = clf.find("-NX", 512, "Horizontal grid points");
  iheight = clf.find("-NY", iwidth, "Vertical grid points");

  int nloops = clf.find("-nloops", 3, "Number of loops over pressure.");
  int oploops = clf.find("-oploops", 1, "Number of orthogonal projection loops.");

  output_path = clf.find("-output_path", "output_images/", "Output path for writing image sequence");

  setNbCores(4);

  string imagename = clf.find("-image", "dali1.jpeg", "Image to drive color");

  clf.usage("-h");
  clf.printFinds();
  PrintUsage();
  cout << "\n\nPROGRAM OUTPUT:\n";

  // initialize a few variables
  scaling_factor = 1.0;
  toggle_animation_on_off = true;
  capture_mode = false;

  // if reading the image fails we need to allocate space for color_source
  if (readOIIOImage(imagename.c_str()) != 0)
    color_source = new float[iwidth*iheight*3]();

  density_source = new float[iwidth*iheight]();

  // create obstruction source and initialize it to 1.0
  obstruction_source = new float[iwidth*iheight];
  for(int i=0;i<iwidth*iheight;i++ ) { obstruction_source[i] = 1.0; }

  divergance_source = new float[iwidth*iheight*3]();

  // initialize fluid
  fluid = new cfd(iwidth, iheight, 1.0, (float)(1.0/24.0), nloops, oploops);
  fluid->setColorSourceField(color_source);

  display_map = new float[iwidth*iheight*3];

  InitializeBrushes(BRUSH_SIZE);

  paint_mode = PAINT_SOURCE;

  // GLUT routines
  glutInit(&argc, argv);

  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
  glutInitWindowSize(iwidth, iheight);

  // Open a window 
  char title[] = "Fluid Simulator";
  glutCreateWindow(title);

  glClearColor(1, 1, 1, 1);

  glutDisplayFunc(&cbDisplay);
  glutIdleFunc(&cbIdle);
  glutKeyboardFunc(&cbOnKeyboard);
  glutMouseFunc(&cbMouseDown);
  glutMotionFunc(&cbMouseMove);

  glutMainLoop();
  return 1;
}