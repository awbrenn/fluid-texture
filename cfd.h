//
// Created by awbrenn on 1/20/16.
//

#ifndef CFD_H
#define CFD_H

class cfd
{
  public:
    // constructors/destructors
    cfd(const int nx, const int ny, const float dx, const float dt, int Nloops, int Oploops);
    ~cfd();

    // public methods
    void advect();
    void sources();

    // getters
    float* getColorPointer()    const { return color1; }

    // setters
    void setDensitySourceField(float* dsrc)     { densitySourceField = dsrc; }
    void setColorSourceField(float* csrc)       { colorSourceField = csrc; }
    void setObstructionSourceField(float* osrc) { obstructionSourceField = osrc; }
    void setDivergenceSourceField(float* dsrc) { divergenceSourceField = dsrc; }

    // indexing
    int dIndex(int i, int j)        const { return i+Nx*j; }
    int pIndex(int i, int j)        const { return i+Nx*j; }
    int oIndex(int i, int j)        const { return i+Nx*j; }
    int vIndex(int i, int j, int c) const { return (i+Nx*j)*2+c; }
    int cIndex(int i, int j, int c) const { return (i+Nx*j)*3+c; }

  private:
    int     Nx, Ny;
    int     nloops; // number of loops for pressure calculation
    int     oploops; // number of orthogonal projection loops
    float   Dx;
    float   dt;
    float   gravityX, gravityY;
    float   *density1, *density2;
    float   *velocity1, *velocity2;
    float   *color1, *color2;
    float   *divergence;
    float   *pressure;
    float   *obstruction;
    float   *densitySourceField;
    float   *colorSourceField;
    float   *obstructionSourceField;
    float   *divergenceSourceField;

    // private methods
    void addSourceColor();
    void addSourceDensity();
    void addSourceObstruction();
    void computeDivergence();
    void computePressure();
    void computePressureForces(int i, int j, float* force_x, float* force_y);
    void computeVelocityBasedOnPressureForces();
    void bilinearlyInterpolate(const int ii, const int jj, const float x, const float y);
    void computeVelocity(float force_x, float force_y);
    void computeObstructedFields();
    const float InterpolateColor(int i, int j, int c, float w1, float w2, float w3, float w4);
    const float InterpolateVelocity(int i, int j, int c, float w1, float w2, float w3, float w4);
    const float InterpolateDensity(int i, int j, float w1, float w2, float w3, float w4);
    const float getDensity(int i, int j);
    const float getVelocity(int i, int j, int c);
    const float getColor(int i, int j, int c);
    const float getPressure(int i, int j);
    const float getDivergence(int i, int j);
};

#endif //CFD_H