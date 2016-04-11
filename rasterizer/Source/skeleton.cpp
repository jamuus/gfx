#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using namespace glm;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 400;
const int SCREEN_HEIGHT = 300;

#define PI 3.14159265358979323846f
SDL_Surface* screen;
int t;
vector<Triangle> triangles;
vec3 currentColor;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];


vec3 lightPos(0, -0.5, -0.7);
vec3 lightPower = 12.f * vec3( 1, 1, 1 );
vec3 indirectLightPowerPerArea = 0.1f * vec3( 1, 1, 1 );

float focalLength = SCREEN_HEIGHT;
vec3 cameraPosition(0, 0, -4);

vec3 currentNormal;
vec3 currentReflectance;

float yaw = 0.0f;
float pitch = 0.0f;

float translateSpeed = 2.0f;
float rotateSpeed = 0.1f;
float lightTranslateDist = 1.0f;
mat3 R;

vec3 image[SCREEN_WIDTH][SCREEN_HEIGHT];

int dx;
int dy;

struct Pixel {
    int x;
    int y;
    float zinv;
    vec3 pos3d;

    inline Pixel operator-(Pixel a)
    {
        return (Pixel) {x - a.x, y - a.y, zinv - a.zinv};
    }
};

struct Vertex {
    vec3 position;
    // vec3 normal;
    // vec3 reflectance;
};
/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result );
void Interpolate( Pixel a, Pixel b, vector<Pixel>& result );

void ComputePolygonRows(
    const vector<Pixel>& vertexPixels,
    vector<Pixel>& leftPixels,
    vector<Pixel>& rightPixels );

void VertexShader( const vec3& v, Pixel& p );



int main( int argc, char* argv[] )
{
    screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
    t = SDL_GetTicks(); // Set start value for timer.
    LoadTestModel(triangles);

    SDL_WarpMouse(SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2);

    while ( NoQuitMessageSDL() ) {
        Update();
        Draw();
        for (int y = 0; y < SCREEN_HEIGHT; ++y)
            for (int x = 0; x < SCREEN_WIDTH; ++x)
                depthBuffer[y][x] = 0;
        fflush(stdout);
    }
    printf("\n");

    SDL_SaveBMP(screen, "screenshot.bmp");
    return 0;
}

float averageFps = 0.0f;

void Update()
{
    // Compute frame time:
    int t2 = SDL_GetTicks();
    float dt = float(t2 - t);
    t = t2;
    float dtsec = dt / 1000.0f;


    SDL_GetMouseState( &dx, &dy );
    SDL_WarpMouse(SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2);
    dx = dx - SCREEN_WIDTH / 2;
    dy = dy - SCREEN_HEIGHT / 2;

    Uint8* keystate = SDL_GetKeyState( 0 );
    // if ( keystate[SDLK_UP] ) {
    if ( keystate[SDLK_w] ) {
        // Move camera forward
        cameraPosition += R * (vec3(0, 0, 1) * dtsec * translateSpeed);
    }
    // if ( keystate[SDLK_DOWN] ) {
    if ( keystate[SDLK_s] ) {
        // Move camera backward
        cameraPosition -= R * (vec3(0, 0, 1) * dtsec * translateSpeed);
    }
    // if ( keystate[SDLK_LEFT] ) {
    if ( keystate[SDLK_a] ) {
        // Move camera to the left
        cameraPosition -= R * (vec3(1, 0, 0) * dtsec * translateSpeed);
        // rotate left
        // yaw -= dtsec * rotateSpeed;
    }
    // if ( keystate[SDLK_RIGHT] ) {
    if ( keystate[SDLK_d] ) {
        // Move camera to the right
        cameraPosition += R * (vec3(1, 0, 0) * dtsec * translateSpeed);
    }
    if ( keystate[SDLK_UP] ) {
        lightPos.z += lightTranslateDist * dtsec;
    }
    if ( keystate[SDLK_LEFT] ) {
        lightPos.x -= lightTranslateDist * dtsec;
    }
    if ( keystate[SDLK_DOWN] ) {
        lightPos.z -= lightTranslateDist * dtsec;
    }
    if ( keystate[SDLK_RIGHT] ) {
        lightPos.x += lightTranslateDist * dtsec;
    }
    if ( keystate[SDLK_q] ) {
        lightPos.y += lightTranslateDist * dtsec;
    }
    if ( keystate[SDLK_e] ) {
        lightPos.y -= lightTranslateDist * dtsec;
    }
    yaw -= dx * dtsec * rotateSpeed;
    pitch += dy * dtsec * rotateSpeed;

    if (yaw > PI * 2)
        yaw -= PI * 2;
    if (yaw < 0)
        yaw += PI * 2;

    if (pitch > PI * 2)
        pitch -= PI * 2;
    if (pitch < 0)
        pitch += PI * 2;

    R[0] = vec3( cos(yaw), 0, sin(yaw));
    R[1] = vec3( 0,        1, 0);
    R[2] = vec3(-sin(yaw), 0, cos(yaw));

    mat3 rpitch;
    rpitch[0] = vec3(1, 0, 0);
    rpitch[1] = vec3(0, cos(pitch), -sin(pitch));
    rpitch[2] = vec3(0, sin(pitch), cos(pitch));

    R = R * rpitch;

    averageFps = averageFps - (averageFps / 30) + ((1 / dtsec) / 30);

    // cout << "Render time: " << dt << " ms." << endl;
    printf("%.4f ms render, %.4f fps, %.4f avg, ", dt, (1 / dtsec), averageFps);
}

void PixelShader(const Pixel& p)
{
    int x = p.x;
    int y = p.y;
    if (x < 0 || x >= SCREEN_WIDTH || y < 0 || y >= SCREEN_HEIGHT) {
        // printf("(%d, %d)\n", x, y);
        return;
    }
    if ( p.zinv > depthBuffer[y][x] ) {
        depthBuffer[y][x] = p.zinv;

        vec3 nhat = currentNormal;
        vec3 r = lightPos - p.pos3d;
        vec3 rhat = r / length(r);
        vec3 t = lightPower * glm::max(dot(rhat, nhat), 0.0f);
        float rlength = length(rhat);
        vec3 D = t / (4 * PI * rlength * rlength);
        vec3 badR = currentReflectance * (D
                                          + indirectLightPowerPerArea
                                         );

        vec3 illumination = badR * currentColor
                            // wew
                            * R;

        image[x][y] = illumination;

        // PutPixelSDL(screen, x, y, illumination);
    }
}

void VertexShader( const Vertex& v, Pixel& p )
{
    vec3 newV = (v.position - cameraPosition) * R;
    p.zinv = 1.0f / newV.z;
    p.x = (int)(focalLength * newV.x * p.zinv) + SCREEN_WIDTH / 2;
    p.y = (int)(focalLength * newV.y * p.zinv) + SCREEN_HEIGHT / 2;

    p.pos3d = v.position;
    //lmao
    // * R;
}
void DrawRows( SDL_Surface* surface,
               const vector<Pixel>& leftPixels,
               const vector<Pixel>& rightPixels )
{
    for (int row = 0; row < leftPixels.size(); row++) {
        Pixel lp = leftPixels[row];
        Pixel rp = rightPixels[row];
        vector<Pixel> utput(rp.x - lp.x + 1);
        Interpolate(lp, rp, utput);
        for (int i = 0; i < utput.size(); ++i) {
            Pixel p = utput[i];
            PixelShader(p);
        }
    }
}

void DrawPolygon( SDL_Surface* surface,
                  const vector<Vertex>& vertices )
{
    int V = vertices.size();
    vector<Pixel> vertexPixels( V );
    for ( int i = 0; i < V; ++i )
        VertexShader( vertices[i], vertexPixels[i] );
    vector<Pixel> leftPixels;
    vector<Pixel> rightPixels;
    ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
    DrawRows( surface, leftPixels, rightPixels );
}

void ComputePolygonRows(
    const vector<Pixel>& vertexPixels,
    vector<Pixel>& leftPixels,
    vector<Pixel>& rightPixels )
{
    // Find max and min y-value of the polygon
    //    and compute the number of rows it occupies.
    int maxy = -numeric_limits<int>::max();
    int miny = +numeric_limits<int>::max();
    for (int i = 0; i < vertexPixels.size(); i++) {
        if (vertexPixels[i].y < miny)
            miny = vertexPixels[i].y;
        if (vertexPixels[i].y > maxy)
            maxy = vertexPixels[i].y;
    }
    int rows = maxy - miny + 1;

    // Resize leftPixels and rightPixels
    //    so that they have an element for each row.
    leftPixels  = vector<Pixel>(rows);
    rightPixels = vector<Pixel>(rows);

    // Initialize the x-coordinates in leftPixels
    //    to some really large value and the x-coordinates
    //    in rightPixels to some really small value.
    for ( int i = 0; i < rows; ++i ) {
        leftPixels[i].x  = +numeric_limits<int>::max();
        rightPixels[i].x = -numeric_limits<int>::max();
    }

    // Loop through all edges of the polygon and use
    //    linear interpolation to find the x-coordinate for
    //    each row it occupies. Update the corresponding
    //    values in rightPixels and leftPixels.
    for (int i = 0; i < vertexPixels.size(); i++) {
        int j = (i + 1) % vertexPixels.size();
        std::vector<Pixel> result(abs(vertexPixels[j].y - vertexPixels[i].y) + 1);
        Interpolate(vertexPixels[i], vertexPixels[j], result);

        for (int iman = 0; iman < result.size(); iman++) {
            if (result[iman].x < leftPixels[result[iman].y - miny].x) {
                leftPixels[result[iman].y - miny] = result[iman];
            }

            if (result[iman].x > rightPixels[result[iman].y - miny].x) {
                rightPixels[result[iman].y - miny] = result[iman];
            }
        }
    }
}


void Interpolate( Pixel p1, Pixel p2, vector<Pixel>& result )
{
    ivec2 a(p1.x, p1.y);
    ivec2 b(p2.x, p2.y);
    int N = result.size();

    vec2 stepPixel = vec2(b.x - a.x, b.y - a.y) /
                     float(glm::max(N - 1, 1));

    vec3 stepPos = vec3(p2.pos3d * p2.zinv - p1.pos3d * p1.zinv ) /
                   float(glm::max(N - 1, 1));

    float stepZinv = (p2.zinv - p1.zinv) /
                     float(glm::max(N - 1, 1));

    vec2 itPix( a );
    vec3 itPos( p1.pos3d * p1.zinv );
    float itZinv =  p1.zinv;
    // printf("(%.4f,%.4f,%.4f)\n", itPos.x, itPos.y, itPos.z);
    // printf("%.1f, %.1f, %.1f %.1f\n", itPix.x, itPix.y, stepPixel.x, stepPixel.y);

    for ( int i = 0; i < N; ++i ) {
        result[i].x = (int)itPix.x;
        result[i].y = (int)itPix.y;
        itPix += stepPixel;

        result[i].zinv = itZinv;
        itZinv += stepZinv;

        result[i].pos3d = itPos / result[i].zinv;
        itPos += stepPos;
    }
}
void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result )
{
    int N = result.size();
    vec2 step = vec2(b - a) / float(glm::max(N - 1, 1));
    vec2 current( a );
    for ( int i = 0; i < N; ++i ) {
        result[i] = current;
        current += step;
    }
}

void DrawLineSDL( SDL_Surface* surface, Pixel a, Pixel b, vec3 color )
{
    Pixel delta = glm::abs( a - b );
    int pixels = glm::max( delta.x, delta.y ) + 1;
    // You can then get the pixel positions of the line by calling the Interpolation function:
    vector<Pixel> line( pixels );
    Interpolate( a, b, line );
    for (int i = 0; i < line.size(); i++) {
        PutPixelSDL( surface, line[i].x, line[i].y, color );
    }
}

void Draw()
{
    SDL_FillRect( screen, 0, 0 );
    if (SDL_MUSTLOCK(screen))
        SDL_LockSurface(screen);

    // Compute frame time:
    // t1 = t2;
    // float dtsec = dt / 1000.0f;
    int t1 = SDL_GetTicks();

    for ( int i = 0; i < triangles.size(); ++i ) {
        vector<Vertex> vertices(3);
        vertices[0].position = triangles[i].v0;
        vertices[1].position = triangles[i].v1;
        vertices[2].position = triangles[i].v2;

        currentNormal = triangles[i].normal;
        currentReflectance = vec3(1, 1, 1);
        currentColor = triangles[i].color;

        DrawPolygon(screen, vertices);
    }
    int t2 = SDL_GetTicks();

    for (int y = 0; y < SCREEN_HEIGHT; y++) {
        for (int x = 0; x < SCREEN_WIDTH; x++) {
            // PutPixelSDL(screen, x, y, vec3(depthBuffer[y][x], 0, 0));
            // printf("%.4f\n", depthBuffer[y][x] * 100);
            PutPixelSDL(screen, x, y, image[x][y]);
            image[x][y] = vec3(0, 0, 0);
        }
    }


    if ( SDL_MUSTLOCK(screen) )
        SDL_UnlockSurface(screen);

    SDL_UpdateRect( screen, 0, 0, 0, 0 );

    int t3 = SDL_GetTicks();
    int computeTime = t2 - t1;
    int renderTime = t3 - t1;

    printf("%d ms compute, %d ms screen, %.4f ratio \r", computeTime, renderTime, computeTime / ((float)renderTime + computeTime));
}
