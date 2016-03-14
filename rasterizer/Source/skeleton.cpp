#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using namespace glm;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 800;
const int SCREEN_HEIGHT = 600;

#define PI 3.14159265358979323846f
SDL_Surface* screen;
int t;
vector<Triangle> triangles;
vec3 currentColor;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];


vec3 lightPos(0, -0.5, -0.7);
vec3 lightPower = 12.f * vec3( 1, 1, 1 );
vec3 indirectLightPowerPerArea = 0.3f * vec3( 1, 1, 1 );

float focalLength = SCREEN_HEIGHT;
vec3 cameraPosition(0, 0, -10);

vec3 currentNormal;
vec3 currentReflectance;


float translateSpeed = 1.0f;
float rotateSpeed = 0.1f;
mat3 R;

int dx;
int dy;

struct Pixel {
    int x;
    int y;
    float zinv;
    // vec3 illumination;
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


    while ( NoQuitMessageSDL() ) {
        Update();
        Draw();
        for (int y = 0; y < SCREEN_HEIGHT; ++y)
            for (int x = 0; x < SCREEN_WIDTH; ++x)
                depthBuffer[y][x] = 0;
    }

    SDL_SaveBMP( screen, "screenshot.bmp" );
    return 0;
}

float yaw = 0.0f;

void Update()
{
    // Compute frame time:
    int t2 = SDL_GetTicks();
    float dt = float(t2 - t);
    t = t2;
    float dtsec = dt / 1000.0f;


    SDL_GetRelativeMouseState( &dx, &dy );

    Uint8* keystate = SDL_GetKeyState( 0 );
    if ( keystate[SDLK_UP] ) {
        // Move camera forward
        cameraPosition += (vec3(0, 0, 1) * dtsec * translateSpeed) * R;
    }
    if ( keystate[SDLK_DOWN] ) {
        // Move camera backward
        cameraPosition -= (vec3(0, 0, 1) * dtsec * translateSpeed) * R;
    }
    if ( keystate[SDLK_LEFT] ) {
        // Move camera to the left
        cameraPosition -= (vec3(1, 0, 0) * dtsec * translateSpeed) * R;
        // rotate left
        // yaw -= dtsec * rotateSpeed;
    }
    if ( keystate[SDLK_RIGHT] ) {
        // Move camera to the right
        cameraPosition += (vec3(1, 0, 0) * dtsec * translateSpeed) * R;
    }
    yaw -= dx * dtsec * rotateSpeed;

    R[0] = vec3( cos(yaw), 0, sin(yaw));
    R[1] = vec3( 0,        1, 0);
    R[2] = vec3(-sin(yaw), 0, cos(yaw));

    // cout << "Render time: " << dt << " ms." << endl;
    printf("Render time: %.4f ms - %.4f fps   \r", dt, (1 / dtsec));
}

void PixelShader(const Pixel& p)
{
    int x = p.x;
    int y = p.y;
    if (x < 0 || x >= SCREEN_WIDTH || y < 0 || y >= SCREEN_HEIGHT) {
        printf("(%d, %d)\n", x, y);
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

        PutPixelSDL(screen, x, y, illumination);
        // PutPixelSDL(screen, x, y, currentColor);
    }
}

void VertexShader( const Vertex& v, Pixel& p )
{
    vec3 newV = R * (v.position - cameraPosition);
    p.x = (int)(focalLength * newV.x / newV.z + SCREEN_WIDTH / 2);
    p.y = (int)(focalLength * newV.y / newV.z + SCREEN_HEIGHT / 2);
    p.zinv = 1.0f / length(newV);
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


void Interpolate( Pixel a1, Pixel b1, vector<Pixel>& result )
{
    vec3 a(a1.x, a1.y, a1.zinv);
    vec3 b(b1.x, b1.y, b1.zinv);
    int N = result.size();

    vec3 step1 = vec3(b - a) / float(glm::max(N - 1, 1));
    vec3 step2 = vec3(b1.pos3d * b1.zinv - a1.pos3d * a1.zinv) / float(glm::max(N - 1, 1));

    vec3 current1( a );
    vec3 current2( a1.pos3d * a1.zinv );
    // printf("(%.4f,%.4f,%.4f)\n", current2.x, current2.y, current2.z);

    for ( int i = 0; i < N; ++i ) {
        result[i] = (Pixel) {current1.x, current1.y, current1.z};
        current1 += step1;

        result[i].pos3d = current2 / result[i].zinv;
        current2 += step2;
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
    if ( SDL_MUSTLOCK(screen) )
        SDL_LockSurface(screen);

    for ( int i = 0; i < triangles.size(); ++i ) {
        vector<Vertex> vertices(3);
        vertices[0].position = triangles[i].v0;
        vertices[1].position = triangles[i].v1;
        vertices[2].position = triangles[i].v2;

        currentNormal = triangles[i].normal;
        currentReflectance = vec3(1, 1, 1);
        currentColor = triangles[i].color;

        // vertices[0].normal = triangles[i].normal;
        // vertices[1].normal = triangles[i].normal;
        // vertices[2].normal = triangles[i].normal;

        // vertices[0].reflectance = vec3(1, 1, 1);
        // vertices[1].reflectance = vec3(1, 1, 1);
        // vertices[2].reflectance = vec3(1, 1, 1);
        DrawPolygon(screen, vertices);
    }


    if ( SDL_MUSTLOCK(screen) )
        SDL_UnlockSurface(screen);

    SDL_UpdateRect( screen, 0, 0, 0, 0 );
}
