#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
// using glm::vec3;
// using glm::mat3;
using namespace glm;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;
vector<Triangle> triangles;
vec3 currentColor;
float depthBuffer[SCREEN_WIDTH][SCREEN_HEIGHT];

struct Pixel {
    int x;
    int y;
    float zinv;

    inline Pixel operator-(Pixel a)
    {
        return (Pixel) {x - a.x, y - a.y, zinv - a.zinv};
    }
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

// void DrawPolygonRows(
//     const vector<Pixel>& leftPixels,
//     const vector<Pixel>& rightPixels );


void VertexShader( const vec3& v, Pixel& p );

// void ComputePolygonRows(
//     const vector<ivec2>& vertexPixels,
//     vector<ivec2>& leftPixels,
//     vector<ivec2>& rightPixels );

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

void Update()
{
    // Compute frame time:
    int t2 = SDL_GetTicks();
    float dt = float(t2 - t);
    t = t2;
    cout << "Render time: " << dt << " ms." << endl;
}

float focalLength = SCREEN_HEIGHT;
vec3 cameraPosition(0, 0, -3.01);
void VertexShader( const vec3& v, Pixel& p )
{
    vec3 newV = v - cameraPosition;
    int x = (int)(focalLength * newV.x / newV.z + SCREEN_WIDTH / 2);
    int y = (int)(focalLength * newV.y / newV.z + SCREEN_HEIGHT / 2);
    p.x = x;
    p.y = y;
    p.zinv = 1.0f / length(newV);
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
            if (p.zinv > depthBuffer[p.x][p.y]) {
                PutPixelSDL(surface, p.x, p.y, currentColor);
                depthBuffer[p.x][p.y] = p.zinv;
            }
        }
        // for (int x = lp.x; x < rightPixels[row].x; x++) {
        //     if()
        //     PutPixelSDL(surface, x, y, currentColor);
        // }
    }
}

void DrawPolygon( SDL_Surface* surface,
                  const vector<vec3>& vertices )
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
// 1. Find max and min y-value of the polygon
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
// 2. Resize leftPixels and rightPixels
//    so that they have an element for each row.
    leftPixels  = vector<Pixel>(rows);
    rightPixels = vector<Pixel>(rows);
// 3. Initialize the x-coordinates in leftPixels
//    to some really large value and the x-coordinates
//    in rightPixels to some really small value.
    for ( int i = 0; i < rows; ++i ) {
        leftPixels[i].x  = +numeric_limits<int>::max();
        rightPixels[i].x = -numeric_limits<int>::max();
    }
// 4. Loop through all edges of the polygon and use
//    linear interpolation to find the x-coordinate for
//    each row it occupies. Update the corresponding
//    values in rightPixels and leftPixels.
    for (int i = 0; i < vertexPixels.size(); i++) {
        int j = (i + 1) % vertexPixels.size();
        std::vector<Pixel> result(abs(vertexPixels[j].y - vertexPixels[i].y) + 1);
        Interpolate(vertexPixels[i], vertexPixels[j], result);

        for (int iman = 0; iman < result.size(); iman++) {
            if (result[iman].x < leftPixels[result[iman].y - miny].x) {
                leftPixels[result[iman].y - miny].x = result[iman].x;
                leftPixels[result[iman].y - miny].y = result[iman].y;
                leftPixels[result[iman].y - miny].zinv = result[iman].zinv;
            }

            if (result[iman].x > rightPixels[result[iman].y - miny].x) {
                rightPixels[result[iman].y - miny].x = result[iman].x;
                rightPixels[result[iman].y - miny].y = result[iman].y;
                rightPixels[result[iman].y - miny].zinv = result[iman].zinv;
            }
        }
    }
}


void Interpolate( Pixel a1, Pixel b1, vector<Pixel>& result )
{
    vec3 a(a1.x, a1.y, a1.zinv);
    vec3 b(b1.x, b1.y, b1.zinv);
    int N = result.size();
    vec3 step = vec3(b - a) / float(glm::max(N - 1, 1));
    vec3 current( a );
    for ( int i = 0; i < N; ++i ) {
        result[i] = (Pixel) {current.x, current.y, current.z};
        current += step;
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

void DrawPolygonEdges( const vector<vec3>& vertices )
{
    int V = vertices.size();
    // Transform each vertex from 3D world position to 2D image position:
    vector<Pixel> projectedVertices( V );
    for ( int i = 0; i < V; ++i ) {
        VertexShader( vertices[i], projectedVertices[i] );
    }
    // Loop over all vertices and draw the edge from it to the next vertex:
    for ( int i = 0; i < V; ++i ) {
        int j = (i + 1) % V; // The next vertex
        vec3 color( 1, 1, 1 );
        DrawLineSDL( screen, projectedVertices[i], projectedVertices[j], color );
    }
}

void Draw()
{
    SDL_FillRect( screen, 0, 0 );
    if ( SDL_MUSTLOCK(screen) )
        SDL_LockSurface(screen);

    // for ( int y = 0; y < SCREEN_HEIGHT; ++y ) {
    //     for ( int x = 0; x < SCREEN_WIDTH; ++x ) {
    //         vec3 color( 1.0, 0.0, 0.0 );
    //         PutPixelSDL( screen, x, y, color );
    //     }
    // }
    for ( int i = 0; i < triangles.size(); ++i ) {
        vector<vec3> vertices(3);
        vertices[0] = triangles[i].v0;
        vertices[1] = triangles[i].v1;
        vertices[2] = triangles[i].v2;
        currentColor = triangles[i].color;
        // for (int v = 0; v < 3; ++v) {
        //     ivec2 projPos;
        //     VertexShader( vertices[v], projPos );
        //     vec3 color(1, 1, 1);
        //     PutPixelSDL( screen, projPos.x, projPos.y, color );
        // }
        DrawPolygon(screen, vertices);
    }


    if ( SDL_MUSTLOCK(screen) )
        SDL_UnlockSurface(screen);

    SDL_UpdateRect( screen, 0, 0, 0, 0 );
}
