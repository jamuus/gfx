#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <stdlib.h>     /* srand, rand */
#include <omp.h>


using namespace std;
using glm::vec3;
using glm::mat3;
#define PI 3.14159265358979323846f
/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 300;
const int SCREEN_HEIGHT = 300;
const int NUM_RAYS = 10000;
const int NUM_SAMPLES = 1;
#define MAXDEPTH 1

vec3 image[SCREEN_WIDTH * SCREEN_HEIGHT * NUM_SAMPLES];

SDL_Surface* screen;
int t;

vec3 cameraPos(0.0f, 0.0f, -2.0f);
/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

struct Intersection {
    vec3 position;
    float distance;
    int triangleIndex;
};

// void SetBucket(vec3 a, vec3 power, int depth);
// vec3 GetBucket(vec3 a, int depth);

vector<vec3> GenerateRays();
void Update();
void Draw();
void Interpolate( float a, float b, vector<float>& result );
void Interpolate( vec3 a, vec3 b, vector<vec3>& result );
bool ClosestIntersection(
    vec3 start,
    vec3 dir,
    const vector<Triangle>& triangles,
    Intersection& closestIntersection );
vec3 DirectLight( const Intersection& i );
vec3 DirectLight(const Intersection& i, vec3 lightSource, vec3 power, float);
float dist(vec3 a, vec3 b);

vector<Triangle> model;

void InitImage()
{
    for (int i = 0; i < SCREEN_WIDTH * SCREEN_HEIGHT; i++) {
        for (int j = 0; j < NUM_SAMPLES; j++) {
            image[i * NUM_SAMPLES + j] = vec3(0, 0, 0);
        }
    }
}

int main( int argc, char* argv[] )
{
    srand (time(NULL));
    InitImage();
    LoadTestModel(model);
    screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
    t = SDL_GetTicks(); // Set start value for timer.


    // while ( NoQuitMessageSDL() ) {
    Update();
    Draw();
    // }
    cout << endl;

    SDL_SaveBMP( screen, "screenshot.bmp" );
    return 0;
}

float translateSpeed = 1.0f;
float rotateSpeed = 1.0f;
float lightTranslateDist = 0.1f;
mat3 R;
float yaw = 0.0f;
vec3 lightPos(0, -0.5, -0.7);
vec3 indirectLight = 0.5f * vec3( 1, 1, 1 );
vec3 lightColor = 5.f * vec3(1, 1, 1);
float reflectedPower = 3.0f;

// void InitBucket()
// {
//     grid = (vec3*) calloc(sizeof(vec3), (NUM_BUCKETS * NUM_BUCKETS * NUM_BUCKETS * NUM_BOUNCES));
//     for (int i = 0; i < NUM_BUCKETS * NUM_BUCKETS * NUM_BUCKETS * NUM_BOUNCES; i++) {
//         grid[i] = vec3(0.0f, 0.0f, 0.0f);
//     }
// }

// float maxdim = 1.01;
// float filterSize = 5.0f;

// void SetBucket(vec3 a, vec3 power, int depth)
// {
//     vec3 temp = (a + maxdim) / (2 * maxdim) * (float)NUM_BUCKETS;
//     temp.x = (int)temp.x, temp.y = (int)temp.y, temp.z = (int)temp.z;
//     int position = (int) ((NUM_BUCKETS * NUM_BUCKETS * NUM_BUCKETS * depth) + (NUM_BUCKETS * NUM_BUCKETS * temp.z) + (NUM_BUCKETS * temp.y) + temp.x);
//     // grid[position] = power;
//     grid[position] = grid[position] - (grid[position] / filterSize) + power / filterSize;
// }

// // void IncBucket(vec3 a, vec3 power)
// // {
// //     vec3 oldval = GetBucket(a);
// //     SetBucket(a, oldval + power);
// // }

// vec3 GetBucket(vec3 a, int depth)
// {
//     vec3 temp = (a + maxdim) / (2 * maxdim) * (float)NUM_BUCKETS;
//     temp.x = (int)temp.x, temp.y = (int)temp.y, temp.z = (int)temp.z;
//     int position = (int) ((NUM_BUCKETS * NUM_BUCKETS * NUM_BUCKETS * depth) + (NUM_BUCKETS * NUM_BUCKETS * temp.z) + (NUM_BUCKETS * temp.y) + temp.x);
//     vec3 spixy = grid[position];
//     return spixy;
// }

vector<vec3> GenerateRays()
{
    // srand (time(NULL));
    vector<vec3> rays;
    for (int i = 0; i < NUM_RAYS; i++) {
        float z = rand() / (RAND_MAX / 2.0f) - 1;
        float rxy = sqrt(1 - z * z);
        float phi = rand() / (RAND_MAX / (2 * PI));
        float x = rxy * cos(phi);
        float y = rxy * sin(phi);
        vec3 ray(x, y, z);
        // ray /= 10.0f;
        rays.push_back(ray);
    }
    return rays;
}

vector<vec3> GenerateRays(vec3 norm)
{
    // srand (time(NULL));

    vector<vec3> rays;
    for (int i = 0; i < NUM_RAYS; i++) {
        float theta = ((rand() / RAND_MAX) * 2 * PI) - PI;
        float azi   = (rand() / RAND_MAX) * PI - PI;
        float x = sin(theta) * cos(azi);
        float y = sin(theta) * sin(azi);
        float z = cos(theta);

        vec3 o(1, 0, 0);
        vec3 v = cross(norm, o);

        mat3 v_x(0, v.z, v.y,
                 v.z, 0, -v.x,
                 -v.y, v.x, 0);

        float s = length(v);
        float c = dot(v, o);

        mat3 R = mat3(1, 0, 0,
                      0, 1, 0,
                      0, 0, 1) +
                 v_x + v_x * v_x *
                 ((1.0f - c) / (s * s));

        vec3 ray(x, y, z);
        ray /= length(ray);

        rays.push_back(ray * R);
    }
    return rays;
}

void Update()
{
    // Compute frame time:
    int t2 = SDL_GetTicks();
    float dt = float(t2 - t);
    t = t2;

    float dtsec = dt / 1000.0f;

    Uint8* keystate = SDL_GetKeyState( 0 );
    if ( keystate[SDLK_UP] ) {
        // Move camera forward
        cameraPos += (vec3(0, 0, 1) * dtsec * translateSpeed) * R;
    }
    if ( keystate[SDLK_DOWN] ) {
        // Move camera backward
        cameraPos -= (vec3(0, 0, 1) * dtsec * translateSpeed) * R;
    }
    if ( keystate[SDLK_LEFT] ) {
        // Move camera to the left
        // cameraPos.x -= dtsec * translateSpeed;
        // rotate left
        yaw -= dtsec * rotateSpeed;
    }
    if ( keystate[SDLK_RIGHT] ) {
        // Move camera to the right
        // cameraPos.x += dtsec * translateSpeed;
        // rotate right
        yaw += dtsec * rotateSpeed;
    }

    if ( keystate[SDLK_w] ) {
        lightPos.z += lightTranslateDist;
    }
    if ( keystate[SDLK_a] ) {
        lightPos.x -= lightTranslateDist;
    }
    if ( keystate[SDLK_s] ) {
        lightPos.z -= lightTranslateDist;
    }
    if ( keystate[SDLK_d] ) {
        lightPos.x += lightTranslateDist;
    }
    if ( keystate[SDLK_q] ) {
        lightPos.y += lightTranslateDist;
    }
    if ( keystate[SDLK_e] ) {
        lightPos.y -= lightTranslateDist;
    }

    R[0] = vec3( cos(yaw), 0, sin(yaw));
    R[1] = vec3( 0,        1, 0);
    R[2] = vec3(-sin(yaw), 0, cos(yaw));

    // printf("Render time: %.4f ms - %.4f fps   \r", dt, (1 / dtsec));
    flush(cout);
}

vec3 topLeft(1, 0, 0);  // red
vec3 topRight(0, 0, 1);  // blue
vec3 bottomRight(0, 1, 0); // green
vec3 bottomLeft(1, 1, 0); // yellow?

float f = SCREEN_HEIGHT / 2;


vec3 shootRay(vec3 pos, vec3 dir, float totDist, int depth)
{
    Intersection i;
    if (depth > 0) {
        // calc where the ray intersects the world
        if (ClosestIntersection(pos,
                                dir,
                                model,
                                i)) {
            Triangle t = model[i.triangleIndex];

            if (t.lightSource) {
                return t.intensity;
            } else {
                vec3 nhat = t.normal;
                nhat = nhat / dist(vec3(0, 0, 0), nhat);

                float r = dist(pos, i.position) + totDist;
                float A = 4.0f * PI * r * r;

                vec3 reflecteddir = dir - 2.0f * nhat * dot(dir, nhat);

                vec3 powOfNextPoint = shootRay(i.position, reflecteddir, r, depth - 1);

                vec3 B = powOfNextPoint / A;

                vec3 rhat = pos - i.position;
                rhat = rhat / dist(vec3(0, 0, 0), rhat);

                vec3 D = B *
                         // fraction of light that is reflected
                         max(dot(rhat, nhat), 0.0f);

                return D;
            }
        } else {
            // no intersection return black
            return vec3(0, 0, 0);
        }
    } else {
        // if run out of depth return black
        return vec3(0, 0, 0);
    }
}


void Draw()
{
    int numthreads = 4;
    omp_set_num_threads(numthreads);
    printf("numthreads: %d\n", numthreads);


    SDL_FillRect( screen, 0, 0 );

    if ( SDL_MUSTLOCK(screen) )
        SDL_LockSurface(screen);

    #pragma omp parallel for
    for (int y = 0; y < SCREEN_HEIGHT; y++) {
        for (int x = 0; x < SCREEN_WIDTH; x++) {
            vec3 rayDir(x - SCREEN_WIDTH / 2, y - SCREEN_HEIGHT / 2, f);
            vec3 rotatedRay = rayDir * R;

            vec3 colour(0, 0, 0);
            Intersection i;
            if (ClosestIntersection(cameraPos,
                                    rotatedRay,
                                    model,
                                    i)) {
                Triangle t = model[i.triangleIndex];
                if (t.lightSource) {
                    colour = t.intensity;
                } else {
                    vector<vec3> rays = GenerateRays();
                    for (int j = 0; j < rays.size(); j++) {
                        colour += shootRay(i.position, rays[j], 0.0f, MAXDEPTH);
                    }
                    colour /= rays.size();
                }
                // image[x * SCREEN_WIDTH + y] = colour;
            }
            PutPixelSDL(screen, x, y, colour);
            if (omp_get_thread_num() == 0)
                printf("%.4f\%\r", numthreads * 100.0 * (y * (float)SCREEN_HEIGHT + x) / (SCREEN_HEIGHT * SCREEN_WIDTH));
        }
    }

    if (SDL_MUSTLOCK(screen))
        SDL_UnlockSurface(screen);

    SDL_UpdateRect(screen, 0, 0, 0, 0);
}


float m = std::numeric_limits<float>::max();

bool ClosestIntersection(
    vec3 start,
    vec3 dir,
    const vector<Triangle>& triangles,
    Intersection& closestIntersection)
{
    start += dir * 0.00001f;
    vec3 closestX(m, 0.0f, 0.0f);
    int index = -1;
    bool found = false;

    for (int i = 0; i < triangles.size(); i++) {
        Triangle tri = triangles[i];
        vec3 e1 = tri.v1 - tri.v0;
        vec3 e2 = tri.v2 - tri.v0;
        vec3 b = start - tri.v0;
        mat3 A(-dir, e1, e2);
        vec3 x = glm::inverse(A) * b;
        if (
            x.y >= 0.0f &&
            x.z >= 0.0f &&
            x.y + x.z <= 1.0f &&
            x.x < closestX.x &&
            x.x >= 0.0f
        ) {
            closestX = x;
            index = i;
            found = true;
        }
    }

    closestIntersection.position = start + dir * closestX.x;
    closestIntersection.distance = closestX.x;
    closestIntersection.triangleIndex = index;

    return found;
}

float dist(vec3 a, vec3 b)
{
    return sqrt(pow(a.x - b.x, 2) +
                pow(a.y - b.y, 2) +
                pow(a.z - b.z, 2));
}


vec3 DirectLight(const Intersection& i)
{
    return DirectLight(i, lightPos, lightColor, 0.0f);
}

vec3 DirectLight(const Intersection& i, vec3 lightSource, vec3 power, float dr)
{
    float r = dist(lightSource, i.position) + dr;

    float A = 4.0f * PI * r * r;
    vec3 B = power / A;
    vec3 nhat = model[i.triangleIndex].normal;
    nhat = nhat / dist(vec3(0, 0, 0), nhat);
    vec3 rhat = lightSource - i.position;
    rhat = rhat / dist(vec3(0, 0, 0), rhat);
    vec3 D = B * max(dot(rhat, nhat), 0.0f);

    return D;
}


void Interpolate(float a, float b, vector<float>& result)
{
    if (result.size() == 1) {
        result[0] = (b - a) / 2;
        return;
    }
    float diff = (b - a) / (result.size() - 1);
    for (int i = 0 ; i < result.size(); i++) {
        result[i] = a + diff * i;
    }
}


void Interpolate(vec3 a, vec3 b, vector<vec3>& result)
{
    vec3 diff(
        (b.x - a.x) / (result.size() - 1),
        (b.y - a.y) / (result.size() - 1),
        (b.z - a.z) / (result.size() - 1)
    );

    for (int i = 0; i < result.size(); i++) {
        // result[i].x = a.x + diff.x * i;
        // result[i].y = a.y + diff.y * i;
        // result[i].z = a.z + diff.z * i;

        result[i] = a + diff * (float)i;
    }
}