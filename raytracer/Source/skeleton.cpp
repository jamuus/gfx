#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <stdlib.h>     /* srand, rand */


using namespace std;
using glm::vec3;
using glm::mat3;
#define PI 3.14159265358979323846f
/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 600;
const int SCREEN_HEIGHT = 600;
const int NUM_RAYS = 10000;
const int NUM_BUCKETS = 200;
const int NUM_BOUNCES = 100;
vec3 *grid;

SDL_Surface* screen;
int t;

vec3 cameraPos(0.0f, 0.0f, -2.1f);
/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

struct Intersection {
    vec3 position;
    float distance;
    int triangleIndex;
};

void SetBucket(vec3 a, vec3 power, int depth);
vec3 GetBucket(vec3 a, int depth);
void InitBucket();
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

int main( int argc, char* argv[] )
{
    InitBucket();
    LoadTestModel(model);
    screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
    t = SDL_GetTicks(); // Set start value for timer.


    while ( NoQuitMessageSDL() ) {
        Update();
        Draw();
    }
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

void InitBucket()
{
    grid = (vec3*) calloc(sizeof(vec3), (NUM_BUCKETS * NUM_BUCKETS * NUM_BUCKETS * NUM_BOUNCES));
    for (int i = 0; i < NUM_BUCKETS * NUM_BUCKETS * NUM_BUCKETS * NUM_BOUNCES; i++) {
        grid[i] = vec3(0.0f, 0.0f, 0.0f);
    }
}

float maxdim = 1.01;
float filterSize = 5.0f;

void SetBucket(vec3 a, vec3 power, int depth)
{
    vec3 temp = (a + maxdim) / (2 * maxdim) * (float)NUM_BUCKETS;
    temp.x = (int)temp.x, temp.y = (int)temp.y, temp.z = (int)temp.z;
    int position = (int) ((NUM_BUCKETS * NUM_BUCKETS * NUM_BUCKETS * depth) + (NUM_BUCKETS * NUM_BUCKETS * temp.z) + (NUM_BUCKETS * temp.y) + temp.x);
    // grid[position] = power;
    grid[position] = grid[position] - (grid[position] / filterSize) + power / filterSize;
}

// void IncBucket(vec3 a, vec3 power)
// {
//     vec3 oldval = GetBucket(a);
//     SetBucket(a, oldval + power);
// }

vec3 GetBucket(vec3 a, int depth)
{
    vec3 temp = (a + maxdim) / (2 * maxdim) * (float)NUM_BUCKETS;
    temp.x = (int)temp.x, temp.y = (int)temp.y, temp.z = (int)temp.z;
    int position = (int) ((NUM_BUCKETS * NUM_BUCKETS * NUM_BUCKETS * depth) + (NUM_BUCKETS * NUM_BUCKETS * temp.z) + (NUM_BUCKETS * temp.y) + temp.x);
    vec3 spixy = grid[position];
    return spixy;
}

vector<vec3> GenerateRays()
{
    srand (time(NULL));
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

    printf("Render time: %.4f ms - %.4f fps   \r", dt, (1 / dtsec));
    flush(cout);
}

vec3 topLeft(1, 0, 0);  // red
vec3 topRight(0, 0, 1);  // blue
vec3 bottomRight(0, 1, 0); // green
vec3 bottomLeft(1, 1, 0); // yellow?

float f = SCREEN_HEIGHT / 2;

void shootRay(vec3 pos, vec3 dir, float totDist, int depth, vec3 raycolor)
{
    Intersection i;
    if (depth < NUM_BOUNCES) {
        if (ClosestIntersection(pos,
                                dir,
                                model,
                                i)) {
            float r = dist(pos, i.position) + totDist;

            float A = 4.0f * PI * r * r;
            vec3 B = raycolor / A;
            vec3 nhat = model[i.triangleIndex].normal;
            nhat = nhat / dist(vec3(0, 0, 0), nhat);

            vec3 rhat = pos - i.position;
            rhat = rhat / dist(vec3(0, 0, 0), rhat);

            vec3 D = B * max(dot(rhat, nhat), 0.0f);

            raycolor = model[i.triangleIndex].color * raycolor;
            raycolor = raycolor / dist(vec3(0, 0, 0), raycolor) * reflectedPower;

            SetBucket(i.position, D * raycolor, depth);

            vec3 reflecteddir = dir - 2.0f * nhat * dot(dir, nhat);
            shootRay(i.position, reflecteddir, r, depth + 1, raycolor);
        } else {
            // do noting
        }
    }
}

void Draw()
{
    SDL_FillRect( screen, 0, 0 );

    if ( SDL_MUSTLOCK(screen) )
        SDL_LockSurface(screen);

    vector<vec3> rays = GenerateRays();
    for (int i = 0; i < NUM_RAYS; i++) {


        shootRay(lightPos, rays[i], 0.0f, 0, lightColor);
    }

//**--------------------------------------------------**//
    for (int y = 0; y < SCREEN_HEIGHT; y++) {
        for (int x = 0; x < SCREEN_WIDTH; x++) {
            vec3 rayDir(x - SCREEN_WIDTH / 2, y - SCREEN_HEIGHT / 2, f);
            vec3 rotatedRay = rayDir * R;

            vec3 colour;

            Intersection cameraWorldInt;
            if (ClosestIntersection(cameraPos,
                                    rotatedRay,
                                    model,
                                    cameraWorldInt)) {

                Triangle t = model[cameraWorldInt.triangleIndex];

                vec3 color = vec3(0, 0, 0);
                for (int i = 0; i < NUM_BOUNCES; i++) {
                    color += GetBucket(cameraWorldInt.position, i);
                }
                colour = color;

                // colour = t.color * DirectLight(cameraWorldInt);
                // cout << colour.x << "\n" << colour.y << "\n" << colour.z << endl;
                // vec3 lightDir = lightPos - cameraWorldInt.position;
                // lightDir /= dist(vec3(0, 0, 0), lightDir);

                // Intersection shadowIntersection;
                // if (ClosestIntersection(cameraWorldInt.position,
                //                         lightDir,
                //                         model,
                //                         shadowIntersection)
                //         &&
                //         // dist(shadowIntersection.position, cameraWorldInt.position)
                //         shadowIntersection.distance
                //         <
                //         dist(cameraWorldInt.position, lightPos)
                //    ) {
                //     colour = vec3(0.0f, 0.0f, 0.0f) + indirectLight * t.color;
                // } else {
                // vec3 D = DirectLight(cameraWorldInt);
                // colour = (D + indirectLight) * t.color;
                // }
            } else {
                colour = vec3(0.0, 0.0, 0.0);
            }
            PutPixelSDL(screen, x, y, colour);
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