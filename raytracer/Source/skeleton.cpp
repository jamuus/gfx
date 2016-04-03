#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <stdlib.h>     /* srand, rand */
#include <omp.h>


using namespace std;
// using glm::vec3;
// using glm::mat3;
using namespace glm;
#define PI 3.14159265358979323846f
/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 400;
const int SCREEN_HEIGHT = 400;
const int NUM_RAYS = 100000;
const int NUM_SAMPLES = 1;
#define MAXDEPTH 1

vec3 image[SCREEN_WIDTH * SCREEN_HEIGHT * NUM_SAMPLES];

SDL_Surface* screen;
int t;

vec3 cameraPos(0.0f, 0.0f, -2.01f);
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
    srand(time(NULL));
    // InitImage();
    LoadTestModel(model);
    screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT);
    t = SDL_GetTicks(); // Set start value for timer.

    Update();
    Draw();

    cout << endl;

    SDL_SaveBMP(screen, "screenshot.bmp");
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
// float reflectedPower = 3.0f;

void printMatrix(mat3 R)
{
    printf("[%.4f,%.4f,%.4f] \n", R[0].x, R[0].y, R[0].z);
    printf("[%.4f,%.4f,%.4f] \n", R[1].x, R[1].y, R[1].z);
    printf("[%.4f,%.4f,%.4f] \n", R[2].x, R[2].y, R[2].z);
}

void printVector(vec3 V)
{
    printf("(%.4f, %.4f, %.4f)\n", V.x, V.y, V.z);
}

mat3 vectorRotation(vec3 a, vec3 b)
{
    vec3 v = cross(a, b);

    float s = length(v);

    mat3 v_x(0, -v.z, v.y,
             v.z, 0, -v.x,
             -v.y, v.x, 0);

    float c = dot(b, a);


    mat3 R = mat3(1, 0, 0,
                  0, 1, 0,
                  0, 0, 1);

    if (s != 0.0f) {
        R += v_x +
             v_x * v_x *
             (1.0f - c) / (s * s);
    } else {
        // if a == -b
        R = -mat3(1, 0, 0,
                  0, 1, 0,
                  0, 0, 1);
    }
    return R;
}

vec3 generateRandomSphereVector()
{
    float z = rand() / (RAND_MAX / 2.0f) - 1; // -1 and 1
    // float z = rand() / (RAND_MAX / 1.0f); // 0 and 1

    float rxy = sqrt(1 - z * z);
    float phi = rand() / (RAND_MAX / (2 * PI));
    float x = rxy * cos(phi);
    float y = rxy * sin(phi);

    vec3 ray(x, y, z);
    ray = normalize(ray);
    return ray;
}

vector<vec3> GenerateRays(vec3 norm)
{
    vec3 hemisphereTop(0.0f, 0.0f, 1.0f);

    mat3 R = vectorRotation(normalize(norm), normalize(hemisphereTop));

    vector<vec3> rays;
    for (int i = 0; i < NUM_RAYS; i++) {
        vec3 ray = generateRandomSphereVector();
        if (ray.z < 0)
            ray.z = -ray.z;

        vec3 directedRay = R * ray;
        rays.push_back(directedRay);
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

    R[0] = vec3( cos(yaw), 0, sin(yaw));
    R[1] = vec3( 0,        1, 0);
    R[2] = vec3(-sin(yaw), 0, cos(yaw));

}

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
                         glm::max(dot(rhat, nhat), 0.0f);

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
    int numthreads = 3;
    omp_set_num_threads(numthreads);
    printf("numthreads: %d\n", numthreads);


    SDL_FillRect( screen, 0, 0 );

    if ( SDL_MUSTLOCK(screen) )
        SDL_LockSurface(screen);

    // vec3 dir = normalize(vec3(0, 1, 0));
    // vector<vec3> rays = GenerateRays(dir);

    int done[numthreads];
    for (int i = 0; i < numthreads; i++) {
        done[i] = 0;
    }

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
                    // vector<vec3> rays = GenerateRays(dir);
                    vector<vec3> rays = GenerateRays(t.normal);
                    for (int j = 0; j < rays.size(); j++) {
                        colour += shootRay(i.position, rays[j], 0.0f, MAXDEPTH);
                    }
                    colour /= rays.size();
                }
                // image[x * SCREEN_WIDTH + y] = colour;
            }
            PutPixelSDL(screen, x, y, colour);
            done[omp_get_thread_num()]++;

            int total = 0;
            for (int i = 0; i < numthreads; i++) {
                total += done[i];
            }
            // if (omp_get_thread_num() == 0) {
            printf("%.2f%%   %d/%d   \r", total / ((float)SCREEN_HEIGHT * SCREEN_WIDTH) * 100.f, total, SCREEN_HEIGHT * SCREEN_WIDTH);
            fflush(stdout);
            // }
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

        float Adet = determinant(A);

        float A0 = determinant(mat2(A[1].y, A[1].z,
                                    A[2].y, A[2].z));

        float A3 = determinant(mat2(A[1].z, A[1].x,
                                    A[2].z, A[2].x));


        float A6 = determinant(mat2(A[1].x, A[1].y,
                                    A[2].x, A[2].y));

        vec3 Abot = vec3(A0, A3, A6) / Adet;
        float t = dot(Abot, b);

        if (t < closestX.x && t >= 0.0f) {
            float A1 = determinant(mat2(A[0].z, A[0].y,
                                        A[2].z, A[2].y));

            float A2 = determinant(mat2(A[0].y, A[0].z,
                                        A[1].y, A[1].z));

            float A4 = determinant(mat2(A[0].x, A[0].z,
                                        A[2].x, A[2].z));

            float A5 = determinant(mat2(A[0].z, A[0].x,
                                        A[1].z, A[1].x));

            float A7 = determinant(mat2(A[0].y, A[0].x,
                                        A[2].y, A[2].x));

            float A8 = determinant(mat2(A[0].x, A[0].y,
                                        A[1].x, A[1].y));

            vec3 Amid = vec3(A1, A4, A7) / Adet;
            float v = dot(Amid, b);

            vec3 Atop = vec3(A2, A5, A8) / Adet;
            float u = dot(Atop, b);
            vec3 x = vec3(t, u, v);

            if (
                u >= 0.0f &&
                v >= 0.0f &&
                u + v <= 1.0f
            ) {
                closestX = x;
                index = i;
                found = true;
            }
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
    vec3 D = B * glm::max(dot(rhat, nhat), 0.0f);

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
        result[i] = a + diff * (float)i;
    }
}