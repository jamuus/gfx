#ifndef TEST_MODEL_CORNEL_BOX_H
#define TEST_MODEL_CORNEL_BOX_H

// Defines a simple test model: The Cornel Box

#include <glm/glm.hpp>
#include <vector>

using glm::vec3;
// Used to describe a triangular surface:
class Triangle
{
public:
    glm::vec3 v0;
    glm::vec3 v1;
    glm::vec3 v2;
    glm::vec3 normal;
    glm::vec3 color;
    glm::vec3 intensity;
    bool lightSource;
    float diffuseK;
    float specularK;
    float alpha;
    float absorbtionFactor = 0.8;


    Triangle( glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, glm::vec3 color, float diffuseK, float specularK, float alpha)
        : v0(v0), v1(v1), v2(v2), color(color)
    {
        ComputeNormal();
        lightSource = false;
        this->diffuseK = diffuseK;
        this->specularK = specularK;
        this->alpha = alpha;
    }

    Triangle( glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, glm::vec3 color , glm::vec3 intensity)
        : v0(v0), v1(v1), v2(v2), color(color)
    {
        ComputeNormal();
        lightSource = true;
        this->intensity = intensity;
    }

    void ComputeNormal()
    {
        glm::vec3 e1 = v1 - v0;
        glm::vec3 e2 = v2 - v0;
        normal = glm::normalize( glm::cross( e2, e1 ) );
    }
};

// vec3 red(    0.75f, 0.15f, 0.15f );
vec3 red(    244 / 255.f, 67 / 255.f, 54 / 255.f );
vec3 yellow( 0.75f, 0.75f, 0.15f );
vec3 green(  0.15f, 0.75f, 0.15f );
vec3 cyan(   0.15f, 0.75f, 0.75f );
vec3 blue(   0.15f, 0.15f, 0.75f );
vec3 purple( 0.75f, 0.15f, 0.75f );
vec3 white(  0.75f, 0.75f, 0.75f );

float diffuseK = 0.4f;
float specularK = 0.6f;
float alpha = 10.f;

// Loads the Cornell Box. It is scaled to fill the volume:
// -1 <= x <= +1
// -1 <= y <= +1
// -1 <= z <= +1
void LoadTestModel( std::vector<Triangle>& triangles )
{
    // triangles.clear();
    // triangles.reserve( 5 * 2 * 3 );

    // ---------------------------------------------------------------------------
    // Room


    float L = 555;          // Length of Cornell Box side.

    vec3 A(L, 0, 0);
    vec3 B(0, 0, 0);
    vec3 C(L, 0, L);
    vec3 D(0, 0, L);

    vec3 E(L, L, 0);
    vec3 F(0, L, 0);
    vec3 G(L, L, L);
    vec3 H(0, L, L);

    // Floor:
    triangles.push_back( Triangle( C, B, A, white , diffuseK, specularK, alpha) );
    triangles.push_back( Triangle( C, D, B, white , diffuseK, specularK, alpha) );

    // Left wall
    triangles.push_back( Triangle( A, E, C, white , diffuseK, specularK, alpha) );
    triangles.push_back( Triangle( C, E, G, white , diffuseK, specularK, alpha) );

    // Right wall
    triangles.push_back( Triangle( F, B, D, white , diffuseK, specularK, alpha) );
    triangles.push_back( Triangle( H, F, D, white , diffuseK, specularK, alpha) );

    // Ceiling
    triangles.push_back( Triangle( E, F, G, white , diffuseK, specularK, alpha) );
    triangles.push_back( Triangle( F, H, G, white , diffuseK, specularK, alpha) );

    // Back wall
    triangles.push_back( Triangle( D, C, G, purple , diffuseK, specularK, alpha) );
    triangles.push_back( Triangle( G, H, D, purple , diffuseK, specularK, alpha) );

    // Light
    vec3 offset(L / 2 - 50, 350, L / 2 - 50);
    float lightIntensity = 500.0f;
    triangles.push_back( Triangle( E / 4.f + offset, F / 4.f + offset, G / 4.f + offset, white,
                                   lightIntensity * vec3(1, 1, 1)) );
    triangles.push_back( Triangle( F / 4.f + offset, H / 4.f + offset, G / 4.f + offset, white,
                                   lightIntensity * vec3(1, 1, 1)) );

    // ---------------------------------------------------------------------------
    // Short block

    A = vec3(290, 0, 114);
    B = vec3(130, 0, 65);
    C = vec3(240, 0, 272);
    D = vec3( 82, 0, 225);

    E = vec3(290, 165, 114);
    F = vec3(130, 165, 65);
    G = vec3(240, 165, 272);
    H = vec3( 82, 165, 225);

    float redBoxDiffuse = 0.4f;

    // Front
    triangles.push_back( Triangle(E, B, A, red, redBoxDiffuse, specularK, alpha) );
    triangles.push_back( Triangle(E, F, B, red, redBoxDiffuse, specularK, alpha) );

    // Front
    triangles.push_back( Triangle(F, D, B, red, redBoxDiffuse, specularK, alpha) );
    triangles.push_back( Triangle(F, H, D, red, redBoxDiffuse, specularK, alpha) );

    // BACK
    triangles.push_back( Triangle(H, C, D, red, redBoxDiffuse, specularK, alpha) );
    triangles.push_back( Triangle(H, G, C, red, redBoxDiffuse, specularK, alpha) );

    // LEFT
    triangles.push_back( Triangle(G, E, C, red, redBoxDiffuse, specularK, alpha) );
    triangles.push_back( Triangle(E, A, C, red, redBoxDiffuse, specularK, alpha) );

    // TOP
    triangles.push_back( Triangle(G, F, E, red, redBoxDiffuse, specularK, alpha) );
    triangles.push_back( Triangle(G, H, F, red, redBoxDiffuse, specularK, alpha) );

    // ---------------------------------------------------------------------------
    // Tall block

    A = vec3(423, 0, 247);
    B = vec3(265, 0, 296);
    C = vec3(472, 0, 406);
    D = vec3(314, 0, 456);

    E = vec3(423, 330, 247);
    F = vec3(265, 330, 296);
    G = vec3(472, 330, 406);
    H = vec3(314, 330, 456);

    // Front
    triangles.push_back( Triangle(E, B, A, blue, diffuseK, specularK, alpha) );
    triangles.push_back( Triangle(E, F, B, blue, diffuseK, specularK, alpha) );

    // Front
    triangles.push_back( Triangle(F, D, B, blue, diffuseK, specularK, alpha) );
    triangles.push_back( Triangle(F, H, D, blue, diffuseK, specularK, alpha) );

    // BACK
    triangles.push_back( Triangle(H, C, D, blue, diffuseK, specularK, alpha) );
    triangles.push_back( Triangle(H, G, C, blue, diffuseK, specularK, alpha) );

    // LEFT
    triangles.push_back( Triangle(G, E, C, blue, diffuseK, specularK, alpha) );
    triangles.push_back( Triangle(E, A, C, blue, diffuseK, specularK, alpha) );

    // TOP
    triangles.push_back( Triangle(G, F, E, blue, diffuseK, specularK, alpha) );
    triangles.push_back( Triangle(G, H, F, blue, diffuseK, specularK, alpha) );


    // ----------------------------------------------
    // Scale to the volume [-1,1]^3

    for ( size_t i = 0; i < triangles.size(); ++i ) {
        triangles[i].v0 *= 2 / L;
        triangles[i].v1 *= 2 / L;
        triangles[i].v2 *= 2 / L;

        triangles[i].v0 -= vec3(1, 1, 1);
        triangles[i].v1 -= vec3(1, 1, 1);
        triangles[i].v2 -= vec3(1, 1, 1);

        triangles[i].v0.x *= -1;
        triangles[i].v1.x *= -1;
        triangles[i].v2.x *= -1;

        triangles[i].v0.y *= -1;
        triangles[i].v1.y *= -1;
        triangles[i].v2.y *= -1;

        triangles[i].ComputeNormal();
    }
}
void LoadEnclosedTestModel( std::vector<Triangle>& triangles )
{
    // using glm::vec3;

    float L = 555;          // Length of Cornell Box side.

    // vec3 white( 1.f, 1.f, 1.f );

    vec3 A(L, 0, 0);
    vec3 B(0, 0, 0);
    vec3 C(L, 0, L);
    vec3 D(0, 0, L);

    vec3 E(L, L, 0);
    vec3 F(0, L, 0);
    vec3 G(L, L, L);
    vec3 H(0, L, L);


    // front:
    triangles.push_back( Triangle( E, F, A, white , diffuseK, specularK, alpha) );
    triangles.push_back( Triangle( B, F, A, white , diffuseK, specularK, alpha) );

    LoadTestModel(triangles);
}

#endif