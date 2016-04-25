#ifndef TEST_MODEL_CORNEL_BOX_H
#define TEST_MODEL_CORNEL_BOX_H

// Defines a simple test model: The Cornel Box

#include <glm/glm.hpp>
#include <vector>

using namespace glm;

vec2 bl(0, 0);
vec2 tl(0, 1);
vec2 tr(1, 1);
vec2 br(1, 0);



// Used to describe a triangular surface:
class Triangle
{
public:
    glm::vec3 v0;
    glm::vec2 uv0;

    glm::vec3 v1;
    glm::vec2 uv1;

    glm::vec3 v2;
    glm::vec2 uv2;

    glm::vec3 normal;
    glm::vec3 color;
    glm::vec3 **texture;
    ivec2 texSize;

    Triangle( glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, glm::vec3 color, vec2 uv0, vec2 uv1, vec2 uv2, unsigned texseed)
        : v0(v0), v1(v1), v2(v2), uv0(uv0), uv1(uv1), uv2(uv2), color(color)
    {
        srand(texseed);
        vec3 w(rand() / (float)RAND_MAX, rand() / (float)RAND_MAX, rand() / (float)RAND_MAX);
        vec3 b(rand() / (float)RAND_MAX, rand() / (float)RAND_MAX, rand() / (float)RAND_MAX);
        texSize.x = 5;
        texSize.y = 10;
        vec3 texture2[texSize.x][texSize.y];
        //  = {
        //     {w, b, w, b, w, b, w, b, w, b},
        //     {b, w, b, w, b, w, b, w, b, w},
        //     {w, b, w, b, w, b, w, b, w, b},
        //     {b, w, b, w, b, w, b, w, b, w},
        //     {w, b, w, b, w, b, w, b, w, b},
        //     {b, w, b, w, b, w, b, w, b, w},
        //     {w, b, w, b, w, b, w, b, w, b},
        //     {b, w, b, w, b, w, b, w, b, w},
        //     {w, b, w, b, w, b, w, b, w, b},
        //     {b, w, b, w, b, w, b, w, b, w}
        // };
        for (int y = 0; y < texSize.y; y++) {
            for (int x = 0; x < texSize.x; x++) {
                texture2[x][y] = (x + y) % 2 ? w : b;
            }
        }

        texture = (vec3**)malloc(sizeof(vec3*) * texSize.y);
        for (int i = 0; i < texSize.y; i++) {
            texture[i] = (vec3*)malloc(sizeof(vec3) * texSize.x);
            for (int j = 0; j < texSize.x; j++) {
                texture[i][j] = vec3(texture2[i][j]);
            }
        }
        ComputeNormal();
    }

    void ComputeNormal()
    {
        glm::vec3 e1 = v1 - v0;
        glm::vec3 e2 = v2 - v0;
        normal = glm::normalize( glm::cross( e2, e1 ) );
    }
};

// Loads the Cornell Box. It is scaled to fill the volume:
// -1 <= x <= +1
// -1 <= y <= +1
// -1 <= z <= +1
void LoadTestModel( std::vector<Triangle>& triangles )
{
    using glm::vec3;

    // Defines colors:
    vec3 red(    0.75f, 0.15f, 0.15f );
    vec3 yellow( 0.75f, 0.75f, 0.15f );
    vec3 green(  0.15f, 0.75f, 0.15f );
    vec3 cyan(   0.15f, 0.75f, 0.75f );
    vec3 blue(   0.15f, 0.15f, 0.75f );
    vec3 purple( 0.75f, 0.15f, 0.75f );
    vec3 white(  0.75f, 0.75f, 0.75f );

    triangles.clear();
    triangles.reserve( 5 * 2 * 3 );

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
    triangles.push_back( Triangle( C, B, A, green, tr, bl, br, 0) );
    triangles.push_back( Triangle( C, D, B, green, tr, tl, bl, 0) );

    // Left wall
    triangles.push_back( Triangle( A, E, C, purple, br, tr, bl, 1) );
    triangles.push_back( Triangle( C, E, G, purple, bl, tr, tl, 1) );

    // Right wall
    triangles.push_back( Triangle( F, B, D, yellow, tl, bl, br, 2) );
    triangles.push_back( Triangle( H, F, D, yellow, tr, tl, br, 2) );

    // Ceiling
    triangles.push_back( Triangle( E, F, G, cyan, tr, tl, br, 3) );
    triangles.push_back( Triangle( F, H, G, cyan, tl, bl, br, 3) );

    // Back wall
    triangles.push_back( Triangle( G, D, C, white, tr, bl, br, 4) );
    triangles.push_back( Triangle( G, H, D, white, tr, tl, bl, 4) );

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

    // Front
    triangles.push_back( Triangle(E, B, A, red, tl, br, bl, 5));
    triangles.push_back( Triangle(E, F, B, red, tl, tr, br, 5));

    // Front
    triangles.push_back( Triangle(F, D, B, red, tr, bl, br, 6) );
    triangles.push_back( Triangle(F, H, D, red, tr, tl, bl, 6) );

    // BACK
    triangles.push_back( Triangle(H, C, D, red, tr, bl, br, 7) );
    triangles.push_back( Triangle(H, G, C, red, tr, tl, bl, 7) );

    // LEFT
    triangles.push_back( Triangle(G, E, C, red, tr, tl, br, 8) );
    triangles.push_back( Triangle(E, A, C, red, tl, bl, br, 8) );

    // TOP
    triangles.push_back( Triangle(G, F, E, red, bl, tr, tl, 9) );
    triangles.push_back( Triangle(G, H, F, red, bl, br, tr, 9) );

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
    triangles.push_back( Triangle(E, B, A, blue, tl, br, bl, 10) );
    triangles.push_back( Triangle(E, F, B, blue, tl, tr, br, 10) );

    // RIGHT
    triangles.push_back( Triangle(F, D, B, blue, tr, bl, br, 11) );
    triangles.push_back( Triangle(F, H, D, blue, tr, tl, bl, 11) );

    // BACK
    triangles.push_back( Triangle(H, C, D, blue, tr, bl, br, 12) );
    triangles.push_back( Triangle(H, G, C, blue, tr, tl, bl, 12) );

    // LEFT
    triangles.push_back( Triangle(G, E, C, blue, tr, tl, br, 13) );
    triangles.push_back( Triangle(E, A, C, blue, tl, bl, br, 13) );

    // TOP
    triangles.push_back( Triangle(G, F, E, blue, bl, tr, tl, 14) );
    triangles.push_back( Triangle(G, H, F, blue, bl, br, tr, 14) );


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

#endif