#ifndef THIRD_IMPACT_OCTOTREE_HPP_
#define THIRD_IMPACT_OCTOTREE_HPP_

#include <vector>
#include "third_impact_macros.hpp"
#include "core/Bodies.hpp"
#include <iostream>

typedef struct Octotree Octotree;

// this struct is for the stack
struct octostack {
    Octotree *d;
    bool v;
};

struct Octotree {
    float cx, cy, cz;          /*   coords - center of octant */
    float xmin, ymin, zmin;    /*       coords - lower corner */
    float xmax, ymax, zmax;    /*       coords - upper corner */
    float size;                /*             size of the box */
    float theta;               /*                   tolerance */
    float cmx, cmy, cmz;       /*     coords - center of mass */
    float m;                   /*       total mass in the box */
    int n;                     /* number of bodies in the box */
    int body;                  /*            body in the leaf */
    struct Octotree *child[2][2][2];  /*       child in octant 0bxyz */

    Octotree(const float xmin, const float ymin, const float zmin, 
             const float xmax, const float ymax, const float zmax, 
             const float theta) /* constructor */
    : xmin(xmin), ymin(ymin), zmin(zmin), 
      xmax(xmax), ymax(ymax), zmax(zmax), theta(theta)
    {
        this->cx = (xmax + xmin) * 0.5f;
        this->cy = (ymax + ymin) * 0.5f;
        this->cz = (zmax + zmin) * 0.5f;
        this->m = 0;
        this->n = 0;

        float dx = xmax - xmin;
        float dy = ymax - ymin;
        float dz = zmax - zmin;
        this->size = (dx * dx + dy * dy + dz * dz) / 1.732050808f; // sqrt(3)
    }

    ~Octotree()
    {
        if (n > 1) {
            destroy_children();
        }
    }

    void destroy_children()
    {
        delete child[0][0][0];
        delete child[0][1][0];
        delete child[1][0][0];
        delete child[1][1][0];

        delete child[0][0][1];
        delete child[0][1][1];
        delete child[1][0][1];
        delete child[1][1][1];
    }

    void create_children()
    {
        child[0][0][0] = new Octotree(xmin, ymin, zmin, cx, cy, cz, theta);
        child[0][1][0] = new Octotree(xmin, cy, zmin, cx, ymax, cz, theta);
        child[1][0][0] = new Octotree(cx, ymin, zmin, xmax, cy, cz, theta);
        child[1][1][0] = new Octotree(cx, cy, zmin, xmax, ymax, cz, theta);
        
        child[0][0][1] = new Octotree(xmin, ymin, cz, cx, cy, zmax, theta);
        child[0][1][1] = new Octotree(xmin, cy, cz, cx, ymin, zmax, theta);
        child[1][0][1] = new Octotree(cx, ymin, cz, xmax, cy, zmax, theta);
        child[1][1][1] = new Octotree(cx, cy, cz, xmax, ymax, zmax, theta);

    }
    void insert(const unsigned long &iBody, const float &im, 
        const float &qix, const float &qiy, const float &qiz)
    {   
        if (n >= 1) {
            int xidx = qix > cx;
            int yidx = qiy > cy;
            int zidx = qiz > cz;

            if (n == 1) {
                create_children();
                child[xidx][yidx][zidx]->insert(iBody, im, qix, qiy, qiz);

                xidx = cmx > cx;
                yidx = cmy > cy;
                zidx = cmz > cz;
                
                child[xidx][yidx][zidx]->insert(body, m, cmx, cmy, cmz);
                m += im;
            } else {
                m += im;
                child[xidx][yidx][zidx]->insert(iBody, im, qix, qiy, qiz);
            }
        } else {
            body = iBody;
            cmx = qix;
            cmy = qiy;
            cmz = qiz;
            m = im;
        }

        n++;
    }


    int computeCM()
    {
        if (n == 0) return 0;

        struct octostack *stack = new struct octostack[n];
        int h = 0; // height of the stack

        stack[h] = {this, false};

        struct octostack curr;

        int flops = 0;

        while (h > 0) {
            curr = stack[--h];
            Octotree *node = curr.d;

            if (node->n == 1) continue;

            if (!curr.v) {
                curr.v = true;
                stack[h++] = curr;

                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        for (int k = 0; k < 2; k++) {
                            Octotree *c = node->child[i][j][k];

                            if (c->n > 0) 
                                stack[h++] = {c, false};
                        }
                    }
                }
            } else {
                node->cmx = 0.f;
                node->cmy = 0.f;
                node->cmz = 0.f;

                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        for (int k = 0; k < 2; k++) {
                            Octotree *c = node->child[i][j][k];

                            if (c->n == 0) continue;
                            
                            node->cmx += c->cmx * c->m; // 2 flops
                            node->cmy += c->cmy * c->m; // 2 flops
                            node->cmz += c->cmz * c->m; // 2 flops

                            flops += 6;
                        }
                    }
                }

                float inv_m = 1.f / node->m; // 1 flops
                node->cmx *= inv_m;          // 1 flops
                node->cmy *= inv_m;          // 1 flops
                node->cmz *= inv_m;          // 1 flops

                flops += 4;
            }
        }
        delete[] stack;
        
        return flops;
    }

    int computeAcc(const float &ix, const float &iy, const float &iz, 
        const float &soft, const float &G, accAoS_t<float> &acc)
    {
        if (n == 0) return 0;

        Octotree **stack = new Octotree*[n];
        int h = 0;
        stack[h++] = this;

        int flops = 0;

        while (h > 0) {
            Octotree *node = stack[--h];

            if (node->n == 0) continue;

            const float rx = node->cmx - ix;  // 1 flops
            const float ry = node->cmy - iy;  // 1 flops
            const float rz = node->cmz - iz;  // 1 flops

            const float r_squared = rx * rx + ry * ry + rz * rz; // 5 flops
            const float soft_squared = soft * soft; // 1 flop

            flops += 9;

            if (node->n == 1) {
                const float ai = G * node->m * POW3(FAST_RSQRT(r_squared + soft_squared)); // 5 flop
                acc.ax += ai * rx; // 2 flops
                acc.ay += ai * ry; // 2 flops
                acc.az += ai * rz; // 2 flops
                flops += 11;
            } else {
                const float r = r_squared * FAST_RSQRT(r_squared); // 2 flops
    
                if (node->size / r < node->theta) {
                    const float ai = G * node->m * POW3(FAST_RSQRT(r_squared + soft_squared)); // 5 flops
                    
                    acc.ax += ai * rx; // 2 flops
                    acc.ay += ai * ry; // 2 flops
                    acc.az += ai * rz; // 2 flops
                    flops += 13;
                } else {
                    for (int i = 0; i < 2; i++){
                        for (int j = 0; j < 2; j++){
                            for (int k = 0; k < 2; k++) {
                                Octotree *c = node->child[i][j][k];

                                if (c->n > 0)
                                    stack[h++] = c;
                            }
                        }
                    }
                }
            }
        }

        delete[] stack;

        return flops;
    }


    int computeAcc(const float &ix, const float &iy, const float &iz, 
        const float &soft, const float &G, float &iax, float &iay, float &iaz)
    {
        if (n == 0) return 0;

        Octotree **stack = new Octotree*[n];
        int h = 0;
        stack[h++] = this;

        int flops = 0;

        while (h > 0) {
            Octotree *node = stack[--h];

            if (node->n == 0) continue;

            const float rx = node->cmx - ix;  // 1 flops
            const float ry = node->cmy - iy;  // 1 flops
            const float rz = node->cmz - iz;  // 1 flops

            const float r_squared = rx * rx + ry * ry + rz * rz; // 5 flops
            const float soft_squared = soft * soft; // 1 flop

            flops += 9;

            if (node->n == 1) {
                const float ai = G * node->m * POW3(FAST_RSQRT(r_squared + soft_squared)); // 5 flop
                iax += ai * rx; // 2 flops
                iay += ai * ry; // 2 flops
                iaz += ai * rz; // 2 flops
                flops += 11;
            } else {
                const float r = r_squared * FAST_RSQRT(r_squared); // 2 flops
    
                if (node->size / r < node->theta) {
                    const float ai = G * node->m * POW3(FAST_RSQRT(r_squared + soft_squared)); // 5 flops
                    
                    iax += ai * rx; // 2 flops
                    iay += ai * ry; // 2 flops
                    iaz += ai * rz; // 2 flops
                    flops += 13;
                } else {
                    for (int i = 0; i < 2; i++){
                        for (int j = 0; j < 2; j++){
                            for (int k = 0; k < 2; k++) {
                                Octotree *c = node->child[i][j][k];

                                if (c->n > 0)
                                    stack[h++] = c;
                            }
                        }
                    }
                }
            }
        }

        delete[] stack;

        return flops;
    }
    /* getters: encapsulation in oop is nonsense
     * just a pedantic feature that overcomplicates
     * software scalation.
     * correcteness of the software is programmer's 
     * responsability.
     */

    int getN() const { return n; }

    float getTotalMass() const { return m; }

    float getCMX() const { return cmx; }

    float getCMY() const { return cmy; }

    float getCMZ() const { return cmz; }
};

#endif