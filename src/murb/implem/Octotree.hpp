#ifndef THIRD_IMPACT_OCTOTREE_HPP_
#define THIRD_IMPACT_OCTOTREE_HPP_

#include <vector>
#include "third_impact_macros.hpp"
#include "core/Bodies.hpp"

template <typename T> class Octotree {
protected:
    T cx;                        /*  coord x - center of octant */
    T cy;                        /*  coord y - center of octant */
    T cz;                        /*  coord z - center of octant */
    T xmin;                      /*      coord x - lower corner */
    T ymin;                      /*      coord y - lower corner */
    T zmin;                      /*      coord z - lower corner */
    T xmax;                      /*      coord x - upper corner */
    T ymax;                      /*      coord y - upper corner */
    T zmax;                      /*      coord z - upper corner */
    T size;                      /*             size of the box */
    T theta;                     /*                   tolerance */
    T cmx;                       /*    coord x - center of mass */
    T cmy;                       /*    coord y - center of mass */
    T cmz;                       /*    coord z - center of mass */
    T total_m;                   /*      total mass of the node */
    bool empty;                  /*      true if noeud is empty */
    int depth;                   /*     depth of the sub-octree */
    int n;                       /*            number of bodies */
    int body;                    /*            body in the leaf */
    T x;                         /* cache for the body position */
    T y;                         /* cache for the body position */
    T z;                         /* cache for the body position */
    T m;                         /*     cache for the body mass */
    Octotree<T> *child[2][2][2]; /*              children x y z */

public :
    Octotree(const T xmin, const T ymin, const T zmin, const T xmax, const T ymax, const T zmax, const int depth, const T theta) /* constructor */
    : xmin(xmin), ymin(ymin), zmin(zmin), xmax(xmax), ymax(ymax), zmax(zmax), depth(depth), theta(theta)
    {
        this->cx = (xmax + xmin) / 2;
        this->cy = (ymax + ymin) / 2;
        this->cz = (zmax + zmin) / 2;
        this->empty = true;
        this->total_m = -1;
        this->n = 0;

        T dx = xmax - xmin;
        T dy = ymax - ymin;
        T dz = zmax - zmin;
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
        child[0][0][0] = new Octotree(xmin, ymin, zmin, cx, cy, cz, depth + 1, theta);
        child[0][1][0] = new Octotree(xmin, cy, zmin, cx, ymax, cz, depth + 1, theta);
        child[1][0][0] = new Octotree(cx, ymin, zmin, xmax, cy, cz, depth + 1, theta);
        child[1][1][0] = new Octotree(cx, cy, zmin, xmax, ymax, cz, depth + 1, theta);
        
        child[0][0][1] = new Octotree(xmin, ymin, cz, cx, cy, zmax, depth + 1, theta);
        child[0][1][1] = new Octotree(xmin, cy, cz, cx, ymin, zmax, depth + 1, theta);
        child[1][0][1] = new Octotree(cx, ymin, cz, xmax, cy, zmax, depth + 1, theta);
        child[1][1][1] = new Octotree(cx, cy, cz, xmax, ymax, zmax, depth + 1, theta);

    }
    void insert(const unsigned long &iBody, const T &im, const T &qix, const T &qiy, const T &qiz)
    {   
        if (n >= 1) {
            int xidx = qix > cx;
            int yidx = qiy > cy;
            int zidx = qiz > cz;

            if (n == 1) {
                create_children();
                child[xidx][yidx][zidx]->insert(iBody, im, qix, qiy, qiz);

                xidx = x > cx;
                yidx = y > cy;
                zidx = z > cz;

                child[xidx][yidx][zidx]->insert(body, m, x, y, z);
            } else {
                child[xidx][yidx][zidx]->insert(iBody, im, qix, qiy, qiz);
            }
        } else {
            empty = false;
            body = iBody;
            x = qix;
            y = qiy;
            z = qiz;
            m = im;
        }

        n++;
    }


    void computeMass()
    {
        if (empty) return;

        if (n == 1) {
            cmx = x;
            cmy = y;
            cmz = z;
            total_m = m;
        } else {
            cmx = cmy = cmz = total_m = 0;

            // TODO: CHECK UNROLLING 
            for (int i = 0; i < 2; i++){
                for (int j = 0; j < 2; j++){
                    for (int k = 0; k < 2; k++) {
                        Octotree *c = child[i][j][k];
                        c->computeMass();
                        float cmass = c->getTotalMass();

                        if (cmass < 0) continue;

                        cmx += c->getCMX() * cmass;
                        cmy += c->getCMY() * cmass;
                        cmz += c->getCMZ() * cmass;
                        total_m += cmass;
                    }
                }
            }

            cmx /= total_m;
            cmy /= total_m;
            cmz /= total_m;
        }
    }

    void computeAcc(const T &ix, const T &iy, const T &iz, const T &soft, const float &G, accAoS_t<T> &acc)
    {
        if (empty) return;

        const float rx = cmx - ix;
        const float ry = cmy - iy;
        const float rz = cmz - iz;

        const float r_squared = rx * rx + ry * ry + rz * rz;
        const float soft_squared = soft * soft;

        if (n == 1) {
            const float ai = G * total_m * POW3(FAST_RSQRT(r_squared + soft_squared));

            acc.ax += ai * rx;
            acc.ay += ai * ry;
            acc.az += ai * rz;
        } else {
            const float r = r_squared * FAST_RSQRT(r_squared);

            if (size / r < theta) {
                const float ai = G * total_m * POW3(FAST_RSQRT(r_squared + soft_squared));
                
                acc.ax += ai * rx;
                acc.ay += ai * ry;
                acc.az += ai * rz;
            } else {
                for (int i = 0; i < 2; i++){
                    for (int j = 0; j < 2; j++){
                        for (int k = 0; k < 2; k++) {
                            child[i][j][k]->computeAcc(ix, iy, iz, soft, G, acc);
                        }
                    }
                }
            }
        }
    }
    /* getters: encapsulation in oop is nonsense
     * just a pedantic feature that overcomplicates
     * software scalation.
     * correcteness of the software is programmer's 
     * responsability.
     */
    int getDepth() const { return depth; }

    int getN() const { return n; }

    bool isEmpty() const { return empty; }

    T getTotalMass() const { return total_m; }

    T getCMX() const { return cmx; }

    T getCMY() const { return cmy; }

    T getCMZ() const { return cmz; }
};

#endif