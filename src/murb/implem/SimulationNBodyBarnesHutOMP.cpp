#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include "SimulationNBodyBarnesHutOMP.hpp"

SimulationNBodyBarnesHutOMP::SimulationNBodyBarnesHutOMP(const unsigned long nBodies, const std::string &scheme, const float soft,
                                           const unsigned long randInit, const float theta)
    : SimulationNBodyInterface(nBodies, scheme, soft, randInit), theta(theta)
{
    this->flopsPerIte = 20.f * (float)this->getBodies().getN() * std::log((float)this->getBodies().getN());
    this->accelerations.resize(this->getBodies().getN());
}

void SimulationNBodyBarnesHutOMP::initIteration()
{
    const std::vector<dataAoS_t<float>> &d = this->getBodies().getDataAoS();
    const unsigned long N = this->getBodies().getN();
    float xmin, ymin, zmin;
    float xmax, ymax, zmax;

    xmin = xmax = d[0].qx;
    ymin = ymax = d[0].qy;
    zmin = zmax = d[0].qz;

    this->accelerations[0].ax = 0.f;
    this->accelerations[0].ay = 0.f;
    this->accelerations[0].az = 0.f;

    for (unsigned long iBody = 1; iBody < N; iBody++) {
        this->accelerations[iBody].ax = 0.f;
        this->accelerations[iBody].ay = 0.f;
        this->accelerations[iBody].az = 0.f;

        xmin = xmin <= d[iBody].qx ? xmin : d[iBody].qx;
        ymin = ymin <= d[iBody].qy ? ymin : d[iBody].qy;
        zmin = zmin <= d[iBody].qz ? zmin : d[iBody].qz;

        xmax = xmax >= d[iBody].qx ? xmax : d[iBody].qx;
        ymax = ymax >= d[iBody].qy ? ymax : d[iBody].qy;
        zmax = zmax >= d[iBody].qz ? zmax : d[iBody].qz;
    }

    this->root = new Octotree(xmin, ymin, zmin, xmax, ymax, zmax, theta);
}

void SimulationNBodyBarnesHutOMP::computeBodiesAcceleration()
{
    const std::vector<dataAoS_t<float>> &d = this->getBodies().getDataAoS();
    const unsigned long N = this->getBodies().getN();
    
    for (unsigned long iBody = 0; iBody < N; iBody++) {
        this->root->insert(iBody, d[iBody].m, d[iBody].qx, d[iBody].qy, d[iBody].qz);
    }

    Octotree **child = new Octotree*[8];

    child[0] = this->root->child[0][0][0];
    child[1] = this->root->child[0][1][0];
    child[2] = this->root->child[1][0][0];
    child[3] = this->root->child[1][1][0];
    child[4] = this->root->child[0][0][1];
    child[5] = this->root->child[0][1][1];
    child[6] = this->root->child[1][0][1];
    child[7] = this->root->child[1][1][1];

    #pragma omp parallel for schedule(runtime)
    for (int i = 0; i < 8; i++) {
        child[i]->computeCM();
    }

    for (int i = 0; i < 8; i++) {
        this->root->m += child[i]->m;
        this->root->cmx += child[i]->cmx * child[i]->m;
        this->root->cmy += child[i]->cmy * child[i]->m;
        this->root->cmz += child[i]->cmz * child[i]->m;
    }

    float inv_m = 1 / this->root->m;
    this->root->cmx *= inv_m;
    this->root->cmy *= inv_m;
    this->root->cmz *= inv_m;

    #pragma omp parallel for
    for (unsigned long iBody = 0; iBody < N; iBody++) {
        this->root->computeAcc(d[iBody].qx, d[iBody].qy, d[iBody].qz, soft, G, this->accelerations[iBody]);
    }
}

void SimulationNBodyBarnesHutOMP::updatePositionsAndVelocities() {
    unsigned long n = this->getBodies().getN();
    dataAoS_t<float> *daos = const_cast<dataAoS_t<float>*>(this->getBodies().getDataAoS().data());
    const dataSoA_t<float> &dsoa = this->getBodies().getDataSoA();
    float *qx = const_cast<float*>(dsoa.qx.data());
    float *qy = const_cast<float*>(dsoa.qy.data());
    float *qz = const_cast<float*>(dsoa.qz.data());
    float *vx = const_cast<float*>(dsoa.vx.data());
    float *vy = const_cast<float*>(dsoa.vy.data());
    float *vz = const_cast<float*>(dsoa.vz.data());
    std::vector<accAoS_t<float>> &lacc = this->accelerations;

    #pragma omp parallel for schedule(runtime) firstprivate(n, dt)
    for (unsigned long i = 0; i < n; i++) {
        float axdt = lacc[i].ax * dt;
        float aydt = lacc[i].ay * dt;
        float azdt = lacc[i].az * dt;
        qx[i] += (vx[i] + axdt * 0.5) * dt;
        qy[i] += (vy[i] + aydt * 0.5) * dt;
        qz[i] += (vz[i] + azdt * 0.5) * dt;
        vx[i] += axdt;
        vy[i] += aydt;
        vz[i] += azdt;

        daos[i].qx = qx[i];
        daos[i].qy = qy[i];
        daos[i].qz = qz[i];
        daos[i].vx = vx[i];
        daos[i].vy = vy[i];
        daos[i].vz = vz[i];

        // update the acc here so we won't use another loop
        lacc[i].ax = 0.f;
        lacc[i].ay = 0.f;
        lacc[i].az = 0.f;
    }
}

void SimulationNBodyBarnesHutOMP::computeOneIteration()
{
    this->initIteration();
    this->computeBodiesAcceleration();
    // time integration
    delete this->root;
    //this->bodies.updatePositionsAndVelocities(this->accelerations, this->dt);
    this->updatePositionsAndVelocities();
}
