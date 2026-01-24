#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include "SimulationNBodyOpenMP.hpp"
#include "third_impact_macros.hpp"

SimulationNBodyOpenMP::SimulationNBodyOpenMP(const unsigned long nBodies, const std::string &scheme, const float soft,
                                           const unsigned long randInit)
    : SimulationNBodyInterface(nBodies, scheme, soft, randInit)
{
    float N = (float)this->getBodies().getN();
    this->flopsPerIte = 0.5f * 27.f * ((N - 1) * (N - 2));
    this->accelerations.resize(this->getBodies().getN());
    initIteration();
}

void SimulationNBodyOpenMP::initIteration()
{
    unsigned long N = this->getBodies().getN();

    #pragma omp parallel for firstprivate(N)
    for (unsigned long iBody = 0; iBody < N; iBody++) {
        this->accelerations[iBody].ax = 0.f;
        this->accelerations[iBody].ay = 0.f;
        this->accelerations[iBody].az = 0.f;
    }
}

void SimulationNBodyOpenMP::computeBodiesAcceleration()
{
    const std::vector<dataAoS_t<float>> &d = this->getBodies().getDataAoS();
    const float softSquared = SQUARE(this->soft); // 1 flop
    const float lG = this->G;
    const unsigned long N = this->getBodies().getN();
    std::vector<accAoS_t<float>> &laccelerations = this->accelerations;

    #pragma omp parallel
    {
    accAoS_t<float> *acc_priv = (accAoS_t<float> *)calloc(N, sizeof(accAoS_t<float>));

    #pragma omp for
    for (unsigned long iBody = 0; iBody < N; iBody++) {
        float im  = d[iBody].m;
        float iqx = d[iBody].qx;
        float iqy = d[iBody].qy;
        float iqz = d[iBody].qz;
        accAoS_t<float> acc = {0.f, 0.f, 0.f};

        for (unsigned long jBody = iBody + 1; jBody < N; jBody++) {
            float rijx = d[jBody].qx - iqx; // 1 flop
            float rijy = d[jBody].qy - iqy; // 1 flop
            float rijz = d[jBody].qz - iqz; // 1 flop

            // compute the || rij ||² + e²
            float rijSquared_softSquared = SQUARE(rijx) + SQUARE(rijy) + SQUARE(rijz) + softSquared; // 7 flops
            // compute G / (|| rij ||² + e²)^{3/2}
            float Gxinvrps = lG * POW3(FAST_RSQRT(rijSquared_softSquared)); // 2 flops
            // compute the acceleration value between body i and body j: || ai || = G.mj / (|| rij ||² + e²)^{3/2}
            float ai = Gxinvrps * d[jBody].m; // 1 flops
            // compute the acceleration value between body j and body i: || aj || = G.mi / (|| rij ||² + e²)^{3/2}
            float aj = Gxinvrps * im; // 1 flops

            // add the acceleration value into the acceleration vector: ai += || ai ||.rij
            acc.ax += ai * rijx; // 2 flops
            acc.ay += ai * rijy; // 2 flops
            acc.az += ai * rijz; // 2 flops

            acc_priv[jBody].ax -= aj * rijx; // 2 flops
            acc_priv[jBody].ay -= aj * rijy; // 2 flops
            acc_priv[jBody].az -= aj * rijz; // 2 flops
        }

        acc_priv[iBody].ax += acc.ax;
        acc_priv[iBody].ay += acc.ay;
        acc_priv[iBody].az += acc.az;
    }

    #pragma omp critical
    for (unsigned long i = 0; i < N; i++) {
        laccelerations[i].ax += acc_priv[i].ax;
        laccelerations[i].ay += acc_priv[i].ay;
        laccelerations[i].az += acc_priv[i].az;
    }

    free(acc_priv);
    }
}

void SimulationNBodyOpenMP::updatePositionsAndVelocities() {
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

    #pragma omp parallel for firstprivate(n, dt)
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

void SimulationNBodyOpenMP::computeOneIteration()
{
    //this->initIteration();
    this->computeBodiesAcceleration();
    // time integration
    //this->bodies.updatePositionsAndVelocities(this->accelerations, this->dt);
    this->updatePositionsAndVelocities();
}
