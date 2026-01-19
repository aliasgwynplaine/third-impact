#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include "SimulationNBodyOptim.hpp"
#include "third_impact_macros.hpp"

SimulationNBodyOptim::SimulationNBodyOptim(const unsigned long nBodies, const std::string &scheme, const float soft,
                                           const unsigned long randInit)
    : SimulationNBodyInterface(nBodies, scheme, soft, randInit)
{
    this->flopsPerIte = 20.f * (float)this->getBodies().getN() * (float)this->getBodies().getN();
    this->accelerations.resize(this->getBodies().getN());
}

void SimulationNBodyOptim::initIteration()
{
    for (unsigned long iBody = 0; iBody < this->getBodies().getN(); iBody++) {
        this->accelerations[iBody].ax = 0.f;
        this->accelerations[iBody].ay = 0.f;
        this->accelerations[iBody].az = 0.f;
    }
}

void SimulationNBodyOptim::computeBodiesAcceleration()
{
    const std::vector<dataAoS_t<float>> &d = this->getBodies().getDataAoS();
    const float softSquared = SQUARE(this->soft);
    const float lG = this->G;
    const unsigned long N = this->getBodies().getN();
    std::vector<accAoS_t<float>> &laccelerations = this->accelerations;

    for (unsigned long iBody = 0; iBody < N; iBody++) {
        float im = d[iBody].m;
        float iqx = d[iBody].qx;
        float iqy = d[iBody].qy;
        float iqz = d[iBody].qz;

        for (unsigned long jBody = iBody + 1; jBody < N; jBody++) {
            float rijx = d[jBody].qx - iqx; // 1 flop
            float rijy = d[jBody].qy - iqy; // 1 flop
            float rijz = d[jBody].qz - iqz; // 1 flop

            // compute the || rij ||² distance between body i and body j
            float rijSquared_softSquared = SQUARE(rijx) + SQUARE(rijy) + SQUARE(rijz) + softSquared; // 5 flops
            // compute the acceleration value between body i and body j: || ai || = G.mj / (|| rij ||² + e²)^{3/2}
            float Gxinvrps = lG / std::pow(rijSquared_softSquared, 3.f / 2.f);
            float ai = Gxinvrps * d[jBody].m;
            float aj = Gxinvrps * im;

            // add the acceleration value into the acceleration vector: ai += || ai ||.rij
            laccelerations[iBody].ax += ai * rijx; // 2 flops
            laccelerations[iBody].ay += ai * rijy; // 2 flops
            laccelerations[iBody].az += ai * rijz; // 2 flops

            // apply Newton's third law. Thanks Ivan :D
            laccelerations[jBody].ax -= aj * rijx; // 2 flops
            laccelerations[jBody].ay -= aj * rijy; // 2 flops
            laccelerations[jBody].az -= aj * rijz; // 2 flops

        }
    }
}

void SimulationNBodyOptim::computeOneIteration()
{
    this->initIteration();
    this->computeBodiesAcceleration();
    // time integration
    this->bodies.updatePositionsAndVelocities(this->accelerations, this->dt);
}
