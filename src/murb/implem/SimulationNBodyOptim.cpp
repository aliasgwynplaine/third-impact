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
    const long N = this->getBodies().getN();
    std::vector<accAoS_t<float>> &laccelerations = this->accelerations;

    for (unsigned long iBody = 0; iBody < N; iBody++) {
        for (unsigned long jBody = iBody + 1; jBody < N; jBody++) {
            const float rijx = d[jBody].qx - d[iBody].qx; // 1 flop
            const float rijy = d[jBody].qy - d[iBody].qy; // 1 flop
            const float rijz = d[jBody].qz - d[iBody].qz; // 1 flop

            // compute the || rij ||² distance between body i and body j
            const float rijSquared = SQUARE(rijx) + SQUARE(rijy) + SQUARE(rijz); // 5 flops
            // compute e²
            // compute the acceleration value between body i and body j: || ai || = G.mj / (|| rij ||² + e²)^{3/2}
            const float rps = std::pow(rijSquared + softSquared, 3.f / 2.f);
            const float ai = lG * d[jBody].m / rps; // 5 flops
            const float aj = lG * d[iBody].m / rps;

            // add the acceleration value into the acceleration vector: ai += || ai ||.rij
            laccelerations[iBody].ax += ai * rijx; // 2 flops
            laccelerations[iBody].ay += ai * rijy; // 2 flops
            laccelerations[iBody].az += ai * rijz; // 2 flops

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
