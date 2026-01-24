#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include "SimulationNBodyBarnesHut.hpp"

SimulationNBodyBarnesHut::SimulationNBodyBarnesHut(const unsigned long nBodies, const std::string &scheme, const float soft,
                                           const unsigned long randInit, const float theta)
    : SimulationNBodyInterface(nBodies, scheme, soft, randInit), theta(theta)
{
    this->flopsPerIte = 20.f * (float)this->getBodies().getN() * (float)this->getBodies().getN();
    this->accelerations.resize(this->getBodies().getN());
}

void SimulationNBodyBarnesHut::initIteration()
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

void SimulationNBodyBarnesHut::computeBodiesAcceleration()
{
    const std::vector<dataAoS_t<float>> &d = this->getBodies().getDataAoS();
    const unsigned long N = this->getBodies().getN();
    
    for (unsigned long iBody = 0; iBody < N; iBody++) {
        this->root->insert(iBody, d[iBody].m, d[iBody].qx, d[iBody].qy, d[iBody].qz);
    }

    this->root->computeCM();

    for (unsigned long iBody = 0; iBody < N; iBody++) {
        this->root->computeAcc(d[iBody].qx, d[iBody].qy, d[iBody].qz, soft, G, this->accelerations[iBody]);
    }
}

void SimulationNBodyBarnesHut::computeOneIteration()
{
    this->initIteration();
    this->computeBodiesAcceleration();
    // time integration
    this->bodies.updatePositionsAndVelocities(this->accelerations, this->dt);
    delete this->root;
}
