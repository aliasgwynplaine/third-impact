#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include "SimulationNBodyOptim.hpp"

SimulationNBodyOptim::SimulationNBodyOptim(const unsigned long nBodies, const std::string &scheme, const float soft,
                                           const unsigned long randInit, const float theta)
    : SimulationNBodyInterface(nBodies, scheme, soft, randInit), theta(theta)
{
    this->flopsPerIte = 20.f * (float)this->getBodies().getN() * (float)this->getBodies().getN();
    this->accelerations.resize(this->getBodies().getN());
}

void SimulationNBodyOptim::initIteration()
{
    const std::vector<dataAoS_t<float>> &d = this->getBodies().getDataAoS();
    float xmin, ymin, zmin;
    float xmax, ymax, zmax;

    xmin = xmax = d[0].qx;
    ymin = ymax = d[0].qy;
    zmin = zmax = d[0].qz;

    this->accelerations[0].ax = 0.f;
    this->accelerations[0].ay = 0.f;
    this->accelerations[0].az = 0.f;

    for (unsigned long iBody = 1; iBody < this->getBodies().getN(); iBody++) {
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

    //std::cout<<"pmin: "<<xmin<<" "<<ymin<<" "<<zmin<<std::endl;
    //std::cout<<"pmax: "<<xmax<<" "<<ymax<<" "<<zmax<<std::endl;
    this->root = new Octotree<float>(xmin, ymin, zmin, xmax, ymax, zmax, 0, theta);
}

void SimulationNBodyOptim::computeBodiesAcceleration()
{
    const std::vector<dataAoS_t<float>> &d = this->getBodies().getDataAoS();
    
    for (unsigned long iBody = 0; iBody < this->getBodies().getN(); iBody++) {
        this->root->insert(iBody, d[iBody].m, d[iBody].qx, d[iBody].qy, d[iBody].qz);
    }

    this->root->computeMass();

    for (unsigned long iBody = 0; iBody < this->getBodies().getN(); iBody++) {
        this->root->computeAcc(d[iBody].qx, d[iBody].qy, d[iBody].qz, soft, G, this->accelerations[iBody]);
    }
}

void SimulationNBodyOptim::computeOneIteration()
{
    this->initIteration();
    this->computeBodiesAcceleration();
    // time integration
    this->bodies.updatePositionsAndVelocities(this->accelerations, this->dt);
    delete this->root;
}
