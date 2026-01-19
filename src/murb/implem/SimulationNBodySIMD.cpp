#include <mipp.h>

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include "SimulationNBodySIMD.hpp"

SimulationNBodySIMD::SimulationNBodySIMD(const unsigned long nBodies, const std::string &scheme, const float soft,
                                           const unsigned long randInit)
    : SimulationNBodyInterface(nBodies, scheme, soft, randInit)
{
    this->flopsPerIte = 20.f * (float)this->getBodies().getN() * (float)this->getBodies().getN();
}

void SimulationNBodySIMD::initIteration()
{
    std::fill(this->vAccelerations.ax.begin(), this->vAccelerations.ax.begin(), 0.f);
    
    std::fill(this->vAccelerations.ay.begin(), this->vAccelerations.ay.begin(), 0.f);
    std::fill(this->vAccelerations.az.begin(), this->vAccelerations.az.begin(), 0.f);
}

void SimulationNBodySIMD::computeBodiesAcceleration()
{
    const dataSoA_t<float> &d = this->getBodies().getDataSoA();

    for (unsigned long iBody = 0; iBody < this->getBodies().getN(); iBody += this->floatN) {
        
	for(unsigned long jBody = 0; jBody < this->getBodies().getN(); jBody++){

		//load and calculate distance for qx	
		mipp::Reg<float> vqxi = d.qx[iBody];
		mipp::Reg<float> vqxj = d.qx[jBody];
		mipp::Reg<float> vqx = vqxj - vqxi;

		//load and calculate distance for qy
		mipp::Reg<float> vqyi = d.qy[iBody];
		mipp::Reg<float> vqyj = d.qy[jBody];
		mipp::Reg<float> vqy = vqyj - vqyi;
		
		//load and calculate distance for qz
		mipp::Reg<float> vqzi = d.qz[iBody];
		mipp::Reg<float> vqzj = d.qz[jBody];
		mipp::Reg<float> vqz = vqzj - vqzi;
	
		mipp::Reg<float> vRijSquared = vqx*vqx + vqy*vqy + vqz*vqz;

		mipp::Reg<float> vSoftSquared = this->soft * this->soft;
	

		// (1/sqrt(vRijSquared + vSoftSquared))^(3) = (vRijSquared + vSoftSquared)^(3/2)
		mipp::Reg<float> rsqrt_cubed = rsqrt(vRijSquared + vSoftSquared);
		rsqrt_cubed = rsqrt_cubed * rsqrt_cubed * rsqrt_cubed;

		mipp::Reg<float> vAi = this->G + d.m[jBody];
	       	vAi = vAi / rsqrt_cubed;	
	
		// this->vAccelerations.ax[iBody].store(vAi * vqx);
		// this->vAccelerations.ay[iBody].store(vAi * vqy);
		// this->vAccelerations.az[iBody].store(vAi * vqz);
	
		std::vector<float> tmp;
		mipp::Reg<float> acceltemp = vAi*vqx;

		acceltemp.store(&tmp);
		this->vAccelerations.ax[iBody] += acceltemp;
	
		
		acceltemp = vAi*vqy;

		acceltemp.store(&tmp);
		this->vAccelerations.ay[iBody] += acceltemp;

			
		acceltemp = vAi*vqz;

		acceltemp.store(&tmp);
		this->vAccelerations.az[iBody] += acceltemp;
	}
    
    }
}

void SimulationNBodySIMD::computeOneIteration()
{
    this->initIteration();
    this->computeBodiesAcceleration();
    // time integration
    this->bodies.updatePositionsAndVelocities(this->vAccelerations, this->dt);
}
