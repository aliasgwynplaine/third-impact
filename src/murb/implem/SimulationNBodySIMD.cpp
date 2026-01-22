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

    //TODO : need to allocate the ax ay and az
}

void SimulationNBodySIMD::initIteration()
{
    std::fill(this->vAccelerations.ax.begin(), this->vAccelerations.ax.end(), 0.f); 
    std::fill(this->vAccelerations.ay.begin(), this->vAccelerations.ay.end(), 0.f);
    std::fill(this->vAccelerations.az.begin(), this->vAccelerations.az.begin(), 0.f);
}

void SimulationNBodySIMD::computeBodiesAcceleration()
{
    const dataSoA_t<float> &d = this->getBodies().getDataSoA();

    mipp::Reg<float> vG = this->G;

    mipp::Reg<float> vSoftSquared = this->soft * this->soft;
	

    //sum var    
    float qx_sum = 0.f, qy_sum = 0.f, qz_sum = 0.f;

    for (unsigned long iBody = 0; iBody < this->getBodies().getN(); iBody += this->floatN) {
       
	//load i positions 
	mipp::Reg<float> vqxi = d.qx[iBody];
	mipp::Reg<float> vqyi = d.qy[iBody];
	mipp::Reg<float> vqzi = d.qz[iBody];
	
	for(unsigned long jBody = 0; jBody < this->getBodies().getN(); jBody++){

		//load j positions and calculate distance	
		mipp::Reg<float> vqxj = d.qx[jBody];
		mipp::Reg<float> vqx = vqxj - vqxi;

		mipp::Reg<float> vqyj = d.qy[jBody];
		mipp::Reg<float> vqy = vqyj - vqyi;
		
		mipp::Reg<float> vqzj = d.qz[jBody];
		mipp::Reg<float> vqz = vqzj - vqzi;
		
		//load j mass
		mipp::Reg<float> vmj = d.m[jBody];

		
		mipp::Reg<float> vRijSquared = vqx*vqx + vqy*vqy + vqz*vqz;


		// (1/sqrt(vRijSquared + vSoftSquared))^(3) = (vRijSquared + vSoftSquared)^(3/2)
		mipp::Reg<float> rsqrt_cubed = rsqrt(vRijSquared + vSoftSquared);
		rsqrt_cubed = rsqrt_cubed * rsqrt_cubed * rsqrt_cubed;

		mipp::Reg<float> vAi = vG * vmj / rsqrt_cubed;
	
		mipp::Reg<float> vAccx = vAi * vqx;
		mipp::Reg<float> vAccy = vAi * vqy;
		mipp::Reg<float> vAccz = vAi * vqz;

		qx_sum = mipp::sum(vAccx);
		qy_sum = mipp::sum(vAccy);
		qz_sum = mipp::sum(vAccz);
	
	}
   
	this->vAccelerations.ax[iBody] = qx_sum;
	this->vAccelerations.ay[iBody] = qy_sum;
	this->vAccelerations.az[iBody] = qz_sum;

    }
}

void SimulationNBodySIMD::computeOneIteration()
{
    this->initIteration();
    this->computeBodiesAcceleration();
    // time integration
    this->bodies.updatePositionsAndVelocities(this->vAccelerations, this->dt);
}
