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
    this->vAccelerations.ax.resize(this->getBodies().getN());
    this->vAccelerations.ay.resize(this->getBodies().getN());
    this->vAccelerations.az.resize(this->getBodies().getN());

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
	

    //const in Reg
    mipp::Reg<float> vG = this->G;
    mipp::Reg<float> vSoftSquared = this->soft * this->soft;
	
    //sum registers    
    float qx_sum = 0.f, qy_sum = 0.f, qz_sum = 0.f;

    //init Registers here for readability and optimisation?
    mipp::Reg<float> vqxi, vqxj, vqx;
    mipp::Reg<float> vqyi, vqyj, vqy;
    mipp::Reg<float> vqzi, vqzj, vqz;
    mipp::Reg<float> vmj;
    mipp::Reg<float> vRijSquared, rsqrt_cubed, vAi;
    mipp::Reg<float> vAccx, vAccy, vAccz;
    
    
    for (unsigned long iBody = 0; iBody < this->getBodies().getN(); iBody += this->floatN) {
       
	//load i positions 
	vqxi = d.qx[iBody];
	vqyi = d.qy[iBody];
	vqzi = d.qz[iBody];
	
	for(unsigned long jBody = 0; jBody < this->getBodies().getN(); jBody++){

		//load j positions and calculate distance	
		vqxj = d.qx[jBody];
		vqx = vqxj - vqxi;

		vqyj = d.qy[jBody];
		vqy = vqyj - vqyi;
		
		vqzj = d.qz[jBody];
		vqz = vqzj - vqzi;
		
		//load j mass
		vmj = d.m[jBody];

		
		vRijSquared = vqx*vqx + vqy*vqy + vqz*vqz;


		// (1/sqrt(vRijSquared + vSoftSquared))^(3) = (vRijSquared + vSoftSquared)^(3/2)
		rsqrt_cubed = rsqrt(vRijSquared + vSoftSquared);
		rsqrt_cubed = rsqrt_cubed * rsqrt_cubed * rsqrt_cubed;

		vAi = vG * vmj / rsqrt_cubed;
	
		vAccx = vAi * vqx;
		vAccy = vAi * vqy;
		vAccz = vAi * vqz;

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
