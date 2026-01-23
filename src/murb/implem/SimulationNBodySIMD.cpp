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

    //TODO : have a better size..
    this->vAccelerations.ax.resize(this->getBodies().getN() + this->floatN);
    this->vAccelerations.ay.resize(this->getBodies().getN() + this->floatN);
    this->vAccelerations.az.resize(this->getBodies().getN() + this->floatN);

}

void SimulationNBodySIMD::initIteration()
{
    std::fill(this->vAccelerations.ax.begin(), this->vAccelerations.ax.end(), 0.f); 
    std::fill(this->vAccelerations.ay.begin(), this->vAccelerations.ay.end(), 0.f);
    std::fill(this->vAccelerations.az.begin(), this->vAccelerations.az.end(), 0.f);
}

void SimulationNBodySIMD::computeBodiesAcceleration()
{
    const dataSoA_t<float> &d = this->getBodies().getDataSoA();
	

    //const in Reg
    mipp::Reg<float> vG = this->G;
    mipp::Reg<float> vSoftSquared = this->soft * this->soft;
	
    //init Registers here for readability and optimisation?
    mipp::Reg<float> vqxi, vqxj, vqx;
    mipp::Reg<float> vqyi, vqyj, vqy;
    mipp::Reg<float> vqzi, vqzj, vqz;
    mipp::Reg<float> vmj;
    mipp::Reg<float> vRijSquared, rsqrt_cubed, vAi;
    
    
    for (unsigned long iBody = 0; iBody < this->getBodies().getN(); iBody += this->floatN) {
       
	//loadu i positions 
	vqxi.loadu(&d.qx[iBody]);	
	vqyi.loadu(&d.qy[iBody]);	
	vqzi.loadu(&d.qz[iBody]);	
	    
	//vqxi = d.qx[iBody];
	//vqyi = d.qy[iBody];
	//vqzi = d.qz[iBody];
	
    	mipp::Reg<float> vAccx = 0.f, vAccy = 0.f, vAccz = 0.f;
	
	for(unsigned long jBody = 0; jBody < this->getBodies().getN(); jBody ++){

		//loadu j positions and calculate distance	
		
		//vqxj = d.qx[jBody];
		vqxj = d.qx[jBody];
		vqx = vqxj - vqxi;

		//vqyj = d.qy[jBody];
		vqyj = d.qy[jBody];
		vqy = vqyj - vqyi;
		
		//vqzj = d.qz[jBody];
		vqzj = d.qz[jBody];
		vqz = vqzj - vqzi;
		
		//loadu j mass
		//vmj = d.m[jBody];
		vmj = d.m[jBody];
		
		vRijSquared = vqx*vqx + vqy*vqy + vqz*vqz;


		// (1/sqrt(vRijSquared + vSoftSquared))^(3) = (vRijSquared + vSoftSquared)^(3/2)
		rsqrt_cubed = mipp::rsqrt(vRijSquared + vSoftSquared);
		rsqrt_cubed = rsqrt_cubed * rsqrt_cubed * rsqrt_cubed;

		vAi = vG * vmj * rsqrt_cubed;
	
		vAccx = vAi * vqx + vAccx;
		vAccy = vAi * vqy + vAccy;
		vAccz = vAi * vqz + vAccz;

	}
  	

	mipp::storeu<float>(&this->vAccelerations.ax[iBody], vAccx);
	mipp::storeu<float>(&this->vAccelerations.ay[iBody], vAccy);
	mipp::storeu<float>(&this->vAccelerations.az[iBody], vAccz);

    }
}

void SimulationNBodySIMD::computeOneIteration()
{
    this->initIteration();
    this->computeBodiesAcceleration();
    // time integration
    this->bodies.updatePositionsAndVelocities(this->vAccelerations, this->dt);
}
