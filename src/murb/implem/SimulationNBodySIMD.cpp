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
	
    for (unsigned long iBody = 0; iBody < this->getBodies().getN(); iBody += this->floatN) {
       
	//load i positions 
	const mipp::Reg<float> vqxi = mipp::load<float>(&d.qx[iBody]);	
	const mipp::Reg<float> vqyi = mipp::load<float>(&d.qy[iBody]);	
	const mipp::Reg<float> vqzi = mipp::load<float>(&d.qz[iBody]);	
	    
    	mipp::Reg<float> vAccx = 0.f, vAccy = 0.f, vAccz = 0.f;
	
	for(unsigned long jBody = 0; jBody < this->getBodies().getN(); jBody ++){

		//set j positions in register and calculate distance with i	
		
		const mipp::Reg<float> vqxj = mipp::set1(d.qx[jBody]);
		const mipp::Reg<float> vqx = vqxj - vqxi; //1 flop

		const mipp::Reg<float> vqyj = mipp::set1(d.qy[jBody]);
		const mipp::Reg<float> vqy = vqyj - vqyi; //1 flop
		
		const mipp::Reg<float> vqzj = mipp::set1(d.qz[jBody]);
		const mipp::Reg<float> vqz = vqzj - vqzi; //1 flop
		
		//set j mass in register
		const mipp::Reg<float> vmj = mipp::set1(d.m[jBody]);
		
		const mipp::Reg<float> vRijSquared = vqx*vqx + vqy*vqy + vqz*vqz; //5 flops


		// (1/sqrt(vRijSquared + vSoftSquared))^(3) = (vRijSquared + vSoftSquared)^(3/2)
		mipp::Reg<float> rsqrt_cubed = mipp::rsqrt(vRijSquared + vSoftSquared); //2 flops
		rsqrt_cubed = rsqrt_cubed * rsqrt_cubed * rsqrt_cubed; //2 flops

		const mipp::Reg<float> vAi = vG * vmj * rsqrt_cubed; //2 flops
			
		vAccx = mipp::fmadd(vAi, vqx, vAccx); // 2 flops
		vAccy = mipp::fmadd(vAi, vqy, vAccy); // 2 flops
		vAccz = mipp::fmadd(vAi, vqz, vAccz); // 2 flops

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
