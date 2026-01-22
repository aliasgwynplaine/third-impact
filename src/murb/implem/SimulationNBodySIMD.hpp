#ifndef SIMULATION_N_BODY_SIMD_HPP_
#define SIMULATION_N_BODY_SIMD_HPP_

#include <string>

#include <mipp.h>

#include "core/SimulationNBodyInterface.hpp"

class SimulationNBodySIMD : public SimulationNBodyInterface {
  protected:
    std::vector<accAoS_t<float>> accelerations; /*!< Array of body acceleration structures. */
    
    //MIPP modifications
    const int floatN = mipp::N<float>();
    accSoA_t<float> vAccelerations; /*!< Structure of arrays of body acceleration. */

  public:
    SimulationNBodySIMD(const unsigned long nBodies, const std::string &scheme = "galaxy", const float soft = 0.035f,
                         const unsigned long randInit = 0);
    virtual ~SimulationNBodySIMD() = default;
    virtual void computeOneIteration();

  protected:
    void initIteration();
    void computeBodiesAcceleration();
};

#endif /* SIMULATION_N_BODY_SIMD_HPP_ */
