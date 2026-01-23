#ifndef SIMULATION_N_BODY_OPENMP_HPP_
#define SIMULATION_N_BODY_OPENMP_HPP_

#include <string>
#include <omp.h>

#include "core/SimulationNBodyInterface.hpp"

class SimulationNBodyOpenMP : public SimulationNBodyInterface {
  protected:
    std::vector<accAoS_t<float>> accelerations; /*!< Array of body acceleration structures. */

  public:
    SimulationNBodyOpenMP(const unsigned long nBodies, const std::string &scheme = "galaxy", const float soft = 0.035f,
                         const unsigned long randInit = 0);
    virtual ~SimulationNBodyOpenMP() = default;
    virtual void computeOneIteration();

  protected:
    void initIteration();
    void computeBodiesAcceleration();
    void updatePositionsAndVelocities();
};

#endif /* SIMULATION_N_BODY_NAIVE_HPP_ */
