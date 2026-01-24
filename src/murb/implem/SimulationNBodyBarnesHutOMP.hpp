#ifndef SIMULATION_N_BODY_BARNES_HUT_OMP_HPP_
#define SIMULATION_N_BODY_BARNES_HUT_OMP_HPP_

#include <string>

#include "core/SimulationNBodyInterface.hpp"
#include "Octotree.hpp"

class SimulationNBodyBarnesHutOMP : public SimulationNBodyInterface {
  protected:
    std::vector<accAoS_t<float>> accelerations; /*!< Array of body acceleration structures. */
    Octotree *root;
    const float theta;

  public:
    SimulationNBodyBarnesHutOMP(const unsigned long nBodies, const std::string &scheme = "galaxy", const float soft = 0.035f,
                         const unsigned long randInit = 0, const float theta = 1.5f);
    virtual ~SimulationNBodyBarnesHutOMP() = default;
    virtual void computeOneIteration();

  protected:
    void initIteration();
    void computeBodiesAcceleration();
    void updatePositionsAndVelocities();
};

#endif /* SIMULATION_N_BODY_BARNES_HUT_OMP_HPP_ */
