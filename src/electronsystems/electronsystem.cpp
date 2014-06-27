#include "electronsystem.h"

/*!
 * \class ElectronSystem
 * \brief Base class for any type of electronic system
 */

ElectronSystem::ElectronSystem() :
    m_nParticlesDownSet(false),
    m_nParticlesDown(0)
{
}

uint ElectronSystem::nParticlesUp() {
    return nParticles() - nParticlesDown();
}

uint ElectronSystem::nParticlesDown() {
    if(m_nParticlesDownSet) {
        if(m_nParticlesDown > nParticles()) {
            std::cerr << "Error: More particles set to down than the number of particles in the system." << std::endl;
            throw std::logic_error("");
        }
        return m_nParticlesDown;
    } else {
        return nParticles() / 2;
    }
}

void ElectronSystem::setNParticlesDown(uint nParticlesDown)
{
    if(nParticlesDown > nParticles()) {
        std::cout << "WARNING: The number of particles set to down is currently more than the number of particles in the system." << std::endl;
    }
    m_nParticlesDown = nParticlesDown;
    m_nParticlesDownSet = true;
}
