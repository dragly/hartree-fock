#include "electronsystem.h"

ElectronSystem::ElectronSystem() :
    m_nParticlesDownSet(false),
    m_nParticlesDown(0)
{
}

void ElectronSystem::setNParticlesDown(int nParticles)
{
    m_nParticlesDown = nParticles;
    m_nParticlesDownSet = true;
}
