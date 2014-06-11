#include "hartreefock.h"

#include <QDebug>
#include <cmath>

#include <electronsystems/gaussian/gaussiancore.h>
#include <electronsystems/gaussian/gaussiansystem.h>
#include "solvers/unrestrictedhartreefocksolver.h"
#include "solvers/restrictedhartreefocksolver.h"

using namespace std;

HartreeFock::HartreeFock(QQuickItem *parent) :
    QQuickItem(parent),
    m_nSampleSteps(0),
    m_voxelData(0),
    m_voxelEdgeMin(0),
    m_voxelEdgeMax(0),
    m_orbital(0)
{
    loadPointsFromFile();
}

HartreeFock::~HartreeFock()
{

}

void HartreeFock::generateRandomPoints() {
    m_filePositions.reset();
    m_filePositions = arma::zeros(1,20*20*20,3);
    int nPerDim = 20;
    int counter = 0;
    for(int i = 0; i < nPerDim; i++) {
        for(int j = 0; j < nPerDim; j++) {
            for(int k = 0; k < nPerDim; k++) {
                m_filePositions(0,counter,0) = -nPerDim / 2 + i;
                m_filePositions(0,counter,1) = -nPerDim / 2 + j;
                m_filePositions(0,counter,2) = -nPerDim / 2 + k;
                counter++;
            }
        }
    }
    m_nSampleSteps = 1;
    emit nSampleStepsChanged(m_nSampleSteps);
    emit dataChanged();
}

uint HartreeFock::voxelDataWidth()
{
    return m_orbitalDensities(0).n_rows;
}

uint HartreeFock::voxelDataHeight()
{
    return m_orbitalDensities(0).n_cols;
}

uint HartreeFock::voxelDataDepth()
{
    return m_orbitalDensities(0).n_slices;
}

const cube &HartreeFock::positions() const
{
    return m_filePositions;
}


void HartreeFock::loadPointsFromFile()
{
    double edge = 10;
    vec x = linspace(-edge, edge, 50);
    vec y = linspace(-edge, edge, 50);
    vec z = linspace(-edge, edge, 50);
    cout << "Solving system with Hartree-Fock..." << endl;
    vector<GaussianCore> cores;
//    cores.push_back(GaussianCore({0,           0.0,   0.0000}, "atom_8_basis_STO-3G.tm"));
    cores.push_back(GaussianCore({ -5,   0.0,    0.0000}, "atom_1_basis_6-311++Gdsds.tm"));
    cores.push_back(GaussianCore({  5,    0.0,    0.0000}, "atom_1_basis_6-311++Gdsds.tm"));
    GaussianSystem system;
    for(const GaussianCore &core : cores) {
        system.addCore(core);
    }
//    system.setNParticlesDown(9);
    mat C;
//    RestrictedHartreeFockSolver solver(&system);
    UnrestrictedHartreeFockSolver solver(&system);
    solver.setNIterationsMax(1e4);
    solver.setConvergenceTreshold(1e-10);
    solver.solve();
    cout << "Energy: " << solver.energy() << endl;
    C = join_rows(solver.coeffcientMatrixUp(), solver.coeffcientMatrixDown());
//    C = solver.coefficientMatrix();
    m_orbitalDensities.set_size(C.n_cols);
    m_totalDensity = zeros(x.n_elem, y.n_elem, z.n_elem);
    for(cube& orbitalDensity : m_orbitalDensities) {
        orbitalDensity = zeros(x.n_elem, y.n_elem, z.n_elem);
    }
    cout << "Calculating density..." << endl;
    double maxDensity = 0.0;
    double densitySum = 0.0;
    for(uint i = 0; i < x.n_elem; i++) {
        for(uint j = 0; j < y.n_elem; j++) {
            for(uint k = 0; k < z.n_elem; k++) {
                rowvec orbitalDensities = system.orbitalDensities(C, Vector3(x(i), y(j), z(k)));
                for(uint orbital = 0; orbital < orbitalDensities.n_elem; orbital++) {
                    double density = orbitalDensities(orbital);
                    densitySum += density;
                    m_totalDensity(i,j,k) += density;
                    m_orbitalDensities(orbital)(i,j,k) = density;
                    maxDensity = fmax(maxDensity, density);
                }
            }
        }
    }
    int elementCount = x.n_elem * y.n_elem * z.n_elem;
    double meanDensity = (densitySum / elementCount);
    m_totalDensity /= (meanDensity * m_orbitalDensities.n_elem);
    for(cube& orbitalDensity : m_orbitalDensities) {
        orbitalDensity /= meanDensity; // Normalize densities to the mean
    }
    m_energy = solver.energy();
    emit energyChanged(m_energy);
    double minValue = x.min();
    double maxValue = x.max();
    double mostMaxValue = max(-minValue, maxValue);
    minValue = -mostMaxValue;
    maxValue = mostMaxValue;
    setupVoxelData();
    m_voxelEdgeMin = minValue;
    m_voxelEdgeMax = maxValue;
    emit voxelEdgeMinChanged(m_voxelEdgeMin);
    emit voxelEdgeMaxChanged(m_voxelEdgeMax);
    emit nSampleStepsChanged(m_nSampleSteps);
    emit orbitalCountChanged(m_orbitalDensities.n_elem);
}

void HartreeFock::setupVoxelData() {
    if(m_orbitalDensities.n_elem < 1) {
        return;
    }
    if(m_voxelData) {
        delete[] m_voxelData;
    }
    uint nElements = m_orbitalDensities(0).n_rows * m_orbitalDensities(0).n_cols * m_orbitalDensities(0).n_slices;
    m_voxelData = new GLuint[nElements];
    cube *densityPointer = &m_totalDensity;
    if(m_orbital != -1) {
        densityPointer = &(m_orbitalDensities(m_orbital));
    }
    cube &density = *densityPointer;
    uint rowCount = density.n_rows;
    uint colCount = density.n_rows;
    uint sliceCount = density.n_rows;

    for(uint i = 0; i < (rowCount); i++) {
        for(uint j = 0; j < (colCount); j++) {
            for(uint k = 0; k < (sliceCount); k++) {
                int index = i
                        + j * rowCount
                        + k * colCount * rowCount;
                double value = density(i,j,k);
//                double value = pow(value / m_contrast, m_contrast);
                value = fmin(1.0, value);
                m_voxelData[index] = value * 2147483647;
            }
        }
    }
    emit dataChanged();
}

GLuint *HartreeFock::voxelData() const
{
    return m_voxelData;
}

void HartreeFock::setOrbital(int arg)
{
    if (m_orbital != arg) {
        m_orbital = arg;
        setupVoxelData();
        emit orbitalChanged(arg);
    }
}
