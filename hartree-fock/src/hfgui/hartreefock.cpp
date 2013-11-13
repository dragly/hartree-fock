#include "hartreefock.h"

#include <QDebug>
#include <cmath>

#include <electronsystems/gaussian/gaussiancore.h>
#include <electronsystems/gaussian/gaussiansystem.h>
#include <hartreefocksolver.h>

using namespace std;

HartreeFock::HartreeFock(QQuickItem *parent) :
    QQuickItem(parent),
    m_nSampleSteps(0),
    m_voxelData(0),
    m_voxelEdgeMin(0),
    m_voxelEdgeMax(0)
{
    loadPointsFromFile();
}

HartreeFock::~HartreeFock()
{

}

void HartreeFock::generateRandomPoints() {
    m_filePositions.reset();
    m_filePositions = arma::zeros(1,20*20*20,3);
    double spacing = 1;
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

int HartreeFock::voxelDataWidth()
{
    return m_densityVoxels.n_rows;
}

int HartreeFock::voxelDataHeight()
{
    return m_densityVoxels.n_cols;
}

int HartreeFock::voxelDataDepth()
{
    return m_densityVoxels.n_slices;
}

const cube &HartreeFock::positions() const
{
    return m_filePositions;
}


void HartreeFock::loadPointsFromFile()
{
    cout << "Setting up Hartree Fock!" << endl;
    vector<GaussianCore> cores;
    cores.push_back(GaussianCore({0,0,0}, "oxygen431g.tm"));
    cores.push_back(GaussianCore({-1.43,1.108,0}, "hydrogen431g.tm"));
    cores.push_back(GaussianCore({1.43,1.108,0}, "hydrogen431g.tm"));
    cores.push_back(GaussianCore({0.0,-1.81,0}, "hydrogen431g.tm"));
    GaussianSystem system;
    for(const GaussianCore &core : cores) {
        system.addCore(core);
    }
    mat C;
    HartreeFockSolver solver(&system);
    for(int i = 0; i < 100; i++) {
        solver.advance();
    }
    cout << "Energy: " << solver.energy() << endl;
    C = solver.coefficientMatrix();
    vec x = linspace(-3, 3, 50);
    vec y = linspace(-3, 3, 50);
    vec z = linspace(-3, 3, 50);
    double dx = x(1) - x(0);
    double dy = y(1) - y(0);
    double dz = z(1) - z(0);
    double densitySum = 0;
    m_densityVoxels = cube(x.n_elem, y.n_elem, z.n_elem);
    for(uint i = 0; i < x.n_elem; i++) {
        cout << "Calculating density for x = " << x(i) << endl;
        for(uint j = 0; j < y.n_elem; j++) {
            for(uint k = 0; k < z.n_elem; k++) {
                double density = system.particleDensity(C, x(i), y(j), z(k));
                densitySum += density * dx * dy * dz;
                m_densityVoxels(i,j,k) = density;
            }
        }
    }
    cout << "Density sum: " << densitySum << endl;
    m_densityVoxels -= m_densityVoxels.min();
    m_densityVoxels += 1e-6;
    m_densityVoxels = sqrt(m_densityVoxels);

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
    emit dataChanged();
}

void HartreeFock::setupVoxelData() {
    uint nElements = m_densityVoxels.n_rows * m_densityVoxels.n_cols * m_densityVoxels.n_slices;
    if(m_voxelData) {
        delete[] m_voxelData;
    }
    m_voxelData = new GLushort[nElements];
    for(int i = 0; i < m_densityVoxels.n_rows; i++) {
        for(int j = 0; j < m_densityVoxels.n_cols; j++) {
            for(int k = 0; k < m_densityVoxels.n_slices; k++) {
                int index = i
                        + j * m_densityVoxels.n_rows
                        + k * m_densityVoxels.n_cols * m_densityVoxels.n_rows;
                m_voxelData[index] = m_densityVoxels(i,j,k) * 65534;
            }
        }
    }
}

GLushort *HartreeFock::voxelData() const
{
    return m_voxelData;
}
