#include "hartreefock.h"

#include <QDebug>
#include <cmath>
#include <yaml-cpp/yaml.h>
#include <QString>
#include <iomanip>
#include <iostream>
#include <fstream>

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
    setupVoxelData();
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
    if(m_orbitalDensities.size() <= 0) {
        return 0;
    }
    return m_orbitalDensities(0).n_rows;
}

uint HartreeFock::voxelDataHeight()
{
    if(m_orbitalDensities.size() <= 0) {
        return 0;
    }
    return m_orbitalDensities(0).n_cols;
}

uint HartreeFock::voxelDataDepth()
{
    if(m_orbitalDensities.size() <= 0) {
        return 0;
    }
    return m_orbitalDensities(0).n_slices;
}

const cube &HartreeFock::positions() const
{
    return m_filePositions;
}


void HartreeFock::loadPointsFromFile()
{
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

void operator >> (const YAML::Node& node, Vector3& v)
{
    double x;
    double y;
    double z;
    node[0] >> x;
    node[1] >> y;
    node[2] >> z;
    v = Vector3(x,y,z);
}

void HartreeFock::openFile(QString fileName)
{
    qDebug() << "Opening file" << fileName;
    ifstream fin(fileName.toStdString());
    if(fin.fail()) {
        qWarning() << "Could not open" << fileName;
        return;
    }

    YAML::Parser parser(fin);

    YAML::Node rootNode;
    parser.GetNextDocument(rootNode);
    unsigned int nAtoms = rootNode["atoms"].size();
    string defaultBasis = "3-21G";

    GaussianSystem system;
    for(YAML::Iterator it=rootNode.begin();it!=rootNode.end();++it) {
        string rootKey;
        it.first() >> rootKey;
        if(rootKey == "atoms") {
            const YAML::Node &atomsNode = it.second();
            for(YAML::Iterator it2=atomsNode.begin();it2!=atomsNode.end();++it2) {
                const YAML::Node &atomNode = *it2;
                string typeAbbreviation;
                atomNode["type"] >> typeAbbreviation;
                Vector3 position;
                atomNode["position"] >> position;
                string basis;
                try {
                    atomNode["basis"] >> basis;
                } catch( YAML::TypedKeyNotFound<std::string> ) {
                    basis = defaultBasis;
                }

                basis = HF::escapeBasis(basis);
                system.addCore(GaussianCore(position, typeAbbreviation, basis));
            }
        }
    }

    double edge = 10;
    vec x = linspace(-edge, edge, 50);
    vec y = linspace(-edge, edge, 50);
    vec z = linspace(-edge, edge, 50);
    cout << "Solving system with Hartree-Fock..." << endl;
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
