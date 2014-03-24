#ifndef POSITIONREADER_H
#define POSITIONREADER_H

#include <armadillo>
#include <QObject>
#include <QQuickItem>
#include <QVector3D>
#include <GL/gl.h>

using namespace arma;

class HartreeFock : public QQuickItem
{
    Q_OBJECT
    Q_PROPERTY(int nSampleSteps READ nSampleSteps NOTIFY nSampleStepsChanged)
    Q_PROPERTY(double voxelEdgeMin READ voxelEdgeMin WRITE setVoxelEdgeMin NOTIFY voxelEdgeMinChanged)
    Q_PROPERTY(double voxelEdgeMax READ voxelEdgeMax WRITE setVoxelEdgeMax NOTIFY voxelEdgeMaxChanged)
    Q_PROPERTY(double energy READ energy NOTIFY energyChanged)
    Q_PROPERTY(QVector3D center READ center WRITE setCenter NOTIFY centerChanged)
public:
    HartreeFock(QQuickItem* parent = 0);
    virtual ~HartreeFock();
    const cube& positions() const;
    void generateRandomPoints();
    int nSampleSteps() const
    {
        return m_nSampleSteps;
    }
    int voxelDataWidth();
    int voxelDataHeight();
    int voxelDataDepth();

    GLushort *voxelData() const;

    double voxelEdgeMax() const
    {
        return m_voxelEdgeMax;
    }

    double voxelEdgeMin() const
    {
        return m_voxelEdgeMin;
    }

    QVector3D center() const
    {
        return m_center;
    }

    double energy() const
    {
        return m_energy;
    }

public slots:
    void setVoxelEdgeMax(double arg)
    {
        if (m_voxelEdgeMax != arg) {
            m_voxelEdgeMax = arg;
            emit voxelEdgeMaxChanged(arg);
        }
    }

    void setVoxelEdgeMin(double arg)
    {
        if (m_voxelEdgeMin != arg) {
            m_voxelEdgeMin = arg;
            emit voxelEdgeMinChanged(arg);
        }
    }

    void setCenter(QVector3D arg)
    {
        if (m_center != arg) {
            m_center = arg;
            emit centerChanged(arg);
        }
    }

signals:
    void nSampleStepsChanged(int arg);
    void dataChanged();

    void voxelEdgeMaxChanged(double arg);

    void voxelEdgeMinChanged(double arg);

    void centerChanged(QVector3D arg);

    void energyChanged(double arg);

protected:
    void loadPointsFromFile();
    void setupVoxelData();
private:
    cube m_filePositions;
    cube m_densityVoxels;
    int m_nSampleSteps;
    GLushort *m_voxelData;
    double m_voxelEdgeMin;
    double m_voxelEdgeMax;
    QVector3D m_center;
    double m_energy;
};

#endif // POSITIONREADER_H
