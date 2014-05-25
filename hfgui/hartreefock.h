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
    Q_PROPERTY(int orbital READ orbital WRITE setOrbital NOTIFY orbitalChanged)
    Q_PROPERTY(int orbitalCount READ orbitalCount NOTIFY orbitalCountChanged)
    Q_PROPERTY(double contrast READ contrast WRITE setContrast NOTIFY contrastChanged)
public:
    HartreeFock(QQuickItem* parent = 0);
    virtual ~HartreeFock();
    const cube& positions() const;
    void generateRandomPoints();
    int nSampleSteps() const
    {
        return m_nSampleSteps;
    }
    uint voxelDataWidth();
    uint voxelDataHeight();
    uint voxelDataDepth();

    GLuint *voxelData() const;

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

    int orbital() const
    {
        return m_orbital;
    }

    int orbitalCount() const
    {
        return m_orbitalDensities.n_elem;
    }

    double contrast() const
    {
        return m_contrast;
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

    void setOrbital(int arg);

    void setContrast(double arg)
    {
        if (m_contrast != arg) {
            m_contrast = arg;
            if(arg != 0) {
                m_contrastInverse = 1.0 / arg;
            } else {
                m_contrastInverse = 1e9;
            }
//            setupVoxelData();
            emit contrastChanged(arg);
        }
    }

signals:
    void nSampleStepsChanged(int arg);
    void dataChanged();

    void voxelEdgeMaxChanged(double arg);

    void voxelEdgeMinChanged(double arg);

    void centerChanged(QVector3D arg);

    void energyChanged(double arg);

    void orbitalChanged(int arg);

    void orbitalCountChanged(int arg);

    void contrastChanged(double arg);

protected:
    void loadPointsFromFile();
    void setupVoxelData();
private:
    cube m_filePositions;
    field<cube> m_orbitalDensities;
    cube m_totalDensity;
    int m_nSampleSteps;
    GLuint *m_voxelData;
    double m_voxelEdgeMin;
    double m_voxelEdgeMax;
    QVector3D m_center;
    double m_energy;
    int m_orbital;
    double m_contrast;
    double m_contrastInverse;
};

#endif // POSITIONREADER_H
