/****************************************************************************
**
** Copyright (C) 2012 Digia Plc and/or its subsidiary(-ies).
** Contact: http://www.qt-project.org/legal
**
** This file is part of the Qt3D module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and Digia.  For licensing terms and
** conditions see http://qt.digia.com/licensing.  For further information
** use the contact form at http://qt.digia.com/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 2.1 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL included in the
** packaging of this file.  Please review the following information to
** ensure the GNU Lesser General Public License version 2.1 requirements
** will be met: http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.
**
** In addition, as a special exception, Digia gives you certain additional
** rights.  These rights are described in the Digia Qt LGPL Exception
** version 1.1, included in the file LGPL_EXCEPTION.txt in this package.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3.0 as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL included in the
** packaging of this file.  Please review the following information to
** ensure the GNU General Public License version 3.0 requirements will be
** met: http://www.gnu.org/copyleft/gpl.html.
**
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef VolumeSHADERPROGRAM_P_H
#define VolumeSHADERPROGRAM_P_H

#include <QtCore/qsharedpointer.h>
#include <QtCore/qpointer.h>
#include <QOpenGLShaderProgram>

#include "qquickeffect.h"
#include "qglshaderprogrameffect.h"

QT_BEGIN_NAMESPACE

class VolumeShaderProgram;
class VolumeShaderProgramEffect;

/*!
  \internal
    This class's purpose is to be a proxy for signals for the main shader
    program.

    This is necessary so the qt_metacall trick can be used to tell apart
    the update signals from QML generated properties.  (In short, each update
    signal is connected to an imaginary slot, but calling the slot is
    intercepted and replaced by a call to the appropriate update function).
*/
class VolumeShaderProgramPropertyListener : public QObject
{
    Q_OBJECT
public:
    VolumeShaderProgramPropertyListener(QObject* parent = 0)
        : QObject(parent)
    {
    }
    virtual ~VolumeShaderProgramPropertyListener()
    {
    }
Q_SIGNALS:
    void propertyChanged();
};

/*!
  \internal
  \sa VolumeShaderProgramPropertyListener
*/
class VolumeShaderProgramPropertyListenerEx : public VolumeShaderProgramPropertyListener
{
public:
    VolumeShaderProgramPropertyListenerEx(VolumeShaderProgram* parent, VolumeShaderProgramEffect* effect);
    ~VolumeShaderProgramPropertyListenerEx();

protected:
    virtual int qt_metacall(QMetaObject::Call c, int id, void **a);
    VolumeShaderProgramEffect* effect;
private:
    int VolumeshaderProgramMethodCount;
};


/*!
  \internal
  The VolumeShaderProgramEffect class underlies the VolumeShaderProgram class in Qml/3d.
  It contains the actual QGLVolumeShaderProgram along with all of the necessary
  parameters to use that program.
*/
class VolumeShaderProgramEffect : public QGLShaderProgramEffect
{
public:
    VolumeShaderProgramEffect(VolumeShaderProgram* parent);
    virtual ~VolumeShaderProgramEffect();

    bool create(const QString& vertexShader, const QString& fragmentShader);

    void update(QGLPainter *painter, QGLPainter::Updates updates);
    bool setUniformForPropertyIndex(int propertyIndex, QGLPainter *painter);

    void setPropertiesDirty();
    void setPropertyDirty(int property);

    void setAttributeFields(QGL::VertexAttribute fields);
protected:
    void processTextureUrl(int uniformLocation, QString urlString);
    void afterLink();
    bool beforeLink();

private:
    void setUniform(int uniformValue, const QImage& image,
                             QGLPainter* painter);
    void setUniform(int uniformValue, const QPixmap pixmap,
                             QGLPainter* painter);
    QGLTexture2D* textureForUniformValue(int uniformLocation);
    int textureUnitForUniformValue(int uniformLocation);

    QPointer<VolumeShaderProgram> parent;
    int nextTextureUnit;
    QMap<int, int> propertyIdsToUniformLocations;
    QMap<int, int> uniformLocationsToTextureUnits;
    QList<int> dirtyProperties;
    QArray<int> propertiesWithoutNotificationSignal;
    VolumeShaderProgramPropertyListener* propertyListener;

    // Thes maps are all referenced by uniform location
    QMap<int, QGLTexture2D*> texture2Ds;
    QMap<int, QImage> images;
    QMap<int, QString> urls;

    // These are sets of uniform locations
    QSet<int> loadingTextures;
    QSet<int> changedTextures;

    int m_texture3DuniformValue;
    int m_eyePositionUniformLocation;
};

QT_END_NAMESPACE

#endif // VolumeSHADERPROGRAM_P_H
