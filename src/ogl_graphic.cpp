#include "ogl_graphic.h"
#include <QDebug>
#include <QMouseEvent>

OGL_graphic::OGL_graphic(QWidget *parent)
    : QOpenGLWidget(parent)
{
    num_L       = 1000;
    num_samples = 6;

    // activa el seguimiento del ratón por la ventana de la gráfica
    // generando eventos de movimiento sin nencesidad de hacer click
    setMouseTracking(true);
}

OGL_graphic::~OGL_graphic()
{

}

// ************************************************************************************************
void OGL_graphic::initializeGL()
{
    initializeOpenGLFunctions();

    cuda_init();

    glClearColor(1,1,1,1);      // fondo de ventana de gráfica en blanco
//    glClearColor(0,0,0,1);      // fondo de ventana de gráfica en negro
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    glGenBuffers(1, &m_buf);
    glBindBuffer(GL_ARRAY_BUFFER, m_buf);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * uint(num_L), NULL, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    m_cudaBufHandle = cuda_registerBuffer(m_buf);

}

// ************************************************************************************************
void OGL_graphic::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBindBuffer(GL_ARRAY_BUFFER, m_buf);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(2, GL_FLOAT, 0, 0); // num axis, type data, stride, start data
    // asigna colores diferentes a cada muestra
    for (int i = 0; i < num_samples; i++)
    {
        if (i * 2 >= num_samples)
            glColor3f(1.0, float(0.25 * (-(num_samples - 2*i))), 0.0);
        else
            glColor3f(0.0, float(0.4 * i), 1.0);
        glDrawArrays(GL_LINE_STRIP, i * 2 * num_L, 2 * num_L);  // type line, start point, num points
    }
    glDisableClientState(GL_VERTEX_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// ************************************************************************************************
void OGL_graphic::registerBuffer()
{
    cuda_unregisterBuffer(m_cudaBufHandle);
    glBindBuffer(GL_ARRAY_BUFFER, m_buf);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 4 * uint(num_samples * num_L), NULL, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    m_cudaBufHandle = cuda_registerBuffer(m_buf);

//    qDebug() <<"num_L: " << num_L;
}

// ************************************************************************************************
void OGL_graphic::mapResource(datos_cuda &cuda_data)
{
    cuda_data.d_glPtr = cuda_map(m_cudaBufHandle);
}

// ************************************************************************************************
void OGL_graphic::unmapResource()
{
    cuda_unmap(m_cudaBufHandle);
}

// ************************************************************************************************
void OGL_graphic::unregisterBuffer()
{
    cuda_unregisterBuffer(m_cudaBufHandle);
    glBindBuffer(GL_ARRAY_BUFFER, m_buf);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float), NULL, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// ************************************************************************************************
void OGL_graphic::setNum_samples(int value)
{
    num_samples = value;
}

// ************************************************************************************************
void OGL_graphic::setNum_L(int value)
{
    num_L = value;
}

// ************************************************************************************************
void OGL_graphic::mousePressEvent(QMouseEvent * event)
{
    xAtPress = event->x();
    yAtPress = event->y();
}

// ************************************************************************************************
void OGL_graphic::mouseReleaseEvent(QMouseEvent * event)
{
    xAtRelease = event->x();
    yAtRelease = event->y();

    emit ogl_coordinates(xAtPress, yAtPress, xAtRelease, yAtRelease);
}

// ************************************************************************************************
void OGL_graphic::mouseMoveEvent(QMouseEvent *event)
{
    emit ogl_coordinates(event->x());
}

