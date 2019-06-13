#ifndef OGL_GRAPHIC_H
#define OGL_GRAPHIC_H

#include <QSize>
#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include "data_pack.h"


using namespace std;

/** ***********************************************************************************************
  *  \brief declaración de funciones externas para compilación con nvcc
  *  \fn    void *cuda_registerBuffer(GLuint buf)
  *  \fn    void cuda_unregisterBuffer(void *res)
  *  \fn    void *cuda_map(void *res)
  *  \fn    void cuda_unmap(void *res)
  *  \fn    void cuda_init()
  * ***********************************************************************************************
  */
//extern
void *cuda_registerBuffer(GLuint buf);
//extern
void  cuda_unregisterBuffer(void *res);
//extern
void *cuda_map(void *res);
//extern
void  cuda_unmap(void *res);
//extern
void  cuda_init();



class OGL_graphic : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    OGL_graphic(QWidget *parent = 0);
    ~OGL_graphic();

    /** ***********************************************************************************************
      * \fn void void initializeGL()
      *  \brief Función responsable de inicializar las variables del openGL y la comunicación con CUDA
      * ***********************************************************************************************
      */
    void initializeGL();

    /** ***********************************************************************************************
      * \fn void paintGL()
      *  \brief Función responsable de dibujar cada update() con los datos residentes en la tarjeta
      * ***********************************************************************************************
      */
    void paintGL();

    /** ***********************************************************************************************
      * \fn void setNum_X(int value)
      *  \brief Función responsable de actualizar los datos de longitud de datos y de muestras
      * ***********************************************************************************************
      */
    void setNum_L(int value);
    void setNum_samples(int value);

    /** ***********************************************************************************************
      * \fn void registerBuffer()
      *  \brief Función responsable de registrar un manejador CUDA con el VBO
      * ***********************************************************************************************
      */
    void registerBuffer();

    /** ***********************************************************************************************
      * \fn void mapResource(datos_cuda &cuda_data)
      *  \brief Función responsable de mapear un puntero a los datos a dibujar con el manejados CUDA
      * ***********************************************************************************************
      */
    void mapResource(datos_cuda &cuda_data);

    /** ***********************************************************************************************
      * \fn void unmapResource()
      *  \brief Función responsable de abrir explorador de archivos y capturar fichero a analizara
      * ***********************************************************************************************
      */
    void unmapResource();

    /** ***********************************************************************************************
      * \fn void unregisterBuffer()
      *  \brief Función responsable de deshacer el mapeado para liberar al manejador y poder dibujar
      * ***********************************************************************************************
      */
    void unregisterBuffer();


signals:
    /** ***********************************************************************************************
      * \fn void ogl_coordinates(int x, int y, int xR, int yR)
      *  \brief señal de atención a pulsación y suelta de botón derecho de ratón sobre gráfica
      *  \param x, y    coordenadas de press
      *  \param xR, yR  coordenadas de release
      * ***********************************************************************************************
      */
    void ogl_coordinates(int x);
    void ogl_coordinates(int x, int y, int xR, int yR);


private:

    /** ***********************************************************************************************
      *  \brief variable vinculada a un VBO ARRAY para dibujar la gráfica
      *  \param m_buf GLuint para vincular a un VBO
      * ***********************************************************************************************
      */
    GLuint m_buf;

    /** ***********************************************************************************************
      *  \brief puntero a manejador de CUDA para establecer el enlace entre el VBO y el puntero CUDA
      *  \param *m_cudaBufHndle puntero manejador CUDA
      * ***********************************************************************************************
      */
    void  *m_cudaBufHandle;

    /** ***********************************************************************************************
      *  \brief variables enteras para almacenar características de datos a dibujar
      *  \param num_L       número de datos por cada muestra
      *  \param num_samples número de muestras a dibujar
      * ***********************************************************************************************
      */
    int    num_L;
    int    num_samples;

    /** ***********************************************************************************************
      *  \brief puntero a la posición de memoria que maneja CUDA y devuelve al manejador.
      *  \param *devPtr puntero a posición de memoria de la GPU donde residen los datos a dibujar
      * ***********************************************************************************************
      */
    void  *devPtr;

    /** ***********************************************************************************************
      * \fn void mouseXXXXeEvent(QMouseEvent *)
      *  \brief puntero a eventos de ratón sobre la gráfica para recoger la posición
      *         al click para abrir página web con zona de cromosoma, release para hacer zoom y
      *         movimiento para conocer la posición en el genoma que se está viendo
      *  \param *event  posición del ratón en la ventana de la gráfica
      * ***********************************************************************************************
      */
    void mousePressEvent(QMouseEvent * event);
    void mouseReleaseEvent(QMouseEvent * event);
    void mouseMoveEvent(QMouseEvent * event);

    /** ***********************************************************************************************
      *  \brief posición del ratón al click y al release
      *  \param xAtPress    posición al click
      *  \param xAtRelease  posición al release
      * ***********************************************************************************************
      */
    int xAtPress, yAtPress;
    int xAtRelease, yAtRelease;



};

#endif // OGL_GRAPHIC_H
