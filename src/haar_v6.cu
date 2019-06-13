/*
*  haar_vX is the well synchronized version of haar transform into de GPU.
*  This version work with all the samples as a matrix into de GPU
*  with dimension SAMPLES x (sample_num + data_adjust) (rows x cols)
*  Copyright (C) 2018 Lisardo Fernández Cordeiro <lisardo.fernandez@uv.es>
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2, or (at your option)
*  any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
*
*/

/** \file
*  \brief Archivo para procesamiento de diferentes muestras metiladas de ADN.
*
*  Este archivo contiene la definición de las funciones para:
*         ..carga de datos en GPU
*         ..lanzamiento de proceso de transformación en GPU
*         ..kernel en GPU para control de transformación en niveles definidos
*         ..kernel de transformación wavelet del vector seleccionado
*         ..kernel para copiar coeficientes desde vector auxiliar de sincronización
*/

#include <stdio.h>
#include <GL/gl.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include "data_pack.h"

#define BLOCK_SIZE  1024		// número de hilos por bloque de GPU
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); } // para gestión de errores en GPU



/** ***********************************************************************************************
  * \fn void gpuAssert(cudaError_t, char*, int, bool)
  *  \brief Función responsable de recoger error en GPU y mostrarlo
  *  \param code	código de error de la GPU
  *  \param *file	fichero donde se produce el error
  *  \param line	línea de código donde se produce el error
  *  \param abort	indica si se sale del programa
  * ***********************************************************************************************
  */
extern "C"
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "GPUassert: %s %s %d\n\n", cudaGetErrorString(code), file, line);
        if (abort)
            exit(code);
    }
}

/** ***********************************************************************************************
  * \fn void copyValuesTotal(float*, float *, int, int)
  *  \brief función "hija" en GPU responsable de la copia de los datos del segmento a transformar
  *         proporcionando sincronización a nivel GRID
  *  \param *haar	puntero a vector de datos original
  *  \param *aux	puntero a vector de datos a transformar
  *  \param num		numero datos totales a transformar
  *  \param pi  	posicion inicial de copia - offset
  * ***********************************************************************************************
  */
extern "C"
__global__
void copyValuesTotal(float *haar, float *aux, int num, int posicion_inicial)
{
    // variables ------------------------------------------------------------------------------
    int index = threadIdx.x + blockIdx.x * blockDim.x;	// índice sobre todo el vector

    // copiar todos los valores de haar en aux
    if (index < num)
            aux[index] = haar[index + posicion_inicial];
}


/** ***********************************************************************************************
  * \fn void copyValues(float*, float *, int)
  *  \brief función "hija" en GPU responsable de la copia de los valores escalados
  *         proporcionando sincronización a nivel GRID
  *  \param *aux	puntero a vector de datos a transformar
  *  \param *temp	puntero a vector de datos temporales a copiar
  *  \param num     numero datos a copiar
  * ***********************************************************************************************
  */
extern "C"
__global__
void copyValues(float *aux, float *temp, int num)
{
    // variables ----------------------------------------------------------------------------------
    int index = threadIdx.x + blockIdx.x * blockDim.x;	// índice sobre todo el vector

    // copiar todos los valores de haar en aux
    if (index < num)
        aux[index] = temp[index];
}


/** ***********************************************************************************************
  * \fn void transform(float*, int, int, int)
  *  \brief función "hija" en GPU responsable de la transformación wavelet de un vector
  *         proporcionando sincronización a nivel GRID.
  *  \param *aux	puntero a vector de datos a transformar
  *  \param *temp	puntero a vector de resultados intermedios
  *  \param num		número de posiciones del vector
  * ***********************************************************************************************
  */
extern "C"
__global__
void transform(float *aux, float *temp, int num)
{
    // variables ----------------------------------------------------------------------------------
    int index = threadIdx.x + blockIdx.x * blockDim.x;	// índice sobre todo el vector
    float f   = 0.7071067811865476;                     // coeficiente haar wavelet
    float aux1;                                         // variables auxiliares de sincronización
    int idx;                                            // indice auxiliar para guardar dato

    // transformada haar en paralelo sobre el vector recibido -------------------------------------
    if (index < num)
    {
        if ((index & 0x01) == 0)	// solo los hilos con índice par (0, 2, 4, ...)
        {
            idx = index * 0.5;

            aux1 = (aux[index] + aux[index + 1]) * f;	// escalado (filtro paso-bajo)

            temp[idx]  = aux1;
        }
    }
}


/** ***********************************************************************************************
  * \fn void wavedec(float*, float**, int, int, int, int, int)
  *  \brief Función principal en GPU responsable de calcular y ordenar las partes del vector
  *         para su transformación wavelet multinivel.
  *  \param *haar	puntero a matriz de datos a transformar
  *  \param *aux	puntero a matriz de coeficiente auxiliares para ayuda a la sincronización
  *  \param pitch	desplazamiento óptimo en memoria GPU para alojar cada muestra
  *  \param pitch_2	desplazamiento óptimo en memoria GPU para alojar cálculo auxiliar
  *  \param n		número total de posiciones del vector
  *  \param l		número de niveles a computar
  *  \param samples número de muestras a analizar
  *  \param pi      posición inicial del segmento de datos a analizar
  * ***********************************************************************************************
  */
extern "C"
__global__
void wavedec(float *haar, float *aux, float *temp,
             size_t pitch, size_t pitch_2, size_t pitch_3,
             int n, int l, int samples, int pi)
{
    // variables ----------------------------------------------------------------------------------
    int index_X = threadIdx.x + blockIdx.x * blockDim.x;	// indice de hilos sobre todo el vector
    int level   = 0;                                        // número de nivel
    int num     = n;                                        // número de posiciones en vector
    int hilo;                                               // guarda el hilo asignado para que se resposabilice de todo el proceso

    // limita el número de hilos al de muestras ---------------------------------------------------
    if (index_X < samples)
    {
        hilo = index_X;		// cada hilo se responsabiliza de una misma muestra

        if (hilo == index_X)
        {
            // separar los datos por muestras - - - - - - - - - - - - - - - - - - - - - - - - - - -
            float *haar_c = (float *)((char *)haar + index_X * pitch);
            float *aux_c  = (float *)((char *)aux  + index_X * pitch_2);
            float *temp_c = (float *)((char *)temp + index_X * pitch_3);

            __syncthreads();

            // llamada a función hija para copiar segmento de vector a transformar
            copyValuesTotal<<<(num + BLOCK_SIZE-1) / BLOCK_SIZE, BLOCK_SIZE>>>(haar_c,
                                                                               aux_c,
                                                                               num,
                                                                               pi);


            // procesamiento multinivel del vector de datos ---------------------------------------
            // repite la transformación tantas veces como niveles se han solicitado
            while (level < l && num >= 2)
            {
                // llamada a función hija para transformación del nivel correspondiente
                // con esta división en padre-hijo, se consigue sincronizar cada nivel
                // \param	<<<((datos_x_muestra + num_hilos_bloque-1) / num_hilos_bloque),
                // 		numero hilos por bloque>>>
                transform<<<(num + BLOCK_SIZE-1) / BLOCK_SIZE, BLOCK_SIZE>>>(aux_c,
                                                                             temp_c,
                                                                             num);


                // actualizar variables de nivel  - - - - - - - - - - - - - - - - - - - - - - - - -
                level += 1;
                num    = ceilf(num * 0.5);


                // llamada a función hija para copiar resultados en vector auxiliar
                copyValues<<<(num + BLOCK_SIZE-1) / BLOCK_SIZE, BLOCK_SIZE>>>(aux_c,
                                                                              temp_c,
                                                                              num);


                // actualiza el número de datos para el siguiente nivel - - - - - - - - - - - - - -
                if ((num & 01) == 1)
                {
                    num++;
                    aux_c[num] = 0;
                }
            }
        }
    }
}

/** ***********************************************************************************************
  * \fn void cuda_send_data(datos_cuda &)
  *  \brief Función para enviar los datos a la GPU
  *  \param &cuda_data  estructura con variables de control de datos
  * ***********************************************************************************************
  */
void cuda_send_data(datos_cuda &cuda_data)
{
    // reserva espacio en GPU para el vector a transformar y copia matriz de datos ----------------
    // devuelve valor de desplazamiento (pitch) óptimo para gestión de memoria adecuada
    // en función de la cantdad de datos a alojar
    // \param 	puntero a posición memoria GPU,
    //          desplazamiento óptimo devuelto por CUDA,
    //          cantidad de bytes a reservar por fila,
    //          número de muestras (filas)
    gpuErrchk(cudaMallocPitch(&cuda_data.d_haar,
                              &cuda_data.pitch,
                              (cuda_data.sample_num + cuda_data.data_adjust) * sizeof(float),
                              cuda_data.samples));


    gpuErrchk(cudaMallocPitch(&cuda_data.d_aux,
                              &cuda_data.pitch_2,
                              (cuda_data.sample_num + cuda_data.data_adjust) * sizeof(float),
                              cuda_data.samples));


    // envío de datos a GPU -----------------------------------------------------------------------
    // \param	puntero a posición de memoria GPU,
    //          desplazamiento óptimo,
    //          puntero a posición de datos en CPU a enviar a GPU,
    //          cantidad de bytes a enviar por muestra,
    //          cantidad de bytes a alojar por muestra,
    //          número de filas (muestras)
    gpuErrchk(cudaMemcpy2D( cuda_data.d_haar,
                            cuda_data.pitch,
                            cuda_data.mc_full[0],
                            cuda_data.sample_num * sizeof(float),
                            cuda_data.sample_num * sizeof(float),
                            cuda_data.samples,
                            cudaMemcpyHostToDevice));

}


/** ***********************************************************************************************
  * \fn void cuda_main(datos_cuda &)
  *  \brief Función para procesar los datos en la GPU
  *  \param &cuda_data  estructura con variables de control de datos
  * ***********************************************************************************************
  */
void cuda_main(datos_cuda &cuda_data)
{
    // reserva TODA la memoria CONTIGUA para la matriz de muestras tranformadas -------------------
    // para trasvase de datos entre GPU y CPU con CUDA, la matriz debe ser contigua completa
    cuda_data.h_haar_C = new float*[cuda_data.samples];                             // reservar punteros a filas
    cuda_data.h_haar_C[0] = new float[cuda_data.samples * (cuda_data.h_haar_L[0])];	// reservar toodos los datos (rows * cols)
    for (int i = 1; i < cuda_data.samples; i++)                                     // asignar valor a punteros de fila
        cuda_data.h_haar_C[i] = cuda_data.h_haar_C[i - 1] + cuda_data.h_haar_L[0];


    // reserva memoria para cálculos temporales en GPU --------------------------------------------
    float *d_temp;
    size_t pitch;
    gpuErrchk(cudaMallocPitch(&d_temp,
                              &pitch,
                              (cuda_data.sample_num + 1) * 0.7 * sizeof(float),
                              cuda_data.samples));


    // transforma el número de muestras elegida ---------------------------------------------------
    // realiza la transformación en la GPU del conjunto de muestras cargado
    // \param	<<< número de bloques a utilizar,
    //          número de hilos por bloque >>> (máximo 1024 para PASCAL GTX 1080)
    // \param	puntero a datos a transformar alojados en GPU,
    //          desplazamiento óptimo de datos por fila,
    //          número de datos por muestra (fila) a transformar,
    //          ajuste de longitud de muestra por número impar al dividir la muestra
    wavedec<<<1, cuda_data.samples>>>(cuda_data.d_haar,
                                      cuda_data.d_aux,
                                      d_temp,
                                      cuda_data.pitch,
                                      cuda_data.pitch_2,
                                      pitch,
                                      cuda_data.sample_num,
                                      cuda_data.levels,
                                      cuda_data.samples,
                                      cuda_data.rango_inferior);

    // espera a que la GPU termine el trabajo - - - - - - - - - - - - - - - - - - - - - - - - - - -
    gpuErrchk(cudaDeviceSynchronize());


    // recupera el resultado de la transformación en memoria GPU a memoria CPU- - - - - - - - - - -
    // \param	puntero a matriz de datos a guardar en CPU,
    //          cantidad de bytes a guardar por muestra,
    //          puntero a datos para copiar de GPU,
    //          desplazamiento óptimo de datos por fila en GPU,
    //          cantidad de bytes en GPU a copiar por muestra,
    //          número de muestras (filas)
    gpuErrchk(cudaMemcpy2D(	cuda_data.h_haar_C[0],
                            cuda_data.h_haar_L[0] * sizeof(float),
                            cuda_data.d_aux,
                            cuda_data.pitch,
                            cuda_data.h_haar_L[0] * sizeof(float),
                            cuda_data.samples,
                            cudaMemcpyDeviceToHost));


    //libera la memoria temporal utilizada para cálculos intemedios
    cudaFree(d_temp);
}

/** ***********************************************************************************************
  * \fn void *cuda_init()
  *  \brief Función para inicializar la gpu
  * ***********************************************************************************************
  */
void cuda_init()
{
    int deviceCount = 0;
    int cudaDevice  = 0;
    char cudaDeviceName [100];
    cudaDeviceProp prop;
    cuInit(0);
    cuDeviceGetCount(&deviceCount);
    cuDeviceGet(&cudaDevice, 0);
    cuDeviceGetName(cudaDeviceName, 100, cudaDevice);
    cudaGetDeviceProperties(&prop, cudaDevice);

    if (cudaChooseDevice(&cudaDevice, &prop) != cudaSuccess)
        puts("failed to choose device");
    if (cudaGLSetGLDevice(cudaDevice) != cudaSuccess)
        puts("failed to set gl device");

    printf("Number of devices: %u \t cuda device: %d\n", deviceCount, cudaDevice);
    printf("Device name: %s\n", cudaDeviceName);
    printf("Warp size: %u\n", prop.warpSize);
}

/** ***********************************************************************************************
  * \fn void cuda_end(data buf)
  *  \brief Función para liberar memoria de la GPU
  *  \param &cuda_data  estructura con variables de control de datos
  * ***********************************************************************************************
  */
void cuda_end(datos_cuda &cuda_data)
{
    //libera la memoria de la gpu utilizada para cálculos intemedios
    cudaFree(cuda_data.d_haar);
    cudaFree(cuda_data.d_aux);
}

/** ***********************************************************************************************
  * \fn void calculo_haar_L(datos_cuda &cuda_data)
  *  \brief Función para calcular el número de datos en el nivel dado y el ajuste por impares
  *  \param &cuda_data  estructura con variables de control de datos
  * ***********************************************************************************************
  */
void cuda_calculo_haar_L(datos_cuda &cuda_data)
{
    // cálculo de número de coeficientes por nivel y del ajuste de paso entre escala y coeficiente
    cuda_data.h_haar_L.push_front(cuda_data.sample_num);	// última posición guarda el total de posiciones por muestra

    // para cada nivel se divide por dos la cantidad de posiciones del nivel anterior -------------
    // redondeando al alza y actualizando el ajuste cuando sea impar
    for (int fila = cuda_data.levels; fila > 0; fila--)
    {
        if (ceil(cuda_data.h_haar_L.front() * 0.5 >= 2))
        {
            cuda_data.h_haar_L.push_front(ceil(cuda_data.h_haar_L.front() * 0.5));
            if (fila > 0 && cuda_data.h_haar_L[1] != cuda_data.sample_num)
                cuda_data.data_adjust += size_t(2 * cuda_data.h_haar_L.front() - cuda_data.h_haar_L[1]);
        }
        else
            break;
    }
    cuda_data.h_haar_L.push_front(cuda_data.h_haar_L.front());	// primera posición coincide con el número de datos de escala
}
