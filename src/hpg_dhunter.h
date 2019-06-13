/*
*  HPG_Dhunter is the main class to control and manipulate data files
*  Copyright (C) 2018 Lisardo Fernández Cordeiro <lisardo.fernandez@uv.es>
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 3, or (at your option)
*  any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
*  or see <https://www.gnu.org/licenses/>.
*
*/

/** \file
*  \brief Programa para procesamiento y visualización de diferentes
*         muestras metiladas de ADN.
*
*  Este archivo contiene la definición de las funciones para:
*         ..declaración de funciones externas de procesamiento en GPU
*         ..control de rango de datos a analizar y visualizar
*         ..selección de fichero a analizaar
*/

#ifndef HPG_DHUNTER_H
#define HPG_DHUNTER_H

#include <QMainWindow>
#include <QTextCursor>
#include <QThread>
#include <QMutex>
#include <cuda_runtime.h>
#include <cuda.h>
#include <chrono>
#include "data_pack.h"
#include "files_worker.h"
#include "refgen.h"

#define TIMING

#ifdef TIMING
#define INIT_TIMER_1        auto start_1 = std::chrono::high_resolution_clock::now();
#define INIT_TIMER_2        auto start_2 = std::chrono::high_resolution_clock::now();
#define INIT_TIMER_3        auto start_3 = std::chrono::high_resolution_clock::now();
#define START_TIMER_1       start_1 = std::chrono::high_resolution_clock::now();
#define START_TIMER_2       start_2 = std::chrono::high_resolution_clock::now();
#define START_TIMER_3       start_3 = std::chrono::high_resolution_clock::now();
#define STOP_TIMER_1(name)  qDebug() << "DURACION de " << name << ": " << \
                          std::chrono::duration_cast<std::chrono::milliseconds>( \
                          std::chrono::high_resolution_clock::now()-start_1 \
                          ).count() << " ms ";
#define STOP_TIMER_2(name)  qDebug() << "DURACION de " << name << ": " << \
                          std::chrono::duration_cast<std::chrono::milliseconds>( \
                          std::chrono::high_resolution_clock::now()-start_2 \
                          ).count() << " ms ";
#define STOP_TIMER_3(name)  qDebug() << "DURACION de " << name << ": " << \
                          std::chrono::duration_cast<std::chrono::milliseconds>( \
                          std::chrono::high_resolution_clock::now()-start_3 \
                          ).count() << " ms ";

#else
#define INIT_TIMER_1
#define START_TIMER_1
#define STOP_TIMER_1(name)
#define INIT_TIMER_2
#define START_TIMER_2
#define STOP_TIMER_2(name)
#define INIT_TIMER_3
#define START_TIMER_3
#define STOP_TIMER_3(name)
#endif


using namespace std;


/** ***********************************************************************************************
  *  \brief declaración de funciones externas para compilación con nvcc
  *  \fn    void cuda_send_data(datos_cuda &)
  *  \fn    void cuda_main(datos_cuda &)
  * ***********************************************************************************************
  */
//extern
void cuda_send_data(datos_cuda &);
//extern
void cuda_main(datos_cuda &);
//extern
void cuda_calculo_haar_L(datos_cuda &);
//extern
void cuda_init();
//extern
void cuda_end(datos_cuda &);

namespace Ui {
class HPG_Dhunter;
}

class HPG_Dhunter : public QMainWindow
{
    Q_OBJECT

public:
    explicit HPG_Dhunter(QWidget *parent = nullptr);
    ~HPG_Dhunter();

    /** ***********************************************************************************************
      *  \brief variable de tipo estructura con variables para control de datos a analizar
      *  \param cuda_data   estrutura con variables de control y datos
      * ***********************************************************************************************
      */
    datos_cuda cuda_data;

private slots:
    /** ***********************************************************************************************
      * \fn void on_case_file_clicked() and four more
      *  \brief Funciones responsables de abrir explorador de archivos y capturar directorio a analizar
      *         resto de funciones para ordenamiento y borrado de directorios
      * ***********************************************************************************************
      */
    void on_case_file_clicked();
    void on_case_files_cursorPositionChanged();
    void on_delete_case_clicked();
    void on_up_case_clicked();
    void on_down_case_clicked();

    /** ***********************************************************************************************
      * \fn void on_control_file_clicked() and four more
      *  \brief Funciones responsables de abrir explorador de archivos y capturar directorio a analizar
      *         resto de funciones para ordenamiento y borrado de directorios
      * ***********************************************************************************************
      */
    void on_control_file_clicked();
    void on_control_files_cursorPositionChanged();
    void on_delete_control_clicked();
    void on_up_control_clicked();
    void on_down_control_clicked();

    /** ***********************************************************************************************
      * \fn void on_out_path_clicked()
      *  \brief Función responsable de seleccionar directorio donde guardar los resultados
      * ***********************************************************************************************
      */
    void on_out_path_clicked();

    /** ***********************************************************************************************
      * \fn void on_mC_clicked() and four more
      *  \brief Funciones responsables de adquirir los parámetros para el análisis
      * ***********************************************************************************************
      */
    void on_mC_clicked();
    void on_hmC_clicked();
    void on_forward_clicked();
    void on_reverse_clicked();
    void on_all_chroms_toggled(bool);

    /** ***********************************************************************************************
      * \fn void on_threshold_valueChanged(int value) and five more
      *  \brief Funciones responsables de adquirir los ajustes sobre los parámetros para el análisis
      *  \param value   valor del parámetro
      * ***********************************************************************************************
      */
    void on_threshold_valueChanged(int value);
    void on_mC_cobertura_sliderMoved(int value);
    void on_hmC_cobertura_sliderMoved(int value);
    void on_dmr_dwt_level_valueChanged(int value);
    void on_mC_min_cov_textEdited(const QString &arg1);
    void on_hmC_min_cov_textEdited(const QString &arg1);

    /** ***********************************************************************************************
      * \fn void on_start_clicked() and one more
      *  \brief Funciones responsables de iniciar el proceso de análisis y paralo
      * ***********************************************************************************************
      */
    void on_start_clicked();
    void on_stop_clicked();

    /** ***********************************************************************************************
      * \fn void fichero_leido(int, int, int, int)
      *  \brief Función responsable de capturar los datos de los hilos de lectura
      *  \param sample   muestra que se ha leído
      *  \param chrom    cromosoma que se ha leído
      *  \param inicio   posición inicial de la muestra leída
      *  \param final    posición final de la muestra leída
      * ***********************************************************************************************
      */
    void fichero_leido(int, int, int, int);

    /** ***********************************************************************************************
      * \fn void cromosoma_leido(int)
      *  \brief Función responsable de controlar la elctura e identificación de DMRs por cromosoma
      *  \param chrom    cromosoma que se ha leído
      * ***********************************************************************************************
      */
    void cromosoma_leido(int);

    /** ***********************************************************************************************
      * \fn void refGen_worker_acabado(ulong)
      *  \brief Función responsable de recibir el número de genes definidos para un cromosoma dato
      *  \param ref_genes   número de genes
      * ***********************************************************************************************
      */
    void refGen_worker_acabado(ulong);

    /** ***********************************************************************************************
      * \fn void on_grouped_samples_stateChanged(int) and one more
      *  \brief Funciones responsables de definir tipo de análisis -> muestras agrupadas o individuales
      *  \param arg1   activado / desactivado
      * ***********************************************************************************************
      */
    void on_grouped_samples_stateChanged(int arg1);
    void on_single_samples_stateChanged(int arg1);

    /** ***********************************************************************************************
      * \fn void on_genome_reference_currentIndexChanged(int index)
      *  \brief Función responsable de seleccionar genoma de referencia
      *  \param index   indice en la tabla de genomas de referencia
      * ***********************************************************************************************
      */
    void on_genome_reference_currentIndexChanged(int index);


private:
    Ui::HPG_Dhunter *ui;

    // control de tiempos de ejecución por bloques
    chrono::system_clock::time_point start_1;
    chrono::system_clock::time_point start_2;
    chrono::system_clock::time_point start_3;

    /** ***********************************************************************************************
      *  \brief variables responsables de capturar nombre del ficheros a analizar
      *  \param fichero             string con nombre directorio seleccionado
      *  \param ficheros_case       stringlist con nombre de todos los directorios seleccionados como casos
      *  \param ficheros_control    stringlist con nombre de los directorios seleccionados como control
      *  \param chrom_sizes         lista de los tamaños de cada uno de los cromosomas analizados
      *  \param fichero_gff         string con nombre de fichero en formato gff
      *  \param region_gff          contador de DMRs para numerar en fichero gff
      * ***********************************************************************************************
      */
    QString     fichero;
    QStringList ficheros_case;
    QStringList ficheros_control;
    QList <int> chrom_sizes;
    QString     fichero_gff;
    uint        region_gff;

    /** ***********************************************************************************************
      *  \brief variables para control de datos de cromosoma y hardware
      *  \param memory_available    cantidad de memoria GPU disponible en el PC para controlar capacidad
      * ***********************************************************************************************
      */
    int memory_available;

    /** ***********************************************************************************************
      *  \brief variables para control ventana de visualización de ficheros a analizar
      *  \param case_files          señal para habilitar el tratamiento de coloreado en ventana de casos
      *  \param *cursor_case        puntero a la línea en la ventana de casos para colorear y capturar su info
      *  \param control_files       señal para habilitar el tratamiento de coloreado en ventana de controles
      *  \param *cursor_control     puntero a la línea en la ventana de controles para colorear y capturar su info
      *  \param color               objeto para asignar color a fondo de línea
      *  \param color_char          objeto para asignar color a texto
      * ***********************************************************************************************
      */
    bool             case_files;
    bool             control_files;
    QTextCursor      *cursor_case;
    QTextCursor      *cursor_control;
    QTextBlockFormat color;
    QTextCharFormat  color_char;

    /** ***********************************************************************************************
      *  \brief variables para selección de parámetros con los que realizar el análisis
      *  \param _mc                 selecciona análisis por metilación
      *  \param _hmc                selecciona análisis por hidroximetilación
      *  \param _forward            selecciona análisis de ficheros forward
      *  \param _reverse            selecciona análisis de ficheros reverse
      *  \param _all_chroms         selecciona análisis de todos los cromosomas
      *  \param _grouped_samples    selecciona análisis por mUestras agrupadas
      *  \param _single_samples     selecciona análisis por muestras individuales
      *  \param chrom-list          listado de cromosomas a analizar
      *  \param _mc_min_coverage    valor de mínima cobertura para análisis por metilación
      *  \param _hmc_min_coverage   valor de mínima cobertura para análisis por hidroximetilación
      *  \param _threshold          valor de umbral para identificación de DMRs
      *  \param _dmr_dwt_level      valor del nivel de transformada para identificación de DMRs
      * ***********************************************************************************************
      */
    bool  _mc;
    bool  _hmc;
    bool  _forward;
    bool  _reverse;
    bool  _all_chroms;
    bool  _grouped_samples;
    bool  _single_samples;
    QList <int> chrom_list;
    int   _mc_min_coverage;
    int   _hmc_min_coverage;
    float _threshold;
    int   _dmr_dwt_level;


    /** ***********************************************************************************************
      *  \brief variables para control de posiciones extremas de fichero a analizar
      *  \param limite_inferior posición menor metilada en fichero de datos
      *  \param limite_superior posición mayor metilada en fichero de datos
      * ***********************************************************************************************
      */
    uint limite_inferior;
    uint limite_superior;


    /** ***********************************************************************************************
      *  \brief variables para búsqueda y muestra de DMRs
      *  \param **dmr_diff      datos de diferencias
      *  \param dmr_diff_cols   número de valores por vector con los que buscar DMRs por columna
      *  \param num_genes       número de genes conocidos en el cromosoma analizado
      *  \param dmrs            lista de todas las posiciones DMRs encontradas
      * ***********************************************************************************************
      */
    float       *dmr_diff;
    uint        dmr_diff_cols;
    ulong       num_genes;
    QStringList dmrs;


    /** ***********************************************************************************************
      *  \brief variables para control ventana de visualización de ficheros a analizar
      *  \param directorio  controla si se ha seleccionado el primer fichero para guardar path
      *  \param path        guarda el último path del que se ha cargado un fichero
      * ***********************************************************************************************
      */
    bool    directorio;
    QString path;


    /** ***********************************************************************************************
      *  \brief variables para control de datos por muestras y resultados de transformación en GPU
      *  \param mc          matriz con datos de metilación, cobertura y conteo por muestra y posición
      *  \param mc_grouped  matriz con datos de metilaxión, cobertura y conteo con muestras agrupadas
      *  \param h_haar_C    matriz de recepción de resultados de transformación wavelet
      *  \param posicion_metilada   acumulación de posiciones metiladas para validar DMR
      * ***********************************************************************************************
      */
    vector<vector<vector<double>>> mc;
    vector<vector<vector<double>>> mc_grouped;
    vector<vector<float>>         h_haar_C;
    vector<vector<uint>>          posicion_metilada;

    /** ***********************************************************************************************
      *  \brief variables para control de directorios y parámetros a analizar
      *  \param lista_casos     listado de los directorios correspondientes a los casos a analizar
      *  \param lista_control   listado de los directorios correspondientes a los controles a analizar
      *  \param parametros      listado de los parámetros con los que realizar el análisis
      *  \param lista_chroms    lista de los cromosomas a analizar por muestra
      * ***********************************************************************************************
      */
    QStringList lista_casos;
    QStringList lista_control;
    QStringList parametros;
    QList<int>  lista_chroms;

    /** ***********************************************************************************************
      *  \brief variables para control de procesos en hilos
      *  \param hilo_files_worker   vector de hilos que albergan la función de lectura y procesamiento previo
      *  \param files_worker        vector de funciones de lectura y procesamiento previo de ficheros
      *  \param *hilo_refGen        hilo que alberga la función de lectura de genes por cromosoma
      *  \param *refgen_worker      función de lectura de genes por cromosoma
      * ***********************************************************************************************
      */
    QVector<QThread*>      hilo_files_worker;
    QVector<Files_worker*> files_worker;
    QThread               *hilo_refGen;
    RefGen                *refGen_worker;

    /** ***********************************************************************************************
      * \fn void lectura_acabada()
      *  \brief función responsable de la identificación de DMRs y guardado en disco de los resultados
      * ***********************************************************************************************
      */
    void lectura_acabada();

    /** ***********************************************************************************************
      *  \brief variable de control de acceso a memoria compartida
      *  \param mutex   control de acceso a memoria compartida por los hilos
      * ***********************************************************************************************
      */
    QMutex mutex;


    /** ***********************************************************************************************
      * \fn void find_dmrs() and two more
      *  \brief Funciones responsables de encontrar DMRs y guardar los resultados
      * ***********************************************************************************************
      */
    void find_dmrs();
    void hallar_dmrs();
    void save_dmr_list(int);

    /** ***********************************************************************************************
      *  \brief variable para control de evolución del programa por barra de progreso
      *  \param contador    aumenta su valor conforme avanza la lectura y procesamiento.
      * ***********************************************************************************************
      */
    int contador;

    /** ***********************************************************************************************
      * \fn void enabling_wodgets (bool arg)
      *  \brief Función responsable de habilitar o no los controles de manera conjunta
      * ***********************************************************************************************
      */
    void enabling_widgets (bool arg);

};

#endif // HPG_Dhunter_H
