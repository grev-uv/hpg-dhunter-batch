#ifndef REFGEN_H
#define REFGEN_H

#include <QObject>
#include <data_pack.h>

using namespace std;

class RefGen : public QObject
{
    Q_OBJECT

public:
    RefGen(QObject *parent = nullptr);

    /**
     * @fn void solicitud_lectura(datos_cuda &, int)
     *  @brief Solicita al worker que comience
     *  @param &cuda-datax   referencia a struct con datos para GPU
     *  @param chrom         número de cromosoma
     */
    void solicitud_lectura(datos_cuda &cuda_datax,
                           int chromx);

    void abort();

signals:
    /**
     * @fn void lectura_solicitada()
     *  @brief Esta señal se emite cuando se le solicita al proceso que se active
     */
    void lectura_solicitada();

    /**
     * @fn finished()
     *  @brief Esta señal se emite cuando el proceso termina o se aborta
     */
    void finished();
    void terminado(ulong);

public slots:
    /**
     * @fn void lectura()
     *  @brief ejecuta el trabajo de mapeado sobre todos los elementos de la lista
     */
    void lectura();

private:
    /**
     * @brief variables internas para control de operaciones y almacenamiento de datos en local
     * @param aborted       señal de control de hilo activo
     * @param working       señal de control de hilo trabajando
     * @param *cuda_data    estructura con datos de control globales donde guardar listado de genes
     * @param chrom         cromosoma que se está analizando
     */
    bool aborted;
    bool working;
    datos_cuda *cuda_data;
    int chrom;
};

#endif // REFGEN_H
