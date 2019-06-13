#ifndef FILES_WORKER_H
#define FILES_WORKER_H

#include <QObject>
#include <QFile>
#include <QVector>
#include <QMutex>

using namespace std;

class Files_worker : public QObject
{
    Q_OBJECT

public:
    Files_worker(QObject *parent = nullptr);

    /**
     * \fn void solicitud_lectura(QStringList, QStringList, QStringList, vector<vector<vector<float>>> &, QMutex &)
     * @brief Solicita al worker que comience
     * @param cases_files   ruta del ejecutable
     * @param control_files opciones para la ejecución
     * @param parameters    ruta para guardas los ficheros mapeados
     * @param &mcx          matriz de datos por posición y muestra
     * @param &mutexx       control de acceso a memoria compartida
     */
    void solicitud_lectura(QStringList cases_files,
                           QStringList control_files,
                           QStringList parameters,
                           vector<vector<vector<double>>> &mcx,
                           QMutex &mutexx);

    /**
     * @brief Solicita al worker que se detenga
     */
    void abort();

signals:
    /**
     * @fn void lectura_solicitada()
     * @brief Esta señal se emite cuando se le solicita al proceso que se active
     */
    void lectura_solicitada();

    /**
     * @fn void fichero_leido(int, int, int, int)
     * @brief Esta señal se emite cuando se acaba de leer un fichero de chromosoma
     * @param  sample   informa del orden que ocupa la muestra leída en la lista. empieza con el '0'
     * @param  chrom    informa del cromosoma que se ha leído
     * @param  caso     flag que indica si el fichero leído es caso '0' o control '1'
     * @param  &fichero pasa referencia del vector donde ha guardado los datos que acaba de leer

     */
    void fichero_leido(int sample, int chrom, int inicio, int final);//, QVector<QString> &fichero);

    /**
     * @fn void finished()
     * @brief Esta señal se emite cuando el proceso termina o se aborta
     */
    void finished();

public slots:
    /**
     * @fn void lectura()
     * @brief ejecuta el trabajo de mapeado sobre todos los elementos de la lista
     */
    void lectura();

private:
    /**
     * @brief variables internas para control de operaciones y almacenamiento de datos en local
     * @param aborted           señal de control de hilo activo
     * @param working           señal de control de hilo trabajando
     * @param lista_casos       listado de muestras etiquetadas como caso
     * @param lista_controles   listado de muestras etiquetadas como control
     * @param argumentos        lista de argumentos desde hilo primcipal para carga de ficheros
     * @param data              acceso a fichero en disco para lectura
     */
    bool aborted;
    bool working;
    QStringList lista_casos;
    QStringList lista_controles;
    QStringList argumentos;
    QFile data;

    /**
     * @brief matriz donde se almacenan los datos leídos por cromosoma
     */
    vector<vector<vector<double>>> *mc;

    /**
     * @brief variable de control de acceso a memoria compartida para todos los hilos
     */
    QMutex *mutex;

};

#endif // FILES_WORKER_H
