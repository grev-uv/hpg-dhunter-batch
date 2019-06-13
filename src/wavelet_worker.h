#ifndef WAVELET_WORKER_H
#define WAVELET_WORKER_H

#include <QObject>
#include <QVector>
#include <QMutex>

using namespace std;

class Wavelet_worker : public QObject
{
    Q_OBJECT

public:
    Wavelet_worker(QObject *parent = nullptr);


    /**
     * @brief Solicita al worker que comience
     * @param cases_files   ruta del ejecutable
     * @param control_files opciones para la ejecución
     * @param parameters    ruta para guardas los ficheros mapeados
     * @param chroms        lista de cromosomas a leer
     */
    void solicitud_lectura(QStringList parametros,
                           vector<vector<vector<float>>> &mcx,
                           QVector<QVector<QString>> &datox);

    /**
     * @brief Solicita al worker que se detenga
     */
    void abort();

signals:
    /**
     * @brief Esta señal se emite cuando se le solicita al proceso que se active
     */
    void lectura_solicitada();

    /**
     * @brief Esta señal se emite cuando el proceso termina o se aborta
     */
    void finished();
    void terminado(int);

public slots:
    /**
     * @brief ejecuta el trabajo de mapeado sobre todos los elementos de la lista
     */
    void lectura();

private:
    bool aborted;
    bool working;
    QStringList argumentos;
    QVector<QVector<QString>> *datos;
    vector<vector<vector<float>>> *mc;

    QMutex mutex;
};

#endif // WAVELET_WORKER_H
