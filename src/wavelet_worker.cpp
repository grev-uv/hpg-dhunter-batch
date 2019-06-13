#include "wavelet_worker.h"
#include <sstream>

Wavelet_worker::Wavelet_worker(QObject *parent)
    : QObject(parent)
{
    aborted = false;
    working = false;
}

void Wavelet_worker::solicitud_lectura(QStringList parametros,
                                       vector<vector<vector<float>>> &mcx,
                                       QVector<QVector<QString>> &datox)
{
    // parametros recibidos:
    //  0   hilo (fila de la matriz de datos)
    //  1   cromosoma
    //  2   sample (posición en la lista)
    //  3   caso/control 0/1
    argumentos = parametros;
    mc         = &mcx;
    datos      = &datox;

    aborted    = false;
    working    = true;

    emit lectura_solicitada();
}

void Wavelet_worker::abort()
{
    if (working)
        aborted = true;
}

void Wavelet_worker::lectura()
{
    string numero = "";                // dato de cada muestra en la posición de línea leida
    vector<int> aux1;                   // vector auxiliar para lectura de fichero
    vector<float> aux2;                 // vector auxiliar para proporcion de metilación por fichero
    vector<vector<float>> aux3;         // vector auxiliar para montar fichero de proporción de metilación mc global
    float cobertura_mC   = 0.0;
    float cobertura_hmC  = 0.0;
    float metilado       = 0.0;
    float h_metilado     = 0.0;
    float proporcion_mC  = 0.0;
    float proporcion_hmC = 0.0;

    // carga las posiciones y los datos del fichero en la matriz de datos particular
    for (int i = 0; i < datos->at(argumentos.at(0).toInt()).size(); i++)
    {
        if (aborted)
            break;

        // lectura de los datos de una línea
        mutex.lock();
        stringstream posicion (datos->at(argumentos.at(0).toInt())[i].toStdString());
        mutex.unlock();

        while (getline (posicion, numero, ' '))
            aux1.push_back(stoi(numero));

        // procesamiento de los datos de la línea en función de si se ha elegido
        // analizar el ratio de metilación o de hidroximetilación
        cobertura_mC  = aux1[1] + aux1[2] + aux1[3];
        metilado      = aux1[3];
        cobertura_hmC = aux1[1] + aux1[2] + aux1[4];
        h_metilado    = aux1[4];

        if (cobertura_mC > 0)
            proporcion_mC = metilado / cobertura_mC;

        if (cobertura_hmC > 0)
            proporcion_hmC = h_metilado / cobertura_hmC;

        // guardado de los datos de un cromosoma de una muestra de un sentido
        // 0    posición en el cromosoma
        // 1    proporción de mC frente a la cobertura
        // 2    cobertura de mC (C + noC + mC)reads
        // 3    número de reads identificando una C
        // 4    número de reads identificando una no C
        // 5    número de reads identificando una mC
        // 6    número de reads identificando una hmC
        // 7    proporción de hmC frente a cobertura
        // 8    cobertura de hmC (C + noC + hmC)reads
        // 9    chromosoma
        // 10   muestra (posición en la lista)
        // 11   caso/control (0/1)
        aux2.push_back(aux1[0]);
        aux2.push_back(proporcion_mC);
        aux2.push_back(cobertura_mC);           // incluye la cobertura para selección de datos posterior
        aux2.push_back(aux1[1]);                // incluye número de C en la posición
        aux2.push_back(aux1[2]);                // incluye número de no C en la posición
        aux2.push_back(aux1[3]);                // incluye número de mC en la posición
        aux2.push_back(aux1[4]);                // incluye número de hmC en la posición
        aux2.push_back(proporcion_hmC);         // incluye proporción hmC
        aux2.push_back(cobertura_hmC);          // incluye cobertura hmC
        aux2.push_back(argumentos[1].toInt());  // incluye número de cromosoma
        aux2.push_back(argumentos[2].toInt());  // incluye indice de muestra en la lista
        aux2.push_back(argumentos[3].toInt());  // incluye caso/control
        aux3.push_back(aux2);

        aux1.clear();
        aux2.clear();
    }

    // carga de datos del fichero en la matriz de cobertura global
    mutex.lock();
    mc->push_back(aux3);
    mutex.unlock();

    emit terminado(argumentos.at(0).toInt());

    // trabajo de lectura de ficheros finalizado
    aborted = true;
    working = false;
    emit finished();
}
