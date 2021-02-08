#include "files_worker.h"
#include <QDebug>
#include <sstream>
#include <iostream>

Files_worker::Files_worker(QObject *parent)
    : QObject(parent)
{
    aborted = false;
    working = false;
}

// ************************************************************************************************
void Files_worker::solicitud_lectura(QStringList cases_files,
                                     QStringList control_files,
                                     QStringList parametros,
                                     vector<vector<vector<double>>> &mcx,
                                     QMutex &mutexx)
{
    lista_casos     = cases_files;
    lista_controles = control_files;
    // parametros recibidos:
    //  0   forward bool
    //  1   reverse bool
    //  2   cromosoma
    //  3   número de fichero asignado (hilo)
    argumentos      = parametros;
    mc              = &mcx;
    mutex           = &mutexx;

    aborted         = false;
    working         = true;

    emit lectura_solicitada();
}

// ************************************************************************************************
void Files_worker::abort()
{
    if (working)
        aborted = true;
}

// ************************************************************************************************
void Files_worker::lectura()
{
    int inicio = 100000000;
    int final  = 0;

    int que_leo;
    QStringList leer = {"forward_", "reverse_", "mix_"};
    if (argumentos[0].toInt() && argumentos[1].toInt())
        que_leo = 2;
    else if (argumentos[0].toInt())
        que_leo = 0;
    else {
        que_leo = 1;
    }

    vector<vector<double>> aux3;        // vector auxiliar para montar fichero de proporción de metilación para forward

    // inicializa la posición inferior y superior
    QString fichero       = "";
    QString linea         = "";
    string numero         = "";         // dato de cada muestra en la posición de línea leida
    vector<int> aux1;                   // vector auxiliar para lectura de fichero
    vector<double> aux2;                // vector auxiliar para proporcion de metilación por fichero
    double cobertura_mC   = 0.0;
    double cobertura_hmC  = 0.0;
    double metilado       = 0.0;
    double h_metilado     = 0.0;
    double proporcion_mC  = 0.0;
    double proporcion_hmC = 0.0;

    // abre el fichero correspondiente para leer y almacenar
    if (argumentos[3].toInt() - lista_casos.size() < 0)
        fichero = lista_casos[argumentos[3].toInt()] +
                  "/methylation_map_" +
                  leer[que_leo] +
                  argumentos.at(2) +
                  ".csv";
    else
        fichero = lista_controles[argumentos[3].toInt() - lista_casos.size()] +
                  "/methylation_map_" +
                  leer[que_leo] +
                  argumentos.at(2) +
                  ".csv";

    data.setFileName(fichero);

    // comprueba que el fichero se ha abierto correctamente
    if (!data.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        qDebug() << "ERROR opening file: " << fichero;
    }
    else
    {
        // lee y guarda todos los datos
        while (!data.atEnd())
        {
            if (aborted)
                break;

            linea = data.readLine();
            stringstream posicion (linea.toStdString());

            while (getline (posicion, numero, ' '))
                if (isdigit(numero[0]))
                    aux1.push_back(stoi(numero));

            // procesamiento de los datos de la línea
            // como cobertura se suma en número de C y mC o hmC
            // se descarta el número de nC
            //cobertura_mC  = aux1[1] + aux1[2] + aux1[3];
            cobertura_mC  = aux1[1] + aux1[3];
            metilado      = aux1[3];
            //cobertura_hmC = aux1[1] + aux1[2] + aux1[4];
            cobertura_hmC = aux1[1] + aux1[4];
            h_metilado    = aux1[4];

            if (cobertura_mC > 0)
                proporcion_mC = metilado / cobertura_mC;
            else
            {
                proporcion_mC = 0.0;
                metilado      = 0.0;
                cobertura_mC  = 0.0;
            }

            if (cobertura_hmC > 0)
                proporcion_hmC = h_metilado / cobertura_hmC;
            else
            {
                proporcion_hmC = 0.0;
                h_metilado     = 0.0;
                cobertura_hmC  = 0.0;
            }

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
            // 10   muestra (posición en la lista de caso o control)
            // 11   caso/control (0/1)
            // 12   forward/reverse (0/1)
            aux2.push_back(aux1[0]);
            aux2.push_back(proporcion_mC);
            aux2.push_back(cobertura_mC);
            aux2.push_back(aux1[1]);
            aux2.push_back(aux1[2]);
            aux2.push_back(aux1[3]);
            aux2.push_back(aux1[4]);
            aux2.push_back(proporcion_hmC);
            aux2.push_back(cobertura_hmC);
            aux2.push_back(argumentos[2].toInt());
            aux2.push_back((argumentos[3].toInt() - lista_casos.size() < 0) ? argumentos[3].toInt() : argumentos[3].toInt() - lista_casos.size());
            aux2.push_back((argumentos[3].toInt() < lista_casos.size()) ? 0 : 1 );
            aux2.push_back(que_leo);

            // si la primera posición es cero no se contempla para preservar la integridad de
            // la identificación de DMRs tal y como está definido
            if (aux1[0] > 0)
                aux3.push_back(aux2);

            aux1.clear();
            aux2.clear();
        }

        // actualiza la posición mínima y máxima
        if (inicio >= int(aux3.front().front()))
            inicio = int(aux3.front().front());
        if (final < int(aux3.back().front()))
            final = int(aux3.back().front());
    }

    // cierra el fichero de datos
    data.close();

    // carga los datos en la matriz principal
    mutex->lock();
    mc->push_back(aux3);
    mutex->unlock();

    // envía señal de lectura de fichero para su procesado en otro hilo
    emit fichero_leido(argumentos[3].toInt(), argumentos.at(2).toInt(), inicio, final);

    // trabajo de lectura de ficheros finalizado
    aborted = true;
    working = false;
    data.close();

    emit finished();
}
