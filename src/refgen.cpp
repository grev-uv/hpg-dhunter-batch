#include "refgen.h"

#include <QFile>
#include <QDebug>
#include <sstream>

RefGen::RefGen(QObject *parent)
    : QObject(parent)
{
    aborted = false;
    working = false;
}

// ************************************************************************************************
void RefGen::solicitud_lectura(datos_cuda &cuda_datax,
                               int chromx)
{
    cuda_data = &cuda_datax;
    chrom     = chromx;

    aborted    = false;
    working    = true;

    emit lectura_solicitada();
}

// ************************************************************************************************
void RefGen::abort()
{
    if (working)
        aborted = true;
}

// ************************************************************************************************
void RefGen::lectura()
{
    // carga la matriz de referencias genómicas por posición correspondiente al cromosoma elegido
    // --------------------------------------------------------------------------------------------
    vector<string> auxRef1;
    vector<vector<string>> auxRef2;
    QFile refGene;
    QString line  = "";
    string number = "";
    refGene.setFileName(":/refGen/genmap/refmap_ucsc_chr" + QString::number(chrom) + ".csv");

    // comprueba que el fichero se ha abierto correctamente
    if (!refGene.open(QIODevice::ReadOnly | QIODevice::Text))
        qDebug() << "No se ha podido abrir el fichero de referencias genómicas";

    // lee los datos del fichero y los carga en la matriz
    while (!refGene.atEnd())
    {
        // lectura de los datos de una línea
        line = refGene.readLine();
        stringstream posicion (line.toStdString());

        while (getline (posicion, number, '\t'))
            auxRef1.push_back(number);

        // guardado de los datos
        auxRef2.push_back(auxRef1);
        auxRef1.clear();
    }
    refGene.close();

    // reserva de memoria de la matriz
    // ..limpieza
    if (cuda_data->refGen != nullptr)
    {
        delete [] cuda_data->refGen[0];
        delete[] cuda_data->refGen;
    }

    // ..reserva de espacio
    cuda_data->refGen = new string*[auxRef2.size()];
    cuda_data->refGen[0] = new string[auxRef2.size() * 6];
    for (uint i = 1; i < auxRef2.size(); i++)
            cuda_data->refGen[i] = cuda_data->refGen[i - 1] + 6;

    // copia de todos los datos a la matriz de referencias genómicas
    for (uint m = 0; m < auxRef2.size(); m++)
        for (uint k = 0; k < 6; k++)
            // [0]    nombre gen 1
            // [1]    nombre gen 2
            // [2]    cromosoma
            // [3]    posición inicial
            // [4]    posición final
            // [5]    distancia desde inicio cromosoma
            cuda_data->refGen [m][k] = auxRef2[m][k];

    // trabajo de lectura de ficheros finalizado
    emit terminado(auxRef2.size());

    aborted = true;
    working = false;
    emit finished();
}
