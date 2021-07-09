#include "hpg_dhunter.h"
#include "ui_hpg_dhunter.h"
#include <QFileDialog>
#include <QDebug>
#include <QFile>
#include <QTimer>
#include <QProcess>
#include <QRegularExpression>
#include <QMessageBox>
#include <QDesktopServices>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <exception>

using namespace std;

HPG_Dhunter::HPG_Dhunter(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::HPG_Dhunter),
    mutex()
{
    ui->setupUi(this);

    // inicialización de variables ------------------------------------------------------------
    fichero            = "";
    num_genes          = 0;
    limite_inferior    = 500000000;
    limite_superior    = 0;
    cuda_data.mc_full  = nullptr;
    cuda_data.h_haar_C = nullptr;
    cuda_data.refGen   = nullptr;
    dmr_diff           = nullptr;

    contador           = 0;
    ui->progressBar->setMinimum(0);

    // cursor para ventana con lista de ficheros de control
    cursor_control = new QTextCursor();
    control_files  = false;

    // cursor para ventana con lista de ficheros de casos
    cursor_case = new QTextCursor();
    case_files  = false;

    // mantiene el último path en el explorador de ficheros
    path       = "";
    directorio = false;

    // incialización variables de selección de análisis
    _mc              = true;
    _hmc             = false;
    _forward         = true;
    _reverse         = true;
    _all_chroms      = true;

    // inicialización de variables para cálculo de DMRs
    _mc_min_coverage   = ui->mC_cobertura->value();
    _hmc_min_coverage  = ui->hmC_cobertura->value();
    _threshold      = float(ui->threshold->value() * 0.01);
    _dmr_dwt_level  = ui->dmr_dwt_level->value();

    // inicialización de opciones
    ui->all_chroms->setEnabled(false);
    ui->selected_chrms->setChecked(true);

    // comprueba la memoria disponible en la tarjeta gráfica para controlar los ficheros a cargar
    // ..captura la información suministrada por el comando "nvidia-smi"
    QProcess p;
    p.start("nvidia-smi");
    p.waitForFinished();
    QString data = p.readAllStandardOutput();
    p.close();

    // ..busca los datos que concuerdan con "[[:digit:]]+MiB" y se queda con el segundo dato
    //   que informa de la capacidad total de memoria de la tarjeta
    QRegularExpression re("(\\d+)MiB");
    QRegularExpressionMatchIterator i = re.globalMatch(data);
    QRegularExpressionMatch match_1 = i.next();
    QRegularExpressionMatch match_2 = i.next();
    qDebug() << "memoria GPU ocupada / total ----> " << match_1.captured(1) << "/" << match_2.captured(1);

    // ..se asigna a la variable el valor en MiB
    memory_available = match_2.captured(1).toInt() - match_1.captured(1).toInt();

    ui->statusBar->showMessage("System available GPU RAM: " + QString::number(memory_available));

}

// ************************************************************************************************
HPG_Dhunter::~HPG_Dhunter()
{
    cuda_end(cuda_data);
    delete [] cuda_data.mc_full[0];
    delete [] cuda_data.mc_full;
    delete ui;
}


// ************************************************************************************************
// **************************************ZONA CARGA CASOS******************************************
// ************************************************************************************************
void HPG_Dhunter::on_case_file_clicked()
{
    // abre ventana de explorador de directorios para seleccionar
    // el directorio donde se encuentran los cromosomas de un caso
    fichero = QFileDialog::getExistingDirectory( this,
                                                 tr("Select a case folder"),
                                                 (directorio) ? path : QDir::homePath(),
                                                 QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks
                                               );

    if(fichero.isEmpty() || fichero.isNull())
        fichero = "";
    else
    {
        case_files = false;

        ui->case_files->appendPlainText(fichero);

        QStringList lista = ui->case_files->toPlainText().split("\n");
        if (lista.size() >= 1)
        {
            if (QStringList(ui->control_files->toPlainText().split('\n'))[0].size() != 0)
            {
                ui->start->setEnabled(true);
                ui->start->setFocus();
            }
            else
            {
                ui->start->setEnabled(false);
                ui->stop->setEnabled(false);
            }

            ui->delete_case->setEnabled(true);
            ui->up_case->setEnabled(true);
            ui->down_case->setEnabled(true);
        }

        directorio = true;
        fichero.chop(fichero.split("/").back().size());
        path = fichero;

        // guarda listado de los ficheros cargados para analizar
        ficheros_case.clear();
        for (int i = 0; i < lista.size(); i++)
            ficheros_case.append(lista.at(i));

        qDebug() << ficheros_case.last();

        case_files = true;

    }
}

//*************************************************************************************************
void HPG_Dhunter::on_case_files_cursorPositionChanged()
{
    if (case_files)
    {
        // desmarca la línea previa quitando color de fondo
        color.setBackground(Qt::white);
        cursor_case->select(QTextCursor::LineUnderCursor);
        cursor_case->setBlockFormat(color);

        // adquiere el cursor de la línea seleccionada con el ratón
        *cursor_case = ui->case_files->textCursor();
        cursor_case->movePosition(QTextCursor::StartOfBlock);
        cursor_case->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);

        // marca con color de fondo la línea seleccionada
        color.setBackground(Qt::gray);
        cursor_case->select(QTextCursor::LineUnderCursor);
        cursor_case->setBlockFormat(color);
    }
}

//*************************************************************************************************
void HPG_Dhunter::on_delete_case_clicked()
{
    // variables internas
    int line_number  = cursor_case->blockNumber();                // número de línea actual del cursor
    QStringList list = ui->case_files->toPlainText().split("\n"); // lista de strings con ficheros

    // deshabilita el cambio de color por acción sobre el cursor
    case_files = false;
    // borra el fichero seleccionado
    list.removeAt(line_number);
    // limpia la lista de visualización
    ui->case_files->clear();
    // si quedan ficheros en la lista, la copia en la lista de visualización
    if (!list.isEmpty())
    {
        ui->case_files->appendPlainText(list.join('\n'));
        case_files = true;
    }
    else
    {
        // si no hay ficheros deshabilita el botón de ejecutar alineamiento
        ui->start->setEnabled(false);
        ui->stop->setEnabled(false);
    }

    // aplica fondo blanco a línea bajo el cursor
    color.setBackground(Qt::white);
    cursor_case->select((QTextCursor::LineUnderCursor));
    cursor_case->setBlockFormat(color);

    // mueve el cursor al principo de la lista de visualización
    cursor_case->movePosition(QTextCursor::Start);
    // traslada el cursor a la posición que tenía antes de borrar
    for (int i = 0; i < line_number; i++)
    {
        cursor_case->movePosition(QTextCursor::StartOfBlock);
        cursor_case->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
        cursor_case->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);

        if (i == list.size() - 2)
            break;
    }

    // actualiza la lista con la posición del cursor adecuada para que se resalte
    ui->case_files->setTextCursor(*cursor_case);

    // guarda listado de los ficheros_case cargados para analizar
    QStringList lista = ui->case_files->toPlainText().split("\n");
    ficheros_case.clear();
    for (int i = 0; i < lista.size(); i++)
        if (lista.at(i).size() > 0)
            ficheros_case.append(lista.at(i).split("/").at(lista.at(i).split("/").length() - 2));
}

//*************************************************************************************************
void HPG_Dhunter::on_up_case_clicked()
{
    if (cursor_case->blockNumber() > 0)
    {
        // variables internas
        int line_number  = cursor_case->blockNumber();                  // número de línea del cursor
        QStringList list = ui->case_files->toPlainText().split("\n");   // lista de ficheros
        QString line     = list.at(line_number);                        // fichero en línea seleccionada

        // deshabilita el cambio de color por acción sobre el cursor
        case_files = false;
        // borra el fichero seleccionado
        list.removeAt(line_number);
        // inserta el fichero selecionado en una posición más arriba
        list.insert(line_number - 1, line);
        // limpia la lista de visualización
        ui->case_files->clear();
        // si quedan ficheros en la lista, la copia en la lista de visualización
        if (!list.isEmpty())
        {
            case_files = true;
            ui->case_files->appendPlainText(list.join('\n'));
        }

        // aplica fondo blanco a línea bajo el cursor
        color.setBackground(Qt::white);
        cursor_case->select((QTextCursor::LineUnderCursor));
        cursor_case->setBlockFormat(color);

        // mueve el cursor al principio de la lista de visualización
        cursor_case->movePosition(QTextCursor::Start);
        // traslada el cursor a la posición que ocupa el fichero trasladado
        for (int i = 0; i < line_number - 1; i++)
        {
            cursor_case->movePosition(QTextCursor::StartOfBlock);
            cursor_case->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
            cursor_case->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);
        }

        // actualiza la lista con la posición del cursor adecuada para que se resalte
        ui->case_files->setTextCursor(*cursor_case);
    }
}

//*************************************************************************************************
void HPG_Dhunter::on_down_case_clicked()
{
    // variables internas
    QStringList list = ui->case_files->toPlainText().split("\n");   // lista de ficheros

    // proceder con el desplazamiento abajo si no es el último elemento de la lista
    if (cursor_case->blockNumber() < list.size())
    {
        // variables internas
        int line_number = cursor_case->blockNumber();   // número de línea del cursor
        QString line    = list.at(line_number);         // fichero en línea seleccionada

        // deshabilita el cambio de color por acción sobre el cursor
        case_files = false;
        // borra el fichero seleccionado
        list.removeAt(line_number);
        // inserta el fichero seleccionado en una posición más abajo
        list.insert(line_number + 1, line);
        // limpia la lista de visualización
        ui->case_files->clear();
        // se quedan ficheros en la lista, la copia en la lista de visualización
        if (!list.isEmpty())
        {
            case_files = true;
            ui->case_files->appendPlainText(list.join('\n'));
        }

        // aplica fondo blanco a línea bajo el cursor
        color.setBackground(Qt::white);
        cursor_case->select((QTextCursor::LineUnderCursor));
        cursor_case->setBlockFormat(color);

        // mueve el cursor al principio de la lista de visualización
        cursor_case->movePosition(QTextCursor::Start);
        // traslada el cursor a la posición que ocupa el fichero trasladado
        for (int i = 0; i <= line_number; i++)
        {
            cursor_case->movePosition(QTextCursor::StartOfBlock);
            cursor_case->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
            cursor_case->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);

            // comprueba que no ha alcanzado el final de la lista
            if (i == list.size() - 2)
                break;
        }

        // actualiza la lista con la posición del cursor adecuada para que se resalte
        ui->case_files->setTextCursor(*cursor_case);
    }
}


// ************************************************************************************************
// **************************************ZONA CARGA CONTROLES**************************************
// ************************************************************************************************
void HPG_Dhunter::on_control_file_clicked()
{
    // abre ventana de explorador de directorios para seleccionar
    // el directorio donde se encuentran los cromosomas de un caso
    fichero = QFileDialog::getExistingDirectory( this,
                                              tr("Select a control folder"),
                                              (directorio) ? path : QDir::homePath() ,
                                              QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks
                                            );

    if(fichero.isEmpty() || fichero.isNull())
        fichero = "";
    else
    {
        control_files = false;

        ui->control_files->appendPlainText(fichero);

        QStringList lista = ui->control_files->toPlainText().split("\n");
        if (lista.size() >= 1)
        {
            if (QStringList(ui->case_files->toPlainText().split('\n'))[0].size() != 0)
            {
                ui->start->setEnabled(true);
                ui->start->setFocus();
            }
            else
            {
                ui->start->setEnabled(false);
                ui->stop->setEnabled(false);
            }

            ui->delete_control->setEnabled(true);
            ui->up_control->setEnabled(true);
            ui->down_control->setEnabled(true);
        }

        directorio = true;
        fichero.chop(fichero.split("/").back().size());
        path = fichero;

        // guarda listado de los ficheros cargados para analizar
        ficheros_control.clear();
        for (int i = 0; i < lista.size(); i++)
            ficheros_control.append(lista.at(i).split("/").back());

        qDebug() << ficheros_control.last();

        control_files = true;
    }
}

//*************************************************************************************************
void HPG_Dhunter::on_control_files_cursorPositionChanged()
{
    if (control_files)
    {
        // desmarca la línea previa quitando color de fondo
        color.setBackground(Qt::white);
        cursor_control->select(QTextCursor::LineUnderCursor);
        cursor_control->setBlockFormat(color);

        // adquiere el cursor de la línea seleccionada con el ratón
        *cursor_control = ui->control_files->textCursor();
        cursor_control->movePosition(QTextCursor::StartOfBlock);
        cursor_control->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);

        // marca con color de fondo la línea seleccionada
        color.setBackground(Qt::gray);
        cursor_control->select(QTextCursor::LineUnderCursor);
        cursor_control->setBlockFormat(color);
    }
}

//*************************************************************************************************
void HPG_Dhunter::on_delete_control_clicked()
{
    // variables internas
    int line_number  = cursor_control->blockNumber();                // número de línea actual del cursor
    QStringList list = ui->control_files->toPlainText().split("\n"); // lista de strings con ficheros

    // deshabilita el cambio de color por acción sobre el cursor
    control_files = false;
    // borra el fichero seleccionado
    list.removeAt(line_number);
    // limpia la lista de visualización
    ui->control_files->clear();
    // si quedan ficheros en la lista, la copia en la lista de visualización
    if (!list.isEmpty())
    {
        ui->control_files->appendPlainText(list.join('\n'));
        control_files = true;
    }
    else
    {
        // si no hay ficheros deshabilita el botón de ejecutar alineamiento
        ui->start->setEnabled(false);
        ui->stop->setEnabled(false);
    }

    // aplica fondo blanco a línea bajo el cursor
    color.setBackground(Qt::white);
    cursor_control->select((QTextCursor::LineUnderCursor));
    cursor_control->setBlockFormat(color);

    // mueve el cursor al principo de la lista de visualización
    cursor_control->movePosition(QTextCursor::Start);
    // traslada el cursor a la posición que tenía antes de borrar
    for (int i = 0; i < line_number; i++)
    {
        cursor_control->movePosition(QTextCursor::StartOfBlock);
        cursor_control->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
        cursor_control->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);

        if (i == list.size() - 2)
            break;
    }

    // actualiza la lista con la posición del cursor adecuada para que se resalte
    ui->control_files->setTextCursor(*cursor_control);

    // guarda listado de los ficheros_control cargados para analizar
    QStringList lista = ui->control_files->toPlainText().split("\n");
    ficheros_control.clear();
    for (int i = 0; i < lista.size(); i++)
        if (lista.at(i).size() > 0)
            ficheros_control.append(lista.at(i).split("/").at(lista.at(i).split("/").length() - 2));
}

//*************************************************************************************************
void HPG_Dhunter::on_up_control_clicked()
{
    if (cursor_control->blockNumber() > 0)
    {
        // variables internas
        int line_number  = cursor_control->blockNumber();                 // número de línea del cursor
        QStringList list = ui->control_files->toPlainText().split("\n");  // lista de ficheros
        QString line     = list.at(line_number);                          // fichero en línea seleccionada

        // deshabilita el cambio de color por acción sobre el cursor
        control_files = false;
        // borra el fichero seleccionado
        list.removeAt(line_number);
        // inserta el fichero selecionado en una posición más arriba
        list.insert(line_number - 1, line);
        // limpia la lista de visualización
        ui->control_files->clear();
        // si quedan ficheros en la lista, la copia en la lista de visualización
        if (!list.isEmpty())
        {
            control_files = true;
            ui->control_files->appendPlainText(list.join('\n'));
        }

        // aplica fondo blanco a línea bajo el cursor
        color.setBackground(Qt::white);
        cursor_control->select((QTextCursor::LineUnderCursor));
        cursor_control->setBlockFormat(color);

        // mueve el cursor al principio de la lista de visualización
        cursor_control->movePosition(QTextCursor::Start);
        // traslada el cursor a la posición que ocupa el fichero trasladado
        for (int i = 0; i < line_number - 1; i++)
        {
            cursor_control->movePosition(QTextCursor::StartOfBlock);
            cursor_control->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
            cursor_control->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);
        }

        // actualiza la lista con la posición del cursor adecuada para que se resalte
        ui->control_files->setTextCursor(*cursor_control);
    }
}

//*************************************************************************************************
void HPG_Dhunter::on_down_control_clicked()
{
    // variables internas
    QStringList list = ui->control_files->toPlainText().split("\n");   // lista de ficheros

    // proceder con el desplazamiento abajo si no es el último elemento de la lista
    if (cursor_control->blockNumber() < list.size())
    {
        // variables internas
        int line_number = cursor_control->blockNumber();  // número de línea del cursor
        QString line    = list.at(line_number);           // fichero en línea seleccionada

        // deshabilita el cambio de color por acción sobre el cursor
        control_files = false;
        // borra el fichero seleccionado
        list.removeAt(line_number);
        // inserta el fichero seleccionado en una posición más abajo
        list.insert(line_number + 1, line);
        // limpia la lista de visualización
        ui->control_files->clear();
        // se quedan ficheros en la lista, la copia en la lista de visualización
        if (!list.isEmpty())
        {
            control_files = true;
            ui->control_files->appendPlainText(list.join('\n'));
        }

        // aplica fondo blanco a línea bajo el cursor
        color.setBackground(Qt::white);
        cursor_control->select((QTextCursor::LineUnderCursor));
        cursor_control->setBlockFormat(color);

        // mueve el cursor al principio de la lista de visualización
        cursor_control->movePosition(QTextCursor::Start);
        // traslada el cursor a la posición que ocupa el fichero trasladado
        for (int i = 0; i <= line_number; i++)
        {
            cursor_control->movePosition(QTextCursor::StartOfBlock);
            cursor_control->movePosition(QTextCursor::EndOfBlock, QTextCursor::KeepAnchor);
            cursor_control->movePosition(QTextCursor::Down, QTextCursor::MoveAnchor);

            // comprueba que no ha alcanzado el final de la lista
            if (i == list.size() - 2)
                break;
        }

        // actualiza la lista con la posición del cursor adecuada para que se resalte
        ui->control_files->setTextCursor(*cursor_control);
    }
}




// ************************************************************************************************
// *************************************ZONA EJECUCION PROCESO*************************************
// ************************************************************************************************

// HILO LECTURA DE DATOS DE FICHEROS
// ************************************************************************************************
void HPG_Dhunter::fichero_leido(int sample, int chrom, int inicio, int final)
{
    // paso de los datos leídos al hilo de procesamiento de datos
    qDebug() << "fichero - datos disponibles: " << contador << sample << " " << chrom;

    mutex.lock();
    if (limite_inferior > uint(inicio))
        limite_inferior = uint(inicio);
//    mutex.unlock();
//    mutex.lock();
    if (limite_superior < uint(final))
        limite_superior = uint(final);
//    mutex.unlock();

    // contador de evolución de lectura y análisis
//    mutex.lock();
    contador ++;
    mutex.unlock();
    ui->progressBar->setValue(ui->progressBar->value() + 1);

    if (contador % (lista_casos.size() + lista_control.size()) == 0)
        cromosoma_leido(chrom);
}

// ************************************************************************************************
void HPG_Dhunter::cromosoma_leido(int chrom)
{
    QEventLoop loop;

    QTimer::singleShot(500, &loop, SLOT(quit()));
    loop.exec();

    STOP_TIMER_2("FICHEROS CROMOSOMA LEIDOS -------------");
    START_TIMER_3 // proceso de cálculo

    ui->statusBar->showMessage("identifying DMRs...");

    // calcula wavelet de cada uno de los ficheros leídos
    lectura_acabada();

    STOP_TIMER_3("DMRs CALCULADOS -------")

    limite_inferior = 500000000;
    limite_superior = 0;

    // índice en la lista de cromosomas, del cromosoma leído
    int idx = lista_chroms.indexOf(chrom);

    // si no ha terminado de leer los ficheros de la lista de cromosomas -> lanza nuevo hilo
    if (idx < lista_chroms.size() - 1)
    {
        START_TIMER_2;

        ui->statusBar->showMessage("freeing memory... it will takes a while");

        // limpia las matrices de datos del cromosoma anterior
        vector<vector<vector<double>>>().swap(mc);

        ui->statusBar->showMessage("reading next chromosome...");

        // borra todos los hilos creados
        foreach(QThread *i, hilo_files_worker)
            delete i;

        hilo_files_worker.clear();
        hilo_files_worker.shrink_to_fit();
        files_worker.clear();
        files_worker.shrink_to_fit();


        // hilo y conexiones para el proceso de lectura de ficheros
        for ( int i = 0; i < (lista_casos.size() + lista_control.size()); i++)
        {
            hilo_files_worker.append(new QThread());
            files_worker.append(new Files_worker());
        }


        for (int i = 0; i < hilo_files_worker.size(); i++)
        {
            files_worker[i]->moveToThread(hilo_files_worker[i]);
            connect(files_worker[i], SIGNAL(fichero_leido(int, int, int, int)), SLOT(fichero_leido(int, int, int, int)));
            connect(hilo_files_worker[i], &QThread::finished, files_worker[i], &QObject::deleteLater);
            files_worker[i]->connect(hilo_files_worker[i], SIGNAL(started()), SLOT(lectura()));
            hilo_files_worker[i]->connect(files_worker[i],SIGNAL(lectura_solicitada()), SLOT(start()));
            hilo_files_worker[i]->connect(files_worker[i], SIGNAL(finished()), SLOT(quit()), Qt::DirectConnection);

            // arranque del hilo de lectura
            if (hilo_files_worker[i]->isRunning())
                hilo_files_worker[i]->wait();

            // se le asigna el número de cromosoma
            parametros[2] = QString::number(lista_chroms.at(idx + 1));
            parametros[3] = QString::number(i);
            qDebug() << "cromosoma a leer:" << parametros[2] << parametros[3];

            // se lanza el primer hilo de lectura de ficheros por cromosoma
            files_worker[i]->solicitud_lectura(lista_casos, lista_control, parametros, mc, mutex);
        }

        // se lanza el hilo de carga de referencias genéticas, si se dispone de ellas
        switch (ui->genome_reference->currentIndex())
        {
            case 0:
                break;
            case 1:
                hilo_refGen   = new QThread();
                refGen_worker = new RefGen();
                refGen_worker->moveToThread(hilo_refGen);
                connect(refGen_worker, SIGNAL(terminado(ulong)), SLOT(refGen_worker_acabado(ulong)));
                connect(hilo_refGen, &QThread::finished, refGen_worker, &QObject::deleteLater);
                refGen_worker->connect(hilo_refGen, SIGNAL(started()), SLOT(lectura()));
                hilo_refGen->connect(refGen_worker,SIGNAL(lectura_solicitada()), SLOT(start()));
                hilo_refGen->connect(refGen_worker, SIGNAL(finished()), SLOT(quit()), Qt::DirectConnection);

                refGen_worker->solicitud_lectura(cuda_data, parametros[2].toInt());
                break;
            default:
                ;
        }
    }
    else
    {
        STOP_TIMER_1("PROCESO TERMINADO---------------")

        // lectura de ficheros acabada
        qDebug() << "se han leído todos los ficheros " << mc.size() << "x" << mc.at(mc.size()-1).size();

        QMessageBox::information(this,
                                 tr("CSV to DMRs app"),
                                 tr("The DMRs identification is finished")
                                );

        ui->start->setEnabled(true);
        ui->start->setFocus();
        ui->stop->setEnabled(false);
        enabling_widgets(true);
    }
}


// ************************************************************************************************
void HPG_Dhunter::refGen_worker_acabado(ulong num_genex)
{
    num_genes = num_genex;

    ui->statusBar->showMessage("chromosome: " + parametros[2] +
                               ", with " + QString::number(num_genes) + " known genes");
}


// GESTOR DE EJECUCIÓN
// ************************************************************************************************
void HPG_Dhunter::on_start_clicked()
{
    lista_casos   = ui->case_files->toPlainText().split("\n");
    lista_control = ui->control_files->toPlainText().split("\n");

    qDebug() << "---control muestras con cobertura--- " << lista_casos.length() <<
                ui->min_covSamples_x_region->value() * 0.01 <<
                lista_casos.length() * (ui->min_covSamples_x_region->value() * 0.01);


    // inicializa contador de regiones DMR detectadas
    region_gff = 0;

    // inicializa los parámetros a enviar a los hilos
    parametros = (QStringList() << QString::number(_forward) << // se informa forward reads 0/1
                                QString::number(_reverse) <<    // se informa reverse reads 0/1
                                "0" <<                          // se informa del número de cromosoma
                                "0"                             // se informa del número de hilo asignado
                 );

    switch (ui->genome_reference->currentIndex())
    {
        case 0:
            lista_chroms = {};
            break;
        case 1:
            lista_chroms = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
            break;
        default:
            lista_chroms = {};
    }

    // genera la lista de cromosomas a analizar filtrando números y eliminando duplicados
    if (!_all_chroms)
    {
        lista_chroms.clear();
        QString lista = ui->chromosomes_list->text();
        lista.replace(QRegularExpression("[a-zA-Z,.-]+"), " ");
        QStringList llista = lista.split(QRegularExpression("\\s+"));
        foreach(auto n, llista)
        {
            if (n.toInt() > 0 && n.toInt() < 25)
                if (!lista_chroms.contains(n.toInt()))
                    lista_chroms << n.toInt();
        }
    }

    qDebug() << lista_chroms.isEmpty() << " : " <<
                lista_casos.front() << " - " <<
                lista_control.front() << " - " <<
                parametros << " - " <<
                lista_chroms;

    // ventanas de precaución ante falta de datos previos al proceso
    // --------------------------------------------------------------------------------------------
    // falta la ruta donde guardar los ficheros mapeados CSV
    if (ui->out_path_label->text().isEmpty())
    {
        QMessageBox::warning(this,
                             tr("CSV to DMRs app"),
                             tr("The path to save the DMRs founded files is empty\n"
                                "Please, select a folder to save them")
                            );

        ui->start->setEnabled(true);
        ui->start->setFocus();
        ui->stop->setEnabled(false);

        return;
    }

    // no hay cromosomas que analizar
    if (lista_chroms.isEmpty())
    {
        QMessageBox::warning(this,
                             tr("CSV to DMRs app"),
                             tr("The chromosome list to analyze is empty\n"
                                "Please, fill the list with desired values")
                            );

        ui->start->setEnabled(true);
        ui->start->setFocus();
        ui->stop->setEnabled(false);

        return;
    }

    // verifica la lista de cromosomas a analizar después de filtrada
    QString lista_c;
    foreach(int n, lista_chroms)
    {
        lista_c.append(QString::number(n) + ", ");
    }
    QMessageBox::StandardButton reply;
    reply = QMessageBox::question(this,
                                  "CSV to DMRs app",
                                  "The chromosome list to analyze is:\n" +
                                  lista_c +
                                  "\nIs it right?",
                                  QMessageBox::Yes|QMessageBox::No);

    // confirmación de inicio de proceso de mapeado
    if (reply == QMessageBox::Yes)
    {
        ui->start->setEnabled(false);
        ui->stop->setEnabled(true);
        enabling_widgets(false);

        // inicialización de barra de progreso de trabajo
        contador = 0;
        ui->progressBar->setValue(0);
        ui->progressBar->setMaximum((lista_casos.size() + lista_control.size()) * lista_chroms.size() + (lista_chroms.size() * (_mc + _hmc)));

        qDebug() << "yastamos panalizar";

        // limpia las matrices de datos
        ui->statusBar->showMessage("freeing memory... it will takes a while");

        vector<vector<vector<double>>>().swap(mc);

        // borra todos los posibles hilos creados anteriormente
        foreach(QThread *i, hilo_files_worker)
            delete i;

        QVector<QThread*>().swap(hilo_files_worker);
        QVector<Files_worker*>().swap(files_worker);

        START_TIMER_1 // total trabajo
        START_TIMER_2 // por cromosoma

        ui->statusBar->showMessage("loading files...");

        // crea nuevos hilos para lectura de ficheros
        for (int i = 0; i < lista_casos.size() + lista_control.size(); i++)
        {
            hilo_files_worker.append(new QThread());
            files_worker.append(new Files_worker());
        }

        // moviendo workers a los hilos y conexiones para el proceso de lectura de ficheros
        for (int i = 0; i < hilo_files_worker.size(); i++)
        {
            files_worker[i]->moveToThread(hilo_files_worker[i]);
            connect(files_worker[i], SIGNAL(fichero_leido(int, int, int, int)), SLOT(fichero_leido(int, int, int, int)));
            connect(hilo_files_worker[i], &QThread::finished, files_worker[i], &QObject::deleteLater);
            files_worker[i]->connect(hilo_files_worker[i], SIGNAL(started()), SLOT(lectura()));
            hilo_files_worker[i]->connect(files_worker[i],SIGNAL(lectura_solicitada()), SLOT(start()));
            hilo_files_worker[i]->connect(files_worker[i], SIGNAL(finished()), SLOT(quit()), Qt::DirectConnection);

            // arranque del hilo de lectura
            if (hilo_files_worker[i]->isRunning())
                hilo_files_worker[i]->wait();

            // se le asigna el número de cromosoma
            parametros[2] = QString::number(lista_chroms.at(0));
            parametros[3] = QString::number(i);
            qDebug() << "cromosoma a leer:" << parametros[2] << parametros[3];

            // se lanza el primer hilo de lectura de ficheros por cromosoma
            files_worker[i]->solicitud_lectura(lista_casos, lista_control, parametros, mc, mutex);
        }


        // se lanza el hilo de carga de referencias genéticas, si se dispone de ellas
        switch (ui->genome_reference->currentIndex())
        {
        case 0:
            break;
        case 1:
            hilo_refGen   = new QThread();
            refGen_worker = new RefGen();
            refGen_worker->moveToThread(hilo_refGen);
            connect(refGen_worker, SIGNAL(terminado(ulong)), SLOT(refGen_worker_acabado(ulong)));
            connect(hilo_refGen, &QThread::finished, refGen_worker, &QObject::deleteLater);
            refGen_worker->connect(hilo_refGen, SIGNAL(started()), SLOT(lectura()));
            hilo_refGen->connect(refGen_worker,SIGNAL(lectura_solicitada()), SLOT(start()));
            hilo_refGen->connect(refGen_worker, SIGNAL(finished()), SLOT(quit()), Qt::DirectConnection);

            refGen_worker->solicitud_lectura(cuda_data, parametros[2].toInt());
            break;
        default:
            ;
        }
    }
    else
        return;

}

// ************************************************************************************************
// ************************************************************************************************
// HILO DE PROCESAMIENTO DE DATOS LEIDOS
// ************************************************************************************************
// ************************************************************************************************
void HPG_Dhunter::lectura_acabada()
{
    // procesar el análisis de GPU aqui,
    //    pasar resultado para identificación de DMRs y guardado de fichero a otro hilo
    //    devolver control a lectura de ficheros de siguiente cromosoma

    // incialización de variables
    // --------------------------------------------------------------------------------------------
    cuda_data.h_haar_L.clear();     // vector con número de datos por nivel
    cuda_data.pitch          = 0;   // ajuste óptimo de memoria GPU para datos de cada muestra
    cuda_data.pitch_2        = 0;   // ajuste óptimo de memoria GPU para auxiliar
    cuda_data.sample_num     = 0;   // número de datos por muestra
    cuda_data.samples        = 0;   // número de muestras a trasnformar
    cuda_data.levels         = 0;   // número de niveles a transformar
    cuda_data.data_adjust    = 0;   // ajuste desfase en división por nivel para número impar de datos
    cuda_data.rango_inferior = 0;   // límite inferior ventana de datos a transformar
    cuda_data.rango_superior = 1;

    // libera la memoria de la GPU
    cuda_end(cuda_data);

    // cálculo de la dimensión total del cromosoma leído
    // la dimensión o número de posiciones totales que sea número par
    uint dimension = limite_superior - limite_inferior + 1;
    if ((dimension & 0x01) == 1)
        dimension++;

    // realiza la operación de carga en GPU y de identificación de DMRs para mC y hmC seleccionadas
    // mh = 0 -> analiza mC si está seleccionado este análisis
    // mh = 1 -> analiza hmC si está seleccionado este análisis
    for (int mh = 0; mh < 2; mh++)
    {
        if ((!mh && _mc) || (mh && _hmc))
        {
            // limpia matriz de resultados de procesamiento en GPU
            vector<vector<float>>().swap(h_haar_C);

            // realiza el cálculo de DWT en GPU
            // selecciona bloques de filas de mc para procesar en GPU hasta que procesa toda la matriz mc
            uint filas_procesadas = 0;
            uint filas_a_GPU      = 1;

            vector<vector<uint>>(uint(mc.size()), vector<uint>()).swap(posicion_metilada);

            while (filas_procesadas < mc.size())
            {
                // borra la memoria utilizada por cuda_data.mc_full
                if (cuda_data.mc_full != nullptr)
                {
                    delete [] cuda_data.mc_full[0];
                    delete [] cuda_data.mc_full;
                }

                // borra la memoria utilizada por cuda_data.h_haar_C
                if (cuda_data.h_haar_C != nullptr)
                {
                    delete [] cuda_data.h_haar_C[0];
                    delete [] cuda_data.h_haar_C;
                }

                // calcula el tamaño de la matriz de datos
                filas_a_GPU  = 1;
                uint tamanyo = dimension * filas_a_GPU * sizeof(float) / (1024 * 1024);  // tamaño en MiB

                // calcula el número de filas de mc que puede procesar en GPU simultaneamente
                while (tamanyo < 0.5 * memory_available)
                {
                    if (filas_a_GPU + filas_procesadas < mc.size())
                        filas_a_GPU++;
                    else
                        break;

                    tamanyo = dimension * filas_a_GPU * sizeof(float) / (1024 * 1024);  // tamaño en MiB
                }

                filas_procesadas += filas_a_GPU;

                qDebug() << "-----  hola 3 " << "- filas procesadas / GPU" << filas_procesadas << "/" << filas_a_GPU;

                // actualiza estructura de datos
                cuda_data.samples        = int(filas_a_GPU);                        // número de ficheros a analizar
                cuda_data.sample_num     = dimension;                               // cantidad de datos por fichero
                cuda_data.rango_inferior = 0;                                       // primer valor cromosoma
                cuda_data.rango_superior = limite_superior - limite_inferior;       // último valor
                cuda_data.levels         = ui->dmr_dwt_level->value();              // número de niveles a transformar
                cuda_data.data_adjust    = 0;                                       // ajuste desfase en división por nivel para número impar de datos
                cuda_data.h_haar_L.clear();                                         // vector con número de datos por nivel

                // crea matriz ampliada -------------------------------------------------------------------
                //      -> vectores con todas las posiciones contiguas
                //      -> con ceros en las posiciones sin metilación
                // reserva TODA la memoria CONTIGUA con todos los datos de todas las muestras
                // para trasvase de datos entre GPU y CPU con CUDA, la matriz debe ser contigua completa
                // reserva la memoria para la matriz de datos extendida
                cuda_data.mc_full = new float*[cuda_data.samples];
                cuda_data.mc_full[0] = new float[uint(cuda_data.samples) * cuda_data.sample_num]();
                for (int i = 1; i < cuda_data.samples; i++)
                        cuda_data.mc_full[i] = cuda_data.mc_full[i - 1] + cuda_data.sample_num;

                // copia de todos los datos a la matriz ampliada
                // --------------------------------------------------------------------------------------------
                //for (uint m = filas_procesadas - filas_a_GPU; m < filas_procesadas; m++)
                //uint cuenta_posiciones_metiladas;
                //uint pos_met;
                for (uint m = 0; m < uint(cuda_data.samples); m++)
                {
                    // rellenar con datos las posiciones metiladas si la cobertura es mayor que el umbral
                    uint posicion = m + filas_procesadas - filas_a_GPU;

                    for (uint k = 0; k < mc[posicion].size(); k++)
                    {
                        if(mc[posicion][k][(mh == 0 ? 2 : 8)] >= (mh == 0 ? _mc_min_coverage : _hmc_min_coverage))
                        {
                            cuda_data.mc_full [m][uint(mc[posicion][k][0]) - limite_inferior] = float(mc[posicion][k][(mh == 0 ? 1 : 7)]);

                            posicion_metilada[posicion].push_back(uint(mc[posicion][k][0] - limite_inferior));
                        }
                    }
                }

                // envía los datos a la memoria global de la GPU
                // --------------------------------------------------------------------------------------------
                // libera la memoria de la GPU
                cuda_end(cuda_data);

                // envía el total de los datos a la GPU
                cuda_send_data(cuda_data);

                // procesado de los datos
                cuda_calculo_haar_L(cuda_data);
                cuda_main(cuda_data);

                // recoge los resultados en una matriz, acumulando todos los resultados
                vector<float> aux(uint(cuda_data.h_haar_L[0]), 0.0);
                for (int i = 0; i < cuda_data.samples; i++)
                {
                    for (size_t j = 0; j < size_t(cuda_data.h_haar_L[0]); j++)
                        aux[j] = cuda_data.h_haar_C[i][j];

                    h_haar_C.push_back(aux);
                }
            }

            // comprueba la memoria disponible en la tarjeta gráfica para controlar los ficheros a cargar
            // ..captura la información suministrada por el comando "nvidia-smi"
            QProcess p;
            p.start("nvidia-smi");
            p.waitForFinished();
            QString data = p.readAllStandardOutput();
            p.close();

            // ..busca los datos que concuerdan con "[[:digit:]]+MiB" y se queda con el segundo dato
            //   que informa de la capacidad de memoria disponible de la tarjeta
            QRegularExpression re("(\\d+)MiB");
            QRegularExpressionMatchIterator i = re.globalMatch(data);
            QRegularExpressionMatch match_1 = i.next();
            QRegularExpressionMatch match_2 = i.next();
            qDebug() << "memoria GPU ocupada / total ----> " << match_1.captured(1) << "/" << match_2.captured(1);

            // ..se asigna a la variable el valor en MiB
            memory_available = match_2.captured(1).toInt() - match_1.captured(1).toInt();

            qDebug() << "tamaño final matriz de datos h_haar_C: " << h_haar_C.size() << "x" << h_haar_C.at(0).size()
                     << " y pos_met:" << posicion_metilada.size() << posicion_metilada.at(0).size();

            // con los resultados completos en la matriz, pasa a la identificación de DMRs
            find_dmrs();

            // con los DMRs identificados, salva el resultado en un fichero
            save_dmr_list(mh);

            ui->progressBar->setValue(ui->progressBar->value() + 1);
        }

        // libera la memoria de la GPU y la memoria RAM
        cuda_end(cuda_data);
        if (cuda_data.mc_full != nullptr)
        {
            delete [] cuda_data.mc_full[0];
            delete [] cuda_data.mc_full;
        }
        cuda_data.mc_full = nullptr;

        if (cuda_data.h_haar_C != nullptr)
        {
            delete [] cuda_data.h_haar_C[0];
            delete [] cuda_data.h_haar_C;
        }
        cuda_data.h_haar_C = nullptr;
    }
}


// ************************************************************************************************
void HPG_Dhunter::find_dmrs()
{
    uint numero_casos;
    uint numero_control;

    // reserva la matriz de diferencias de medias por grupos de control y casos
    // y la llena de ceros
    if (dmr_diff != nullptr)
        delete[] dmr_diff;
    dmr_diff = new float[cuda_data.h_haar_L[0]];
    for (int i = 0; i < cuda_data.h_haar_L[0]; i++)
        dmr_diff[i] = 0.0;

    uint paso = uint(pow(2, ui->dmr_dwt_level->value()));

    qDebug() << "----------- buscando DMRs por muestras individuales ----------";
    uint contador = 0;
    uint cont_diff = 0;

    int ultimo_m = 0;

    vector<uint> idx_pos_met (uint(mc.size()), 0);

    // realiza el cálculo de medias de las muestras de control y de los casos con cobertura sobre umbral
    for (uint m = 0; m < uint(cuda_data.h_haar_L[0]); m++) // ...en cada posición
    {
        float media_casos   = 0.0;
        float media_control = 0.0;
        numero_casos        = 0;
        numero_control      = 0;

        for (uint i = 0; i < h_haar_C.size(); i++)
        {
            uint aux_1 = 0;
            uint aux_2 = 0;

            while (m * paso >= posicion_metilada[i][idx_pos_met[i] + aux_1] && idx_pos_met[i] + aux_1 < posicion_metilada[i].size() - 1)
                aux_1++;

            aux_2 = aux_1;

            while ((m + 1) * paso > posicion_metilada[i][idx_pos_met[i] + aux_2] && idx_pos_met[i] + aux_2 < posicion_metilada[i].size() - 1)
                aux_2++;

            idx_pos_met[i] += aux_2;

            if (aux_2 - aux_1 >= paso * uint(ui->min_CpG_x_region->value()) * 0.01)
            {
                if (int(mc[i][0][11]) == 0)
                {
                    media_control += h_haar_C[i][m];
                    numero_control++;
                }
                else
                {
                    media_casos += h_haar_C[i][m];
                    numero_casos++;
                }
            }
        }

        //if (numero_casos > 0 && numero_control > 0)                                                                    // al menos una muestra por grupo tiene cobertura
        if (numero_casos   >= (uint(lista_casos.length()   * (ui->min_covSamples_x_region->value() * 0.01))) &&
            numero_control >= (uint(lista_control.length() * (ui->min_covSamples_x_region->value() * 0.01))))             // al menos el XX% por grupo tienen cobertura
        //if (numero_casos == uint(lista_casos.length()) && numero_control == uint(lista_control.length()))              // todas las muestras tienen cobertura
        {
            // guada diferencias solo si caso y control han resultado diferentes de cero -> hay cobertura mínima en, al menos, una muestra de caso y control
            dmr_diff[m] = (media_casos / numero_casos) - (media_control / numero_control);

            // para control de programa (a borrar)
            ultimo_m = int(m);
            contador++;


        //    if (contador > 4000)
        //        qDebug() << dmr_diff[m];


            if (dmr_diff[m] < -_threshold || dmr_diff[m] > _threshold)
            {
                cont_diff++;
            }
        }
    }
    qDebug() << "número de ventanas dwt con valor > 0 " << contador
             << " número de ventanas con diff mayor que umbral: " << cont_diff
             << " ultimo eme: " << ultimo_m
             << " diferencia último: " << dmr_diff[ultimo_m];


    dmr_diff_cols = uint(cuda_data.h_haar_L[0]);

    // llama a función de cálculo de diferencias entre muestras para ponerlas en la tabla
    hallar_dmrs();

}

// ************************************************************************************************
void HPG_Dhunter::hallar_dmrs()
{
    QString linea = "";
    // calcula el número de posiciones de cromosoma que hay por cada posición del vector DWT calculado
    uint paso     = uint(pow(2, _dmr_dwt_level));
    // crea array con las posiciones iniciales de cada tramo en el cromosoma con DM superior al umbral
    vector<uint> posicion_dmr (uint(cuda_data.h_haar_L[0]), 0);

    // encuentra DMRs en función del threshold establecido -----------------------------
    // rellena las posiciones con diferencias válidas
    // si el valor de la difercia es menor que el umbral, la posición se queda con valor 0
    for (int m = 0; m < cuda_data.h_haar_L[0]; m++)
        if (dmr_diff[m] < -_threshold || dmr_diff[m] > _threshold)
            posicion_dmr[uint(m)] = uint(m) * paso + limite_inferior;

    // buscar y rellenar la lista de DMRs
    dmrs.clear();


    int inicio = 0;
    int fin    = int(num_genes);
    for (uint p = 0; p < uint(cuda_data.h_haar_L[0]); p++)
    {
        if (posicion_dmr[p] >= limite_inferior)
        {
            linea.clear();
            uint q = p;     // para ayuda en la zona de detección de referencia de genoma

            // busca las posiciones inicial y final de la DMR
            //---------------------------------------------------------------------------
            linea.append(QString::number(posicion_dmr[q]));

            while (p + 1 < dmr_diff_cols && posicion_dmr[p + 1] >= limite_inferior)
               p++;

            linea.append("-" + QString::number(posicion_dmr[p] + paso));


            // búsqueda del nombre del GEN implicado o más cercano a los DMRs encontrados
            //---------------------------------------------------------------------------
            // se realiza sobre datos de la genome.ucsc.edu data base sobre genes conocidos
            // ..previamente se han cargado los nombres y posiciones de los genes correspondientes
            // al cromosoma que se está analizando
            // ..por búsqueda binaria sobre este fichero se determina el nombre del gen.
            switch (ui->genome_reference->currentIndex())
            {
                case 0:
                    break;

                case 1:
                    int mitad        = inicio;
                    bool match       = false;
                    uint gen_ini     = 0;
                    uint gen_ant_fin = uint(stoul(cuda_data.refGen[0][4]));

                    while (uint(stoul(cuda_data.refGen[mitad][3])) < posicion_dmr[q] && mitad < fin - 1)
                        mitad++;

                    gen_ini = uint(stoul(cuda_data.refGen[mitad][3]));
                    if (mitad > 0)
                        gen_ant_fin = uint(stoul(cuda_data.refGen[mitad - 1][4]));

                    // el inicio dmr es igual que inicio del gen
                    if (gen_ini == posicion_dmr[q])
                    {
                        match = true;
                        linea.append(" " + QString::fromStdString(cuda_data.refGen[mitad][0]) +
                                     " " + QString::fromStdString(cuda_data.refGen[mitad][1]) +
                                     " 0");
                    }
                    // el inicio dmr es menor que inicio del gen pero el final dmr es mayor que el inicio del gen
                    else if (gen_ini <= posicion_dmr[p] + paso)
                    {
                        match = true;
                        linea.append(" " + QString::fromStdString(cuda_data.refGen[mitad][0]) +
                                     " " + QString::fromStdString(cuda_data.refGen[mitad][1]) +
                                     " -" + QString::number(stoul(cuda_data.refGen[mitad][3]) > posicion_dmr[q] ?
                                                            stoul(cuda_data.refGen[mitad][3]) - posicion_dmr[q] :
                                                            posicion_dmr[q] - stoul(cuda_data.refGen[mitad][3])));
                    }
                    // el inicio dmr es mayor que inicio del gen anterior pero es menor que el fin del gen anterior
                    else if (gen_ant_fin > posicion_dmr[q])
                    {
                        match = true;
                        linea.append(" " + QString::fromStdString(cuda_data.refGen[mitad > 0? mitad - 1 : 0][0]) +
                                     " " + QString::fromStdString(cuda_data.refGen[mitad > 0? mitad - 1 : 0][1]) +
                                     " +" + ((posicion_dmr[q] > stoul(cuda_data.refGen[mitad > 0? mitad - 1 : 0][3])) ?
                                            QString::number(posicion_dmr[q] - stoul(cuda_data.refGen[mitad > 0? mitad - 1 : 0][3])) :
                                            QString::number(stoul(cuda_data.refGen[mitad > 0? mitad - 1 : 0][3]) - posicion_dmr[q])));
                    }


                    // si se encuentra entre genes, ver de qué gen está más cerca
                    // se elige la distancia más pequeña entre:
                    // ..distancia inicio dmr y fin gen anterior
                    // ..distancia fin dmr e inicio gen posterior
                    if (!match)
                    {
                        ulong dif1 = posicion_dmr[q] - gen_ant_fin;
                        ulong dif2 = gen_ini - posicion_dmr[p] + paso;

                        if (dif1 >= dif2)
                        {
                            linea.append(" " + QString::fromStdString(cuda_data.refGen[mitad][0]) +
                                         " " + QString::fromStdString(cuda_data.refGen[mitad][1]) +
                                         " --" + QString::number(dif2));
                        }
                        else
                        {
                            linea.append(" " + QString::fromStdString(cuda_data.refGen[mitad > 0? mitad - 1 : 0][0]) +
                                         " " + QString::fromStdString(cuda_data.refGen[mitad > 0? mitad - 1 : 0][1]) +
                                         " ++" + QString::number(dif1));
                        }
                    }

                    if (mitad > 0)
                        inicio = mitad - 1;
                    break;
            }

            // define si está hipermetilado o hipometilado el caso frente al control
            //-----------------------------------------------------------------------
            linea.append((dmr_diff[q] < 0)? " hipo":" hiper");

            // añade resultado de análisis DWT
            //-----------------------------------------------------------------------
            linea.append(" " + QString::number(double(dmr_diff[q])));
            linea.append("//" + QString::number(q) + " " + QString::number(p));

            // añade la información a la lista de DMRs
            //-----------------------------------------------------------------------
            dmrs.append(linea);
        }
    }

    // informa de proceso en barra inferior -----------------------------------------
    ui->statusBar->showMessage(QString::number(dmrs.size()) + " DMRs found");

    // pasa la lista de dmrs a la ventana de visualización
    qDebug() << "tamaño del fichero que guarda los dmrs localizados: " << dmrs.size() << "x" << (dmrs.size() ? dmrs.at(0).size() : 0);
}


// ************************************************************************************************
void HPG_Dhunter::save_dmr_list(int mh)
{
    // prepara nombre de fichero y directorio para guardar la lista de dmrs
    fichero = ui->out_path_label->text() +
              "/chromosome_" + QString::number(int(mc[0][0][9])) + "_" +
              (mh ? "hmc_thr0" : "mc_thr0") + QString::number(ui->threshold->value()) +
              "_dwt" + QString::number(ui->dmr_dwt_level->value()) +
              "_cov" + (mh ? ui->hmC_min_cov->text() : ui->mC_min_cov->text()) + ".csv";

    QFile data;
    data.setFileName(fichero);

    // prepara nombre de fichero para guardar los datos con formato GFF
    fichero_gff = ui->out_path_label->text() + "/" + ui->out_path_label->text().split("/").last() + "_" +
                  (mh ? "hmc_thr0" : "mc_thr0") + QString::number(ui->threshold->value()) +
                  "_dwt" + QString::number(ui->dmr_dwt_level->value()) +
                  "_cov" + (mh ? ui->hmC_min_cov->text() : ui->mC_min_cov->text()) + ".gff";

    QFile data_gff;
    data_gff.setFileName(fichero_gff);
    bool gff_open;


    // comprueba que el fichero se ha abierto correctamente
    if (!data.open(QIODevice::WriteOnly))
    {
        QMessageBox::warning(this,
                             "ERROR Opening files",
                             "An error occurred opening the file: " + fichero +
                             "\nPlease, check the file for corrupted"
                            );
        qDebug() << "ERROR opening file: " << fichero;
        return;
    }
    else
    {
        QTextStream s(&data);
        if (dmrs.size() > 0)
        {
            QString linea, linea_detail;
            uint pos_inf;
            uint pos_sup;
            uint ancho_dmr;
            vector <uint> posicion_muestra (mc.size(), 0);

            // comprueba que el fichero para guardar información en formato GFF se abre correctamente
            if (!data_gff.open(QIODevice::WriteOnly | QIODevice::Append))
            {
                gff_open = false;
                qDebug() << "ERROR al abrir el fichero con formato GFF";
            }
            else
                gff_open = true;


            // encabezado de la información del dmr
            switch (ui->genome_reference->currentIndex())
            {
            case 0:
                s << "pos_init-pos_end methylation dwt_diff\n";
                break;
            case 1:
                s << "pos_init-pos_end name_1 name_2 distance methylation dwt_diff\n";
                break;
            }

            qDebug() << "guardando datos en ficheros";

            // añade una línea por dmr detectaado y línea de características por muestra en cada dmr
            for (int i = 0; i < dmrs.size(); i++)
            {
//                qDebug() << "empieza escritura en fichero dmr" << i;
                // información del dmr para obtener las características de cada muestra
                linea     = dmrs.at(i);
                pos_inf   = linea.split("-")[0].toUInt();
                pos_sup   = linea.split("-")[1].split(" ")[0].toUInt();
                ancho_dmr = linea.split("-")[1].split(" ")[0].toUInt() - linea.split("-")[0].toUInt();

                //**************************************************************************************************************
                // escribe información del DMR en el fichero GFF
                if (gff_open)
                {
                    QTextStream gff(&data_gff);

                    // columna 1 y 2 (sequence , source)
                //    gff << "chr" << (mc[0][0][9] < 10 ? "0" : "") << QString::number(int(mc[0][0][9])) << "\t" << "HPG-Dhunter\t";
                    gff << "chr" << QString::number(int(mc[0][0][9])) << "\t" << "HPG-Dhunter\t";

                    region_gff++;

                    // columna 3 (feature)
                    switch (ui->genome_reference->currentIndex())
                    {
                    case 0:
                        gff << dmrs.at(i).split(" ")[1] << "\t";
                        break;
                    case 1:
                        gff << dmrs.at(i).split(" ")[4] << "\t";
                        break;
                    }

                    // columna 4 y 5 (start , end)
                    gff << QString::number(pos_inf) << "\t" <<
                           QString::number(pos_sup) << "\t";

                    // columna 6 (score)
                    switch (ui->genome_reference->currentIndex())
                    {
                    case 0:
                        gff << dmrs.at(i).split("//")[0].split(" ")[2] << "\t";
                        break;
                    case 1:
                        gff << dmrs.at(i).split("//")[0].split(" ")[5] << "\t";
                        break;
                    }

                    // columna 7 y 8 (strand , phase)
                    gff << ".\t" << ".\t";

                    // columna 9 (#samples coverage threshold dwt_level)
                    switch (ui->genome_reference->currentIndex())
                    {
                    case 0:
                        gff << "Note=DMR_Region:" << QString::number(region_gff) << ",";
                        break;
                    case 1:
                        gff << "Name=" << dmrs.at(i).split("//")[0].split(" ")[2] << ";" <<
                               "Note=Distance:" << dmrs.at(i).split("//")[0].split(" ")[3] << ",";
                        break;
                    }
                    gff << "Samples:" << QString::number(mc.size()) << "," <<
                           "Coverage:" << QString::number(mh ? ui->hmC_cobertura->value() : ui->mC_cobertura->value()) << "," <<
                           "Threshold:" << ui->threshold_label->text() << "," <<
                           "DWT_level:" << QString::number(ui->dmr_dwt_level->value()) << "," <<
                           "Density:" << QString::number(ui->min_CpG_x_region->value()) << "%,"
                           "Samples/region w/cov:" << QString::number(ui->min_covSamples_x_region->value()) << "%\n";
                }



                //**************************************************************************************************************

                // escribe zona dmr detectada en fichero particular
                s << dmrs.at(i).split("//")[0] << '\n';

                // posición inicial y final de zona dwt en h_haar_C correspondiente al DMR identificado
                uint pos_dwt_ini = dmrs.at(i).split("//")[1].split(" ")[0].toUInt();
                uint pos_dwt_fin = dmrs.at(i).split("//")[1].split(" ")[1].toUInt();
            //    qDebug() << "posición dwt inicial: " << pos_dwt_ini << "posición dwt final: " << pos_dwt_fin;

                // encabezado de las características por fichero dentro de la zona dmr
                s << " sample dwt_value ratio C_positions cov_min cov_mid cov_max sites_Cnm sites_Cnh sites_mC sites_hmC dist_min dist_mid dist_max\n";

                // rellena el fichero
                // guarda información de cada muestra de la zona dmr detectada
                //***************************************************************************************
                // busca la posición inical
                for (uint j = 0; j < mc.size(); j++)
                {
                    uint posicion = posicion_muestra[j];
                    while (pos_inf > mc[j][posicion][0] && posicion < mc[j].size() - 1)
                        posicion++;
                    posicion_muestra[j] = posicion;
                }

                for (uint j = 0; j < uint(mc.size()); j++)
                {
                    if (int(mc[j][0][11]) == 0)
                    {
                        s << " " << lista_casos.at(int(mc[j][0][10])).split("/").back() << " ";

                        int cobertura_minima = 500000000;
                        int cobertura_maxima = 0;
                        int cobertura_media  = 0;
                        int distancia_minima = 500000000;
                        int distancia_maxima = 0;
                        int distancia_media  = 0;
                        int sites_C          = 0;
                        int sites_nC         = 0;
                        int sites_mC         = 0;
                        int sites_hmC        = 0;
                        int posiciones       = 0;
                        float dwt_valor      = 0.0;
                        float ratio_medio    = 0.0;

                        linea_detail.clear();

                        // busca la posición inical
                        uint posicion = posicion_muestra[j];
                    //    while (pos_inf > mc[j][posicion][0])
                    //        posicion++;
                    //    posicion_muestra[j] = posicion;

                        // búsqueda de valores a lo largo del DMR
                        while (pos_sup > mc[j][posicion][0] && posicion < mc[j].size() - 2)
                        {
                            // cobertura
                            if (cobertura_minima >= mc[j][posicion][mh ? 8 : 2])
                                cobertura_minima = int(mc[j][posicion][mh ? 8 : 2]);
                            if (cobertura_maxima < mc[j][posicion][mh ? 8 : 2])
                                cobertura_maxima = int(mc[j][posicion][mh ? 8 : 2]);
                            cobertura_media += mc[j][posicion][mh ? 8 : 2];
                            if (mc[j][posicion][mh ? 8 : 2] > 0)
                                ratio_medio     += float(mc[j][posicion][mh ? 6 : 5] / mc[j][posicion][mh ? 8 : 2]);

                            // distancia
                            if (posicion + 2 < mc[j].size() && ancho_dmr > mc[j][posicion + 1][0] - mc[j][posicion][0])
                            {
                                if (distancia_minima >= mc[j][posicion + 1][0] - mc[j][posicion][0])
                                    distancia_minima = int(mc[j][posicion + 1][0] - mc[j][posicion][0]);
                                if (distancia_maxima < mc[j][posicion + 1][0] - mc[j][posicion][0])
                                    distancia_maxima = int(mc[j][posicion + 1][0] - mc[j][posicion][0]);
                                distancia_media += mc[j][posicion + 1][0] - mc[j][posicion][0];
                            }

                            // número de posiciones detectadas por tipo de mononucleótico
                            sites_C   += (mc[j][posicion][3] > 0) ? 1 : 0;
                            sites_nC  += (mc[j][posicion][4] > 0) ? 1 : 0;
                            sites_mC  += (mc[j][posicion][5] > 0) ? 1 : 0;
                            sites_hmC += (mc[j][posicion][6] > 0) ? 1 : 0;

                            // número de posiciones detectadas con algún tipo de nucleótido sensible
                            if (mc[j][posicion][mh ? 8 : 2] > 0)
                                posiciones++;
                        //    posiciones++;

                            posicion++;
                        //    qDebug() << j << " -> " << mc[j].size() << " - " << posicion;
                        }

                        // valor medio dwt en la región identificada
                        for (uint i = pos_dwt_ini; i <= pos_dwt_fin; i++)
                            dwt_valor += h_haar_C[j][i];                    // (h_haar_C[f][i] / (pos_dwt_fin - pos_dwt_ini + 1));
                        dwt_valor = pos_dwt_fin - pos_dwt_ini + 1 != 0 ? dwt_valor / (pos_dwt_fin - pos_dwt_ini + 1) : dwt_valor;
                    //    dwt_valor = h_haar_C[j][pos_dwt_ini];

                        // carga de resultado en línea de texto para mostrar
                        //if ((posiciones > 1 ? cobertura_media / posiciones : cobertura_media) >= (mh ? _hmc_min_coverage : _mc_min_coverage))
                        if (cobertura_maxima >= (mh ? _hmc_min_coverage : _mc_min_coverage))
                        {
                            s << QString("%1").arg(double(dwt_valor)) << " " << //::number(double(dwt_valor), 'f', 3) << " " <<
                                 QString("%1").arg(posiciones > 1 ? double(ratio_medio / posiciones) : double(ratio_medio)) << " " <<
                                 QString::number(posiciones) << " " <<
                                 QString::number(cobertura_minima >= 500000000 ? 0 : cobertura_minima) << " " <<
                                 QString::number((posiciones > 1) ? cobertura_media / posiciones : cobertura_media) << " " <<
                                 QString::number(cobertura_maxima) << " " <<
                                 QString::number(sites_C) << " " <<
                                 QString::number(sites_nC) << " " <<
                                 QString::number(sites_mC) << " " <<
                                 QString::number(sites_hmC) << " " <<
                                 QString::number(distancia_minima >= 500000000 ? 0 : distancia_minima) << " " <<
                                 QString::number((posiciones > 1) ? distancia_media / posiciones : distancia_media) << " " <<
                                 QString::number(distancia_maxima) << "\n";
                        }
                        else
                        {
                            s << QString("%1").arg(double(dwt_valor)) << " " << //::number(double(dwt_valor), 'f', 3) << " " <<
                                 QString("%1").arg(posiciones > 1 ? double(ratio_medio / posiciones) : double(ratio_medio)) << " " <<
                                 "0 0 0 0 0 0 0 0 0 0 0\n";
                        }
                    }
                }

                for (uint j = 0; j < uint(mc.size()); j++)
                {
                    if (int(mc[j][0][11]) != 0)
                    {
                        s << " " << lista_control.at(int(mc[j][0][10])).split("/").back() << " ";

                        int cobertura_minima = 500000000;
                        int cobertura_maxima = 0;
                        int cobertura_media  = 0;
                        int distancia_minima = 500000000;
                        int distancia_maxima = 0;
                        int distancia_media  = 0;
                        int sites_C          = 0;
                        int sites_nC         = 0;
                        int sites_mC         = 0;
                        int sites_hmC        = 0;
                        int posiciones       = 0;
                        float dwt_valor      = 0.0;
                        float ratio_medio    = 0.0;

                        linea_detail.clear();

                        // busca la posición inical
                        uint posicion = posicion_muestra[j];
                    //    while (pos_inf > mc[j][posicion][0])
                    //        posicion++;
                    //    posicion_muestra[j] = posicion;

                        // búsqueda de valores a lo largo del DMR
                        while (pos_sup > mc[j][posicion][0] && posicion < mc[j].size() - 2)
                        {
                            // cobertura
                            if (cobertura_minima >= mc[j][posicion][mh ? 8 : 2])
                                cobertura_minima = int(mc[j][posicion][mh ? 8 : 2]);
                            if (cobertura_maxima < mc[j][posicion][mh ? 8 : 2])
                                cobertura_maxima = int(mc[j][posicion][mh ? 8 : 2]);
                            cobertura_media += mc[j][posicion][mh ? 8 : 2];
                            if (mc[j][posicion][mh ? 8 : 2] > 0)
                                ratio_medio     += float(mc[j][posicion][mh ? 8 : 2] != 0.0 ?
                                                         mc[j][posicion][mh ? 6 : 5] / mc[j][posicion][mh ? 8 : 2] : 0);

                            // distancia
                            if (posicion + 2 < mc[j].size() && ancho_dmr > mc[j][posicion + 1][0] - mc[j][posicion][0])
                            {
                                if (distancia_minima >= mc[j][posicion + 1][0] - mc[j][posicion][0])
                                    distancia_minima = int(mc[j][posicion + 1][0] - mc[j][posicion][0]);
                                if (distancia_maxima < mc[j][posicion + 1][0] - mc[j][posicion][0])
                                    distancia_maxima = int(mc[j][posicion + 1][0] - mc[j][posicion][0]);
                                distancia_media += mc[j][posicion + 1][0] - mc[j][posicion][0];
                            }

                            // número de posiciones detectadas por tipo de mononucleótico
                            sites_C   += (mc[j][posicion][3] > 0) ? 1 : 0;
                            sites_nC  += (mc[j][posicion][4] > 0) ? 1 : 0;
                            sites_mC  += (mc[j][posicion][5] > 0) ? 1 : 0;
                            sites_hmC += (mc[j][posicion][6] > 0) ? 1 : 0;

                            // número de posiciones detectadas con algún tipo de nucleótido sensible
                            if (mc[j][posicion][mh ? 8 : 2] > 0)
                                posiciones++;
                            //posiciones++;

                            posicion++;
                        //    qDebug() << j << " -> " << mc[j].size() << " - " << posicion;
                        }

                        // valor medio dwt en la región identificada
                        for (uint i = pos_dwt_ini; i <= pos_dwt_fin; i++)
                            dwt_valor += h_haar_C[j][i];                    // (h_haar_C[f][i] / (pos_dwt_fin - pos_dwt_ini + 1));
                        dwt_valor = pos_dwt_fin - pos_dwt_ini + 1 ? dwt_valor / (pos_dwt_fin - pos_dwt_ini + 1) : dwt_valor;
                    //    dwt_valor = h_haar_C[j][pos_dwt_ini];

                        // carga de resultado en línea de texto para mostrar
                        //if ((posiciones > 1 ? cobertura_media / posiciones : cobertura_media) >= (mh ? _hmc_min_coverage : _mc_min_coverage))
                        if (cobertura_maxima >= (mh ? _hmc_min_coverage : _mc_min_coverage))
                        {
                            s << QString("%1").arg(double(dwt_valor)) << " " << //::number(double(dwt_valor), 'f', 3) << " " <<
                                 QString("%1").arg(posiciones > 1 ? double(ratio_medio / posiciones) : double(ratio_medio)) << " " <<
                                 QString::number(posiciones) << " " <<
                                 QString::number(cobertura_minima >= 500000000 ? 0 : cobertura_minima) << " " <<
                                 QString::number((posiciones > 1) ? cobertura_media / posiciones : cobertura_media) << " " <<
                                 QString::number(cobertura_maxima) << " " <<
                                 QString::number(sites_C) << " " <<
                                 QString::number(sites_nC) << " " <<
                                 QString::number(sites_mC) << " " <<
                                 QString::number(sites_hmC) << " " <<
                                 QString::number(distancia_minima >= 500000000 ? 0 : distancia_minima) << " " <<
                                 QString::number((posiciones > 1) ? distancia_media / posiciones : distancia_media) << " " <<
                                 QString::number(distancia_maxima) << "\n";
                        }
                        else
                        {
                            s << QString("%1").arg(double(dwt_valor)) << " " << //::number(double(dwt_valor), 'f', 3) << " " <<
                                 QString("%1").arg(posiciones > 1 ? double(ratio_medio / posiciones) : double(ratio_medio)) << " " <<
                                 "0 0 0 0 0 0 0 0 0 0 0\n";
                        }
                    }
                }
//               }

                s << "\n";
                //qDebug() << "fin escritura en fichero dmr: " << i;
            }
        }
        else
            s << "no DMRs were found\n";

        data.close();
        data_gff.close();
        qDebug() << "cerrando ficheros";
    }


    qDebug() << fichero;
}


// ************************************************************************************************
// ********************************ZONA MODIFICACION OPCIONES**************************************
// ************************************************************************************************
void HPG_Dhunter::on_mC_clicked()
{
    _mc = ui->mC->isChecked();

    if (!_mc && !_hmc)
    {
        ui->hmC->setChecked(true);
        _hmc = true;
    }


    ui->mC_cobertura->setEnabled(_mc);
    ui->hmC_cobertura->setEnabled(_hmc);
}

// ************************************************************************************************
void HPG_Dhunter::on_hmC_clicked()
{
    _hmc = ui->hmC->isChecked();

    if (!_mc && !_hmc)
    {
        ui->mC->setChecked(true);
        _mc = true;
    }

    ui->mC_cobertura->setEnabled(_mc);
    ui->hmC_cobertura->setEnabled(_hmc);
}

// ************************************************************************************************
void HPG_Dhunter::on_forward_clicked()
{
    _forward = ui->forward->isChecked();

    if (!_forward && !_reverse)
    {
        ui->reverse->setChecked(true);
        _reverse = true;
    }
}

// ************************************************************************************************
void HPG_Dhunter::on_reverse_clicked()
{
    _reverse = ui->reverse->isChecked();

    if (!_forward && !_reverse)
    {
        ui->forward->setChecked(true);
        _forward = true;
    }
}

void HPG_Dhunter::on_out_path_clicked()
{
    // abre ventana de explorador de directorios para seleccionar
    // el directorio donde se encuentran los cromosomas de un caso
    fichero = QFileDialog::getExistingDirectory( this,
                                              tr("Select an output folder"),
                                              (directorio) ? path : QDir::homePath() ,
                                              QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks
                                            );

    if(fichero.isEmpty() || fichero.isNull())
        fichero = "";
    else
    {
        ui->out_path_label->setText(fichero);
    }
}

// ************************************************************************************************
void HPG_Dhunter::on_all_chroms_toggled(bool checked)
{
    _all_chroms = checked;

    if (!_all_chroms)
    {
        ui->chromosomes_list->setEnabled(true);
    }
    else
    {
        ui->chromosomes_list->setEnabled(false);
        ui->chromosomes_list->clear();
    }
}

// ************************************************************************************************
void HPG_Dhunter::on_mC_cobertura_sliderMoved(int value)
{
    _mc_min_coverage = value;

    ui->mC_min_cov->setText(QString::number(value));
}

// ************************************************************************************************
void HPG_Dhunter::on_hmC_cobertura_sliderMoved(int value)
{
    _hmc_min_coverage = value;

    ui->hmC_min_cov->setText(QString::number(value));
}

// ************************************************************************************************
void HPG_Dhunter::on_threshold_valueChanged(int value)
{
    _threshold = float(value * 0.01);
    ui->threshold_label->setText(QString::number(double(_threshold), 'f', 2));
}

// ************************************************************************************************
void HPG_Dhunter::on_dmr_dwt_level_valueChanged(int value)
{
    _dmr_dwt_level = value;
}

// ************************************************************************************************
void HPG_Dhunter::on_stop_clicked()
{
    QEventLoop loop;

    foreach (Files_worker *h, files_worker)
        h->abort();

    ui->start->setEnabled(true);
    ui->start->setFocus();
    ui->stop->setEnabled(false);
    enabling_widgets(true);
    lista_chroms.clear();

    ui->progressBar->setValue(ui->progressBar->maximum());

    QTimer::singleShot(500, &loop, SLOT(quit()));
    loop.exec();

    qDebug() << "hilo de lectura de ficheros abortado";
}

// ************************************************************************************************
void HPG_Dhunter::on_genome_reference_currentIndexChanged(int index)
{
    switch (index)
    {
        case 0:
            ui->all_chroms->setEnabled(false);
            ui->selected_chrms->setChecked(true);
            break;
        case 1:
            ui->all_chroms->setEnabled(true);
            ui->all_chroms->setChecked(true);
            break;
        default:
            ui->all_chroms->setEnabled(false);
            ui->selected_chrms->setChecked(true);
    }
}

// ************************************************************************************************
void HPG_Dhunter::enabling_widgets (bool arg)
{
    ui->genome_reference->setEnabled(arg);

    if (ui->genome_reference->currentIndex() != 0 && arg)
        ui->all_chroms->setEnabled(arg);
    else
        ui->all_chroms->setEnabled(false);

    ui->selected_chrms->setEnabled(arg);
    ui->mC->setEnabled(arg);
    ui->hmC->setEnabled(arg);
    ui->forward->setEnabled(arg);
    ui->reverse->setEnabled(arg);
    ui->mC_cobertura->setEnabled(arg & _mc);
    ui->hmC_cobertura->setEnabled(arg & _hmc);
    ui->min_CpG_x_region->setEnabled(arg);
    ui->min_covSamples_x_region->setEnabled(arg);
    ui->threshold->setEnabled(arg);
    ui->dmr_dwt_level->setEnabled(arg);
}

// ************************************************************************************************
void HPG_Dhunter::on_mC_min_cov_textEdited(const QString &arg1)
{
    _mc_min_coverage = arg1.toInt();

    ui->mC_cobertura->setValue(arg1.toInt());
}

// ************************************************************************************************
void HPG_Dhunter::on_hmC_min_cov_textEdited(const QString &arg1)
{
    _hmc_min_coverage = arg1.toInt();

    ui->hmC_cobertura->setValue(arg1.toInt());
}
