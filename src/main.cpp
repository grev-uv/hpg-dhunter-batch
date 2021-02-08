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


#include "hpg_dhunter.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    HPG_Dhunter w;

    w.setWindowTitle("DMR identification with DWT - vFDD");
//    w.setStyleSheet("QMainWindow {background: 'grey';}");
    w.setWindowIcon(QIcon(":/images/icon.png"));

    w.show();

    return a.exec();
}
