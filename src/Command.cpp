#include "Command.h"
#include <iostream>

bool readFileName(int &index, int argc, char** argv, const char*& file_name, const char *second_command, void (*secondCommand) ()) 
{
    if (++index < argc)
    {
        if (second_command)
        {
            if ((std::string) argv[index] == (std::string) second_command) 
            {
                if (secondCommand)
                    secondCommand();
                return false;   
            }
        }
        if (file_name)
            file_name = argv[index];
    }
    return true;
}


void help() 
{
    std::cout << "presicion=[=10 >0] - число знаков после запятой, выводимых в файл\n"
    << "normaN=[=1 >0.] - нормировка плотности {n*normaN см^-3}\n"
    << "mesh - сетка\n"
    << "\tz-axis - сетка по z\n"
    << "\t\t[1] array N - равномерный массив, N - число ячеек\n"
    << "\t\t\tmin - начальное значение\n" 
    << "\t\t\tmax - конечное значение\n"
    << "\t\t[2] n N - указать границы\n"
    << "\t\t\tvalues - значение границы\n"
    << "\tr-axis - сетка по r\n"
    << "\t\t[1] array N - равномерный массив, N - число ячеек\n"
    << "\t\t\tmin - начальное значение [=0 обязательно указать!]\n" 
    << "\t\t\tmax - конечное значение\n"
    << "\t\t[2] n N - указать границы\n"
    << "\t\t\tvalues - значение границы\n"
    << "\tNphi=[>=1] - число ячеек по углу равномерно разбвивает отрезок 0 до 2pi\n"
    << "\tni - плотность ионов от r"
    << "\t\tvalues - значение плотности в ячейки ir\n"
    << "mesh end\n"
    << "count - настройка подсчета профиля захвата\n"
    << "\tt-epsilon=[=0 >=0] - шаг по лучу, который можно считать одинаковым\n"
    << "\tt-epsilon-first=[=0 >=0] - минимальный шаг по лучу, чтобы определить следующию границу\n"
    << "\tBcenter=[=0 >=0] - магнитное поле в Гс для расчета ларморовского радиуса\n"
    << "\tinjector - параметры инжектора\n"
    << "\t\ttheta=[>0 <180] - угол инжекции (питч угол быстрых ионов)\n"
    << "\t\tsigma=[>0] - сечение захвата см2\n"
    << "\t\tparticles=[>0] - число запусков частиц\n"
    << "\t\tr0=[=0 >=0] - ширина гауссовского профиля пучка\n"
    << "\t\tdtheta=[=0 >=0] угловой разброс в градусах\n"
    << "\t\tr-injection=[-1] радиус входа луча, если <0 то по границе сетки\n"
    << "\t\t[plus-direction] - если установлен флаг то направляющий вектор (-sin(theta)*sin(phi), sin(theta)*cos(phi), cos(theta)), установлен по умолчанию\n"
    << "\t\t[minus-direction] - если установлен флаг то направляющий вектор (sin(theta)*sin(phi), -sin(theta)*cos(phi), cos(theta))\n"
    << "\t\tposition - точка через которую проходит луч\n"
    << "\t\t\trho=[>=0] - минимальное растояние между лучом и осью z (прицельный параметер)\n"
    << "\t\t\tphi=[>=0 <360] - азимутальный угол при минимальном сближеине\n"
    << "\t\t\tz=[] - положение по оси z при минимальном смещении\n"
    << "\t\tparticle - параметры быстрого иона\n"
    << "\t\t\tE=[>=0] - энергия в кэВ\n"
    << "\t\t\tZ=[>=1] - заряд иона в ед. элеметарного заряда\n"
    << "\t\t\tM=[>=1] - масса иона в массах протона\n"
    << "\tinjector end\n"
    << "count end\n";
}

bool commands(std::string &inputFile, std::string &outputFile, int argc, char **argv)
{
    
    bool command = true;

    for (int i = 1; i < argc; i++) 
    {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-H")
        {
            help();
            inputFile = "";
        }
        else if (arg == "-i" || arg == "--input")
        {
            const char *file = "";
            command = readFileName(i, argc, argv, file, nullptr, nullptr);
            if (command)
                inputFile = file;
        }
        else if (arg == "-o" || arg == "--out")
        {
            const char *file = "";
            command = readFileName(i, argc, argv, file, nullptr, nullptr);
            if (command)
                outputFile = file;
        }
        else
        {
            std::cout << "Команда не найдена" << std::endl;
            command = false;
        }

    }

    return command;
}
