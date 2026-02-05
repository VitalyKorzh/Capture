#include "InputReader.h"
#include "StringReader.h"

#include <cmath>

bool InputReader::readInjector(std::istream &in)
{
    double theta = -1;
    uint nParticles = 0;
    double r0 = 0;
    double sigma = -1;

    double rho = -1;
    double z = 0;
    double phi = -1;

    double deltaTheta = 0;

    uint Z = 1;
    uint M = 1;
    double E = 0;

    bool plusDirection = true;

    std::string line;
    if (!getline(in, line, true, true))
        return false;

    while(line.find("injector end") == std::string::npos)
    {
        StringReader::getDoubleParameter(line, "theta ", theta);
        StringReader::getDoubleParameter(line, "dtheta ", deltaTheta);
        StringReader::getUnsignedParameter(line, "particles ", nParticles);
        StringReader::getDoubleParameter(line, "r0 ", r0);
        StringReader::getDoubleParameter(line, "sigma ", sigma);

        if (line.find("minus-direction") != std::string::npos)
            plusDirection = false;

        if (line.find("plus-direction") != std::string::npos)
            plusDirection = true;

        if (line.find("position") != std::string::npos)
        {
            if (!readPosition(in, rho, z, phi))
                return false;
        }

        if (line.find("particle") != std::string::npos && line.find("particles") == std::string::npos)
        {
            if (!readParticle(in, E, Z, M))
                return false;
        }

        skip(in, line, true);
        if (in.fail()) 
        {
            errorMessage("не найдено закрытие injector end");
            return false;
        }
    }


    if (theta <= 0. || theta > 90.)
    {
        errorMessage("угол в неправильном диапозоне 0 < theta <= 90");
        return false;
    }

    if (deltaTheta < 0)
    {
        errorMessage("угловой разброс dtheta >= 0");
        return false;
    }

    if (nParticles == 0)
    {
        errorMessage("не правильное число частиц particles > 0");
        return false;
    }

    if (sigma <= 0.)
    {
        errorMessage("сечение захват > 0");
        return false;
    }

    if (r0 < 0.)
    {
        errorMessage("разброс пучка r0 >= 0");
        return false;
    }

    theta *= M_PI / 180.;

    injectors.push_back(Injector(rho, phi, z, sigma, r0, theta, deltaTheta, nParticles, E, M, Z, plusDirection));

    return true;
}

InputReader::InputReader(std::istream &in)
{
    std::string line = "";
    work = true;
    error_message = "";
    numberLine = 0;
    
    {
        precision = 10;
        normaDensity = 1.;
        t_epsilon = 0.;
        t_epsilon_first = 0.;
        Bcenter = 0.;
    }
    
    bool findMesh = false;
    bool findCount = false;

    while (getline(in, line, true))
    {
        if (line.empty() || isComment(line))
            continue;

        StringReader::getUnsignedParameter(line, "precision=", precision);
        StringReader::getDoubleParameter(line, "normaN=", normaDensity);

        if (normaDensity <= 0.)
            normaDensity = 1.;

        if (!findMesh && isLine(line, "mesh")) 
        {
            findMesh = true;
            work = readMesh(in);
            if (!work)
                return;
        }
        else if (isLine(line, "count")) {
            findCount = true;
            work = readCount(in);
            if (!work)
                return;
        }
    }

    if (!findMesh) 
    {
        work = false;
        errorMessage("не указан mesh");
    }
    if (!findCount)
    {
        work = false;
        errorMessage("не указан count");
    }

}

void InputReader::errorMessage(std::string error)
{
    error_message = "# " + std::to_string(numberLine) + ": " + error + "\n";
}

void InputReader::errorConfigConstNumberPar(std::string part1, const std::vector<std::string> PAR_NAMES, const bool *array, const uint N_STEP)
{
    std::string part2 = "";
    for (uint ik = 0; ik < N_STEP; ik++) 
    {
        if (!array[ik]) 
        {
            if (part2 != "")
                part2 += ", ";
            
            part2 = part2 + PAR_NAMES[ik];   
        }
    }
    part2 = part2 + "]";
    errorMessage(part1+part2);
}

uint InputReader::countSpace(std::string line) const
{
    uint space = 0;
    while (line[space] == ' ')
        space++;
    return space;
}

void InputReader::skip(std::istream &in, std::string &line, bool ignoreEqual)
{
    getline(in, line, true, ignoreEqual);
    while ((line == " " || line == "" || isComment(line)) && !in.fail())
    {
        getline(in, line, true, ignoreEqual);
    }

    if (in.fail())
        line = "";
}

bool InputReader::readMesh(std::istream &in)
{
    std::string line;
    if (!getline(in, line, true))
        return false;

    nphi = 0;

    while(line.find("mesh end") == std::string::npos)
    {
        if (!line.empty() && !isComment(line))
        {
            if (line.find("z-axis") != std::string::npos) 
            {
                if (!readAxis(in, zArray, nz, "z")) //добавить на условие больше нуля
                return false;
            }
            else if (line.find("r-axis") != std::string::npos)
            {
                if (!readAxis(in, rArray, nr, "r"))
                    return false;

                if (rArray.front() != 0.)
                {
                    errorMessage("сетка по r должна начинаться с нуля");
                    return false;
                }

            }
            else if (line.find("ni") != std::string::npos && nz > 0)
            {
                ni.clear();
                ni.resize(nz*nr);
                for (uint i = 0; i < nr; i++)
                {
                    double val;
                    in >> val;
                    if (val < 0)
                    {
                        errorMessage("не правильное значение плотности ионов [ni >= 0]");
                        return false;
                    }

                    for (uint j = 0; j < nz; j++)
                        ni[i+nr*j] = val;
                }
                
                if (in.fail())
                {
                    errorMessage("не удалось прочитать ni");
                    return false;
                }

            }

            StringReader::getUnsignedParameter(line, "Nphi=", nphi);
                
            if(in.fail()) {
                errorMessage("ошибки при чтения данных");
                return false;
            }
        }
        skip(in, line);
        if (in.fail()) 
        {
            errorMessage("не найдено закрытие mesh end");
            return false;
        }
    }

    phiArray.clear();
    phiArray.resize(nphi+1);
    for (uint iphi = 0; iphi < nphi+1; iphi++)
        phiArray[iphi] = 2.*M_PI * iphi / nphi;

    if (nphi == 0)
    {
        errorMessage("Nphi не указан");
        return false;
    }

    if (zArray.empty())
    {
        errorMessage("z-axis не указан");
        return false;
    }
    if (ni.empty())
    {
        errorMessage("ni не задан");
        return false;
    }
    if (rArray.empty())
    {
        errorMessage("r-axis не указан");
        return false;
    }
    
    return true;
}

bool InputReader::readPosition(std::istream &in, double &rho, double &z, double &phi)
{
    std::string line;

    const uint N_PAR = 3;
    bool array[] = {false, false, false};
    z = 0;

    for (uint i = 0; i < N_PAR; i++)
    {
        skip(in, line, true);
        arrayBit(array[0], StringReader::getDoubleParameter(line, "z ", z));
        arrayBit(array[1], StringReader::getDoubleParameter(line, "rho ", rho));
        arrayBit(array[2], StringReader::getDoubleParameter(line, "phi ", phi));
    }

    if (rho < 0)
    {
        errorMessage("не правельно указано минимальное расстояние до оси rho >= 0");
        return false;
    }
    if (phi < 0 || phi > 360)
    {
        errorMessage("не правильный угол 0<=phi<360");
        return false;
    }

    if (!checkArray(array, N_PAR))
    {
        errorConfigConstNumberPar(
            "не указаны все параметры position [",
            {"z", "rho", "phi"},
            array,
            N_PAR
        );
        return false;
    }

    phi *= M_PI/180.;

    if (in.fail()) 
    {
        errorMessage("не удалось прочитать position");
        return false;
    }

    return true;
}


bool InputReader::readParticle(std::istream &in, double &E, uint &Z, uint &M)
{
    std::string line;

    const uint N_PAR = 3;
    bool array[] = {false, false, false};

    for (uint i = 0; i < N_PAR; i++)
    {
        skip(in, line, true);
        arrayBit(array[0], StringReader::getDoubleParameter(line, "E ", E));
        arrayBit(array[1], StringReader::getUnsignedParameter(line, "Z ", Z));
        arrayBit(array[2], StringReader::getUnsignedParameter(line, "M ", M));
    }

    if (E < 0)
    {
        errorMessage("не правельно указана энергия E >= 0");
        return false;
    }
    if (Z == 0)
    {
        errorMessage("не правильный заряд Z > 0");
        return false;
    }
    if (M == 0)
    {
        errorMessage("не правильная масса M > 0");
        return false;
    }

    if (!checkArray(array, N_PAR))
    {
        errorConfigConstNumberPar(
            "не указаны все параметры position [",
            {"E", "Z", "M"},
            array,
            N_PAR
        );
        return false;
    }

    if (in.fail()) 
    {
        errorMessage("не удалось прочитать particle");
        return false;
    }

    return true;
}


bool InputReader::readAxis(std::istream &in, darray &axis, uint &size, const std::string &name)
{
    std::string line;
    skip(in, line, true);
    if (StringReader::getUnsignedParameter(line, "n ", size)) 
    {
        if (size == 0) 
        {
            errorMessage("число разбиений должно n[>=1]");
            return false;
        }
        axis.clear();
        axis.reserve(size+1);
        for (uint i = 0; i < size+1; i++) 
        {
            double val;
            in >> val;
            if (i != 0 && val <= axis.back()) 
            {
                errorMessage("сетка по " + name + " задается по возрастанию");
                return false;
            }
            axis.push_back(val);
        }
    }
    else if (StringReader::getUnsignedParameter(line, "array ", size))
    {
        if (size == 0) 
        {
            errorMessage("число разбиений должно n[>=1]");
            return false;
        }
        axis.clear();
        axis.reserve(size+1);

        double min = 0;
        double max = 0;

        const uint N_PAR = 2;
        bool array[] = {false, false};

        for (uint i = 0; i < N_PAR; i++)
        {
            skip(in, line, true);
            arrayBit(array[0], StringReader::getDoubleParameter(line, "min ", min));
            arrayBit(array[1], StringReader::getDoubleParameter(line, "max ", max));
        }

        if (checkArray(array, N_PAR))
        {
            if (max <= min)
            {
                errorMessage("не правельные границы min < max");
                return false;
            }

            for (uint i = 0; i < size+1; i++)
                axis.push_back(min + (max-min)/size*i);

        }
        else {
            errorConfigConstNumberPar("указаны не все праметры [", {"min", "max"}, array, N_PAR);
            return  false;
        }
    }
    else{
        errorMessage("не известен способ разбиения интервала");
        return false;
    }


    if (in.fail())
    {
        errorMessage("не удалось прочитать разбиение по " + name);
        return false;
    }

    return true;
}

bool InputReader::readCount(std::istream &in)
{
    std::string line;
    if (!getline(in, line, true, true))
        return false;

    while(line.find("count end") == std::string::npos)
    {
        if (!line.empty() && !isComment(line))
        {

            StringReader::getDoubleParameter(line, "t-epsilon ", t_epsilon);
            StringReader::getDoubleParameter(line, "t-epsilon-first", t_epsilon_first);
            StringReader::getDoubleParameter(line, "Bcenter ", Bcenter);

            if (isLine(line, "injector"))
            {
                if (!readInjector(in))
                    return false;
            }

            if(in.fail()) {
                errorMessage("ошибки при чтения данных");
                return false;
            }
        }
        skip(in, line, true);
        if (in.fail()) 
        {
            errorMessage("не найдено закрытие count end");
            return false;
        }
    }

    if (t_epsilon < 0)
    {
        errorMessage("указан не правильный порог t-epsilon >= 0");
        return false;
    }

    if (Bcenter < 0)
    {
        errorMessage("указано не верное магнитное поле Bcenter>=0");
        return false;
    }

    return true;
}