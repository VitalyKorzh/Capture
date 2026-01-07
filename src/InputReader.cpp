#include "InputReader.h"
#include "StringReader.h"

#include <cmath>

InputReader::InputReader(std::istream &in)
{
    std::string line = "";
    work = true;
    error_message = "";
    numberLine = 0;
    
    {
        precision = 10;
        sigma = 0.;
        normaDensity = 1.;
        nParticles = 0;
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

    if (!generateInjectionLine())
    {
        work = false;
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
            }
            else if (line.find("ni") != std::string::npos && nz > 0)
            {
                ni.clear();
                ni.resize(nz);
                for (uint i = 0; i < nz; i++)
                {
                    double val;
                    in >> val;
                    if (val < 0)
                    {
                        errorMessage("не правильное значение плотности ионов [ni >= 0]");
                        return false;
                    }
                    ni[i] = val;
                }
                
                if (in.fail())
                {
                    errorMessage("не удалось прочитать ni");
                    return false;
                }

            }
                
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

bool InputReader::readPosition(std::istream &in, std::pair<double, double> &p)
{
    std::string line;
    
    double p1, p2;
    const uint N_PAR = 2;
    bool array[] = {false, false};

    for (uint i = 0; i < N_PAR; i++)
    {
        skip(in, line, true);
        arrayBit(array[0], StringReader::getDoubleParameter(line, "z ", p1));
        arrayBit(array[1], StringReader::getDoubleParameter(line, "r ", p2));
    }

    if (checkArray(array, N_PAR))
    {
        p.first = p1;
        p.second = p2; //потом добавить условие больше нуля
    }
    else
    {
        errorConfigConstNumberPar(
            "не указаны все параметры position [",
            {"z", "r"},
            array,
            N_PAR
        );
        return false;
    }

    if (in.fail()) 
    {
        errorMessage("не удалось прочитать position");
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
    sigma = -1.;
    theta = 0;
    position.first = 0;
    position.second = 0;
    std::string line;
    if (!getline(in, line, true, true))
        return false;

    while(line.find("count end") == std::string::npos)
    {
        if (!line.empty() && !isComment(line))
        {
            StringReader::getDoubleParameter(line, "sigma ", sigma);
            StringReader::getUnsignedParameter(line, "particles ", nParticles);
            StringReader::getDoubleParameter(line, "theta ", theta);

            if (line.find("position") != std::string::npos)
            {
                if (!readPosition(in, position))
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

    if (nParticles == 0)
    {
        errorMessage("указано не правильное число частиц particles[>0]");
        return false;
    }
    if (sigma < 0)
    {
        errorMessage("указано не правильная сечение sigma [>0]");
        return false;
    }
    if (theta < 0. || theta > 90.)
    {
        errorMessage("не указан правильный угол инжекции theta [>=0 <=90]");
        return false;
    }

    theta *= M_PI/180.;

    return true;
}


bool InputReader::generateInjectionLine()
{
    ns = 0;
    sArray.clear();
    index.clear();

    double cosTheta = cos(theta);
    double sinTheta = sin(theta);
    double z0 = position.first;
    double r0 = position.second;

    // z = z0 + t*cos(theta)
    // r = r0 + t*sin(theta)
    uint iz0 = 0;
    uint ir0 = 0;

    bool found = false;

    for (uint iz = 0; iz < nz; iz++)
    {
        double z1 = zArray[iz];
        double z2 = zArray[iz+1];

        if (z0 >= z1 && z0 < z2)
        {
            iz0 = iz;
            found = true;
            break;
        }
    }

    if (!found) {
        errorMessage("начальная точка не найдена по z");
        return false;
    }

    found = false;
    for (uint ir = 0; ir < nr; ir++)
    {
        double r1 = rArray[ir];
        double r2 = rArray[ir+1];

        if (r0 >= r1 && r0 < r2)
        {
            ir0 = ir;
            found = true;
            break;
        }

    }

    if (!found) {
        errorMessage("начальная точка не найдена по r");
        return false;
    }

    if (fabs(sinTheta) < 1e-10)
    {
        for (uint iz = 0; iz < nz; iz++) {
            index.emplace_back(iz, ir0);
            sArray.push_back(zArray[iz+1] - zArray[iz]);
            ns++;
        }
    }
    else if (fabs(cosTheta) < 1e-10)
    {
        for (uint ir = 0; ir < nr; ir++) {
            index.emplace_back(iz0, ir);
            sArray.push_back(rArray[ir+1] - rArray[ir]);
            ns++;
        }
    }
    else
    {
        int stepR = (sinTheta > 0) ? 1 : -1;
        int stepZ = (cosTheta > 0) ? 1 : -1;

        double tMaxR = (stepR > 0) ? (rArray[ir0+1]-r0)/sinTheta : (rArray[ir0]-r0)/sinTheta;
        double tMaxZ = (stepZ > 0) ? (zArray[iz0+1]-z0)/cosTheta : (zArray[iz0]-z0)/cosTheta;

        double tDeltaR = fabs((rArray[1] - rArray[0]) / sinTheta);
        double tDeltaZ = fabs((zArray[1] - zArray[0]) / cosTheta);
        
        int ir = ir0;
        int iz = iz0;
        
        // Проходим по всем ячейкам вдоль луча
        while (ir >= 0 && ir < (int)nr && iz >= 0 && iz < (int)nz) {
            // Добавляем текущую ячейку
            index.emplace_back(ir, iz);
            ns++;
            
            // Определяем, какую границу пересечем следующей
            if (tMaxR < tMaxZ) {
                // Пересекаем границу по R
                // Длина пути в текущей ячейке = расстояние до границы по R
                double segment_length = tMaxR - ((ns == 1) ? 0.0 : sArray.back());
                sArray.push_back(segment_length);
                
                // Переходим к следующей ячейке
                tMaxR += tDeltaR;
                ir += stepR;
                
                // Обновляем tDeltaR для новой ячейки (если сетка неравномерная)
                if (ir >= 0 && ir < (int)nr) {
                    int cellR = (stepR > 0) ? ir : ir + 1;
                    if (cellR >= 0 && cellR < (int)nr) {
                        tDeltaR = fabs((rArray[cellR + 1] - rArray[cellR]) / sinTheta);
                    }
                }
            } else {
                // Пересекаем границу по Z
                // Длина пути в текущей ячейке = расстояние до границы по Z
                double segment_length = tMaxZ - ((ns == 1) ? 0.0 : sArray.back());
                sArray.push_back(segment_length);
                
                // Переходим к следующей ячейке
                tMaxZ += tDeltaZ;
                iz += stepZ;
                
                // Обновляем tDeltaZ для новой ячейки
                if (iz >= 0 && iz < (int)nz) {
                    int cellZ = (stepZ > 0) ? iz : iz + 1;
                    if (cellZ >= 0 && cellZ < (int)nz) {
                        tDeltaZ = fabs((zArray[cellZ + 1] - zArray[cellZ]) / cosTheta);
                    }
                }
            }
            
            // Проверяем, не вышли ли за пределы сетки
            if (ir < 0 || ir >= (int)nr || iz < 0 || iz >= (int)nz) {
                break;
            }
        }
        
        // Для последней ячейки длина может быть неполной
        if (!sArray.empty() && sArray.size() < index.size()) {
            // Оцениваем оставшуюся длину в последней ячейке
            double last_length = std::min(tMaxR, tMaxZ) - sArray.back();
            if (last_length > 0) {
                sArray.push_back(last_length);
            }
        }

    }

    if (sArray.empty() || index.empty())
    {
        errorMessage("линия инжекции не пересекает сетку");
        return false;
    }

    return true;
}