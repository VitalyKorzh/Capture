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

bool commands(std::string &inputFile, std::string &outputFile, int argc, char **argv)
{
    
    bool command = true;

    for (int i = 1; i < argc; i++) 
    {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-H")
        {
        }
        else if (arg == "--clear" || arg == "-C")
        {

        }
        else if (arg == "-i" || arg == "--input")
        {
            const char *file;
            command = readFileName(i, argc, argv, file, nullptr, nullptr);
            if (command)
                inputFile = file;
        }
        else if (arg == "-o" || arg == "--out")
        {
            const char *file;
            command = readFileName(i, argc, argv, file, nullptr, nullptr);
            if (command)
                outputFile = file;
        }
        else
        {
            std::cout << "Команда не найдена, посмотреть возможные команды можно введя --help" << std::endl;
            command = false;
        }

    }

    return command;
}
