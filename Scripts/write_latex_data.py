import os

def write_latex_data(filename,argument,arg_value):

    if type(filename) is not str:
        raise ValueError("The first argument must be a string")

    if type(argument) is not str:
        raise ValueError("The second argument must be a string")

    if (type(arg_value) is not str):
        raise ValueError("The third argument must be a string")

    #If the file does not exist yet, I create it
    if os.path.isfile(filename) is False:
        a_file = open(filename, "w")
        line2write='%s = %s\n' % (argument,arg_value)
        a_file.write(line2write)
        a_file.close()
    #If the file does exist yet, I check if it contains already the line with the argument and, in case, I replace it,
    #otherwise I write it
    else:
        #get list of lines
        a_file = open(filename, "r")
        lines = a_file.readlines()
        a_file.close()
        a_file = open(filename, "w")
        tmp=0
        line=lines[0]
        for line in lines:
            if line.split(' = ')[0] != argument:
                a_file.write(line)
            else:
                line2write = '%s = %s\n' % (argument, arg_value)
                a_file.write(line2write)
                tmp=1

        if tmp==0:
            line2write = '%s = %s\n' % (argument, arg_value)
            a_file.write(line2write)

        a_file.close()


if __name__ == "__main__":
    from pathlib import Path
    home = str(Path.home())

    filename='%s/GIT/AC_Agulhas_eddy_2021/Scripts/prova.dat' % home
    argument = 'ciao4'
    arg_value=30.456

    arg_value = '%0.2f' % arg_value
    write_latex_data(filename,argument,arg_value)


