from subprocess import STDOUT, check_output, TimeoutExpired, CalledProcessError
from os import listdir
from os.path import isfile, join
import subprocess

dataFolder="data/"
binFolder="bin/"
srcFolder="src/"

datasets=[4,8]

def check_output(resultFile, output) :
    resultat=""
    with open(resultFile) as f:
        for line in f:
            # print("output :", output)
            #print("line :", line)
            resultat+=line.replace('\n', '').replace(' ', '')
        # print("output :", output)
        # print("result :", resultat)
        return str(resultat==output)

def runDataset(name, i) :
    output=""
    try:
        binFile=binFolder+name
        if not isfile(binFile) :
            return "Error"
        output = subprocess.check_output(["mpirun", "--oversubscribe", "-np", "4", binFolder+name, dataFolder+"a_"+str(i), dataFolder+"b_"+str(i)], stderr=STDOUT, timeout=5,universal_newlines=True).strip().replace('\n', '').replace(' ', '')
        # print(output)
        return check_output(dataFolder+"/result_"+str(i), output)

    except  CalledProcessError as e :
        #we still have to check the output because
        #some student call returns non 0 value in case of success !

        output=e.output.strip().replace('\n', '')
        #print(output)
        if (e.returncode < 0) :
            return "Error"
        return check_output(dataFolder+"/result_"+str(i), output)
    except (TimeoutExpired) as e :
        return "Error"



def runProject(name) :
    result=""
    for d in datasets :
        result = result + ";" + runDataset(name, d)
    return result


def runAllProjects() :
    onlyfiles = [f for f in listdir(binFolder) if isfile(join(binFolder, f))]
    for f in onlyfiles:
        print(runProject(f))

def compileAndRunProjects() :
    onlyfiles = [f for f in listdir(srcFolder) if isfile(join(srcFolder, f))]
    for f in onlyfiles:
        name=binFolder + f[:-2]
        #print(name)
        try:
            output = subprocess.check_output(["mpicc","-std=c99", "-o" , name,  srcFolder+ f ,"-lm", "-fopenmp"], stderr=STDOUT, universal_newlines=True).strip().replace('\n', '')
            resultCompilation=(output=="")
        #print(output)
        except  CalledProcessError as e :
            resultCompilation=False
        print(f[:-2] + ";" + str(resultCompilation)  + runProject(f[:-2]))


compileAndRunProjects()
