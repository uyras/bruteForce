#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
//#include <mpi.h>
#include "PartArray.h"
#include "Part.h"
#include <string>
#include <string.h>
#include <bitset>

#define L 24
#define N 24
#define RADIUS 100

void printhelp(FILE* f){
    const string text = 
        "Program for calculating DOS of the magnetic system with dipole-dipole\n" 
        "interaction with mpi-interconnect\n" 
        "\n"
        "Usage:\n"
        "mpirun -np <N> ./bruteForce <filename>\n"
        "        <E_min> <E_max> [<precision> [<range> [<output>]]]\n"
        "\n"
        "where mpirun might be any mpi execution program, even the srun\n"
        "\n"
        "avaliable params:\n"
        "    <N>         -   Number of MPI processes\n"
        "    <filename>  -   Path to text file with structure of the system. \n"
        "                    File format is described below.\n"
        "    <E_min>     -   Approximage minimal energy of the system. May be setted\n"
        "                    lower than actual value with a small margin.\n"
        "    <E_max>     -   Approx maximal energy level of the system. May be setted\n"
        "                    higher than actual value with a small margin.\n"
        "    <precision> -   Number of digits after decimal point for DOS storage.\n"
        "                    All numbers after this value will be rounded.\n"
        "                    Default is 0. The value may be below 0.\n"
        "    <range>     -   Maximal distance between particles which interaction\n"
        "                    counts in total energy. Partices located further are \n"
        "                    considered non-interacting. \n"
        "                    Default value is 0 which means all-to-all interaction.\n"
        "    <output>    -   Resulting file name where the program writes DOS.\n"
        "                    By default program prints DOS to STDOUT (to screen).\n"
        "\n"
        "Parameters <E_min> <E_max> and <precision> mostly affects on the RAM consumption.\n"
        "The total requred memory for each MPI process is calculated as^\n"
        "M=((<E_max>-<E_min>)*10^<precision>+1)*sizeof(long long).\n"
        "Typically sizeof(long) is 8 bytes.\n"
        "\n"
        "Input file format:\n"
        "It should be the text file, csv-like. Each dipole is defined in one line.\n"
        "Total lines number in file equals the number of dipoles in system.\n"
        "\n"
        "Fieds in single line:\n"
        "<X> <Y> <Mx> <MY>\n"
        "<X> and <Y> are coordinates of dipole,\n"
        "<Mx> and <My> are magnetic vector coordinates. This vector starts\n" 
        "at <X>,<Y> position, so it is relative.\n"
        "Feld delimiter is space or tab symbol.\n"
        "Empty lines and lines started from '#' are skipped.\n";
    fprintf(f,text.c_str());
}

int main(int argc, char* argv[])
{
    bool dbg=true;

    double irange=0;
    double emin,emax;
    int precision=0;
    FILE * ofile = stdout;

    /////////// read command-line params
    if (argc<4){
        printhelp(stdout);
        return 0;
    }

    ifstream f(argv[1]);
    if (!f.is_open()) {
        cout<<"file read error: "<<argv[1]<<endl;
        printhelp(stdout);
        return 0;
    }

    emin=atof(argv[2]);
    emax=atof(argv[3]);

    if (argc>=5){
        precision=atoi(argv[4]);
    }

    if (argc>=6){
        irange=atof(argv[5]);
    }

    if (argc>=7){
        ofile = fopen(argv[6],"w");
        if (ofile==NULL){
            cout<<"error open the file: "<<argv[6];
            return 0;
        }
    }


    /////////// import the system
    PartArray sys;
    sys.setInteractionRange(irange);

    if(!strcmp(argv[1] + strlen(argv[1]) - strlen(".mfsys"), ".mfsys"))
    {
        sys.load(string(argv[1]));
        sys.state.hardReset();
    } else {
        string tmp;
        double a,b,c,d;

        while (getline(f,tmp))
        {
            if (tmp.length()<4 || tmp[0] == '#')
                continue;
            else {
                istringstream(tmp)>>a>>b>>c>>d;
                Part p;
                p.pos.setXYZ(a,b,0);
                p.m.setXYZ(c,d,0);
                sys.insert(p);
            }
        }
    }
    f.close();


    if (dbg){
        cout<<"file '"<<argv[1]<<"': imported "<<sys.size()<<" parts"<<endl;
        cout<<"Energy: "<<sys.E()<<" range: All"<<endl;
    }

    if (sys.size()>64){
        cout<<"System size is more than 64. It is impossible task! Exiting."<<endl;
        return 0;
    }

    /////////// reserve the memory
    const double divider=pow(10,precision);
    const int memsize = (emax-emin)*divider+1;
    unsigned long long * dos = new unsigned long long [memsize];
    for (unsigned i=0;i<memsize;++i){
        dos[i]=0;
    }

    /////////// enumerate the states

    /**
     * gray code iteration algorithm (from wiki): 
     * https://en.wikipedia.org/wiki/Gray_code#Constructing_an_n-bit_Gray_code
     * To construct the binary-reflected Gray code iteratively, at step 0 start with the 
     * {\displaystyle \mathrm {code} _{0}=0}{\displaystyle \mathrm {code} _{0}=0}, 
     * and at step {\displaystyle i>0}i>0 find the bit position of the least significant 
     * 1 in the binary representation of {\displaystyle i}i and flip the bit at that 
     * position in the previous code 
     * {\displaystyle \mathrm {code} _{i-1}}\mathrm {code} _{i-1} 
     * to get the next code 
     * {\displaystyle \mathrm {code} _{i}}{\displaystyle \mathrm {code} _{i}}. 
     * The bit positions start 0, 1, 0, 2, 0, 1, 0, 3,...
     * */

    /*do{
        cout<<sys.state.toString()<<" | "<<sys.E()<<endl;
    } while(sys.state.next());*/

    unsigned long long binState;
    unsigned long long states=1;
    states = states << (sys.size()-1);
    unsigned partnum = 0;
    unsigned long long every = states/100;
    unsigned long long everyCounter = every;
    unsigned percent = 0;

    printf("# ");
    for (binState=1;binState<=states;++binState){
        if (sys.E()<emin)
            dos[0]+=2;
        else if (sys.E()>emax)
            dos[memsize-1]+=2;
        else
            dos[int((sys.E()-emin)*divider)]+=2;
        partnum = ffsll(binState)-1;// ffsll() - find first set bit in word
        sys.parts[partnum]->rotate();
        --everyCounter;
        if (everyCounter==0){
            everyCounter = every;
            printf("%d%%..",percent);
            fflush(stdout);
            ++percent;
        }
    }
    printf("\n");

    /////////// print the dos
    for (unsigned i=0;i<memsize;++i){
        if (dos[i]>0){
            fprintf(ofile,"%f\t%d\n",(i/divider)+emin,dos[i]);
        }
    }

    delete[] dos;

}