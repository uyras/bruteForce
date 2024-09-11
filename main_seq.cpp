#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#include <string.h>
#include <bitset>
#include "PartArray.h"
#include "Part.h"

double sizex = 816*4;
double sizey = 816*2;
double sizez = 816*2;
bool PBC = true;
double showPercentEvery=0.01; // 1.0 - every percent, 0.1 - ten times per percent

void printhelp(FILE* f){
    const string text = 
        "Program for calculating DOS of the magnetic system with dipole-dipole\n" 
        "interaction, non-paralle code, but with PBC\n" 
        "\n"
        "Usage:\n"
        "./bruteForce_seq <filename>\n"
        "        <E_min> <E_max> [<precision> [<range> [<output>]]]\n"
        "\n"
        "\n"
        "avaliable params:\n"
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
    fprintf(f,"%s",text.c_str());
}

Vect radiusPBC(const Vect& a, const Vect& b){
    Vect dist = a - b;
    if (fabs(dist.x) > sizex/2){
        if (a.x < b.x){
            dist.x += sizex;
        } else {
            dist.x -= sizex;
        }
    }

    if (fabs(dist.y) > sizey/2){
        if (a.y < b.y){
            dist.y += sizey;
        } else {
            dist.y -= sizey;
        }
    }

    if (fabs(dist.z) > sizez/2){
        if (a.z < b.z){
            dist.z += sizez;
        } else {
            dist.z -= sizez;
        }
    }
    
    return dist;
}

double hamiltonianDipolarPBC(Part *a, Part *b)
{
    Vect rij = radiusPBC(b->pos,a->pos);
    double r2, r, r5,E;
    r2 = rij.x * rij.x + rij.y * rij.y + rij.z * rij.z;
    r = sqrt(r2);
    r5 = r2 * r2 * r; //радиус в пятой
    
    E = //энергия считается векторным методом, так как она не нужна для каждой оси
            (( (a->m.x * b->m.x + a->m.y * b->m.y + a->m.z * b->m.z) * r2)
                -
                (3 * (b->m.x * rij.x + b->m.y * rij.y + b->m.z * rij.z) * (a->m.x * rij.x + a->m.y * rij.y + a->m.z * rij.z)  )) / r5;
    return E;
}

void setPBCEnergies(PartArray & sys)
{
    sys.E();
    // first update all neighbours
    sys.neighbours.clear();

    //определяем соседей частицы
    if (sys.interactionRange() != 0.){ //только если не все со всеми
        sys.neighbours.resize(sys.size());
        Part *part, *temp;
        for (unsigned i=0; i<sys.size(); i++){
            sys.neighbours[i].clear();
            part = sys[i];
            vector<Part*>::iterator iter = sys.parts.begin();
            while(iter!=sys.parts.end()){
                temp = *iter;
                if (temp != part && radiusPBC(part->pos,temp->pos).length() < sys.interactionRange()){
                    sys.neighbours[i].push_front(temp);
                }
                iter++;
            }
        }
    }
    sys.changeSystem();

    //then set the hamiltonian
    sys.setHamiltonian(hamiltonianDipolarPBC);
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

    sys.setInteractionRange(irange);
    if (PBC) setPBCEnergies(sys);


    if (dbg){
        cout<<"file '"<<argv[1]<<"': imported "<<sys.size()<<" parts"<<endl;
        cout<<"Energy: "<<sys.E()<<" range: "<<irange<<endl;
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
    unsigned long long every = states*showPercentEvery/100;
    unsigned long long everyCounter = every;
    unsigned percent = 0;

    double eOld = sys.E();

    printf("# ");
    int filenumber1 = 1;
    int filenumber2 = 1;
    int filenumber3 = 1;
    for (binState=1;binState<=states;++binState){

        if (false && eOld < -0.1){ //save states
            sys.save("gssr1_"+std::to_string(filenumber1)+".mfsys");
            printf("\ngs1 %d %.8f\n",filenumber1,eOld);
            ++filenumber1;
        }
/*
        if (true && eOld >= -0.1122 && eOld < -0.1121){
            sys.save("gs2_"+std::to_string(filenumber2)+".mfsys");
            printf("\ngs2 %d %.8f\n",filenumber2,eOld);
            ++filenumber2;
        }

        if (true && eOld >= -0.1121 && eOld < -0.110){
            sys.save("gs3_"+std::to_string(filenumber3)+".mfsys");
            printf("\ngs3 %d %.8f\n",filenumber3,eOld);
            ++filenumber3;
        }*/

        if (eOld<emin)
            dos[0]+=2;
        else if (eOld>emax)
            dos[memsize-1]+=2;
        else
            dos[int((eOld-emin)*divider)]+=2;
        partnum = ffsll(binState)-1;// ffsll() - find first set bit in word
        
        { // вращаем спин partnum
            unsigned j=0;
            for (Part* neigh : sys.neighbours[partnum]){
                if (neigh->state==sys.parts[partnum]->state) //assume it is rotated, inverse state in mind
                    eOld -=  2. * sys.eAt(partnum,j);
                else
                    eOld += 2. * sys.eAt(partnum,j);
                ++j;
            }
            sys.parts[partnum]->rotate();
        }
        --everyCounter;
        if (everyCounter==0){
            everyCounter = every;
            printf("%f%%..",percent*showPercentEvery);
            fflush(stdout);
            ++percent;
            eOld = sys.E();
        }
    }
    printf("\n");

    /////////// print the dos
    for (unsigned i=0;i<memsize;++i){
        if (dos[i]>0){
            fprintf(ofile,"%.*f\t%llu\n",precision,(i/divider)+emin,dos[i]);
        }
    }

    delete[] dos;

}