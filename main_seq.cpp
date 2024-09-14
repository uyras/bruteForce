#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#include <cstring>
#include <bitset>
#include <PartArray.h>
#include <Part.h>
#include "common.h"


double showPercentEvery=1; // 1.0 - every percent, 0.1 - ten times per percent

FILE * ofile;
int main(int argc, char* argv[])
{

    auto parser = argument_parser{};
    parser.config().program( argv[0] ).description( "Program for calculating DOS of the magnetic system with dipole-dipole "
        "interaction in sequential mode" );
    configParser(parser);

    /////////// read command-line params
    if ( !parser.parse_args( argc, argv, 1 ) )
      return 1;

    if (!outFilename.empty()){
        ofile = fopen(outFilename.c_str(),"w");
        if (ofile==NULL){
            cout<<"error open the file: "<<outFilename;
            return 1;
        }
    } else {
        ofile = stdout;
    }

    PBC = sizePBC.x || sizePBC.y || sizePBC.z;


    /////////// import the system
    PartArray sys;

    ifstream f(inFilename);
    if (!f.is_open()) {
        cout<<"file read error: "<<inFilename<<endl;
        return 1;
    }

    const char* inFilename_c = inFilename.c_str();
    if(!strcmp(inFilename_c + strlen(inFilename_c) - strlen(".mfsys"), ".mfsys"))
	{
        sys.load(inFilename);
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

    emax = calcEmax(sys);
    emax += emax*0.01; //add 1% to the value
	emin = -emax;

    if (sys.size()>64){
        cout<<"System size is more than 64. It is impossible task! Exiting."<<endl;
        return 0;
    }
	
	cout << "# file: " << inFilename<<endl;
    cout << "#    N: " << sys.size()<<endl;
    cout << "#    E: " << sys.E()<<endl;
	cout << "# emin: " << emin << endl;
	cout << "# emax: " << emax << endl;
    if (PBC){
        config::Instance()->set3D();
        cout<<"# PBC enabled; size: "<<sizePBC<<endl;
    } else {
        cout<<"# PBC disabled"<<endl;
    }
	cout << "# precision: " << precision << endl;	
	cout << "# range: " << irange << endl;

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
            cerr<<percent*showPercentEvery<<"%..";
            cerr.flush();
            ++percent;
            eOld = sys.E();
        }
    }
    cerr<<"\n";
    printf("\n");

	unsigned int chek = 0;
    /////////// print the dos
    for (unsigned i=0;i<memsize;++i){
        if (dos[i]>0){
            fprintf(ofile,"%.*f\t%llu\n",precision,(i/divider)+emin,dos[i]);
			chek += dos[i];
        }
    }

	if (chek != pow(2, sys.size()))
		cout << "ERROR sum(dos) != 2^N" << endl;

    delete[] dos;

}