#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#include <cstring>
#include <bitset>
#include <omp.h>
#include <PartArray.h>
#include <Part.h>
#include "common.h"

using namespace argumentum;
using namespace std;


unsigned int grayencode(unsigned int g) 
{
    return g ^ (g >> 1);
}

FILE * ofile;
int main(int argc, char* argv[])
{
    auto parser = argument_parser{};
    parser.config().program( argv[0] ).description( "Program for calculating DOS of the magnetic system with dipole-dipole "
        "interaction with OMP parallelism" );
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
	
    PartArray sys_;
    
    {
        ifstream f(inFilename);
        if (!f.is_open()) {
            cout<<"file read error: "<<inFilename<<endl;
            return 1;
        }
        f.close();
    }
    const char* inFilename_c = inFilename.c_str();
	if(!strcmp(inFilename_c + strlen(inFilename_c) - strlen(".mfsys"), ".mfsys"))
	{ 
		sys_.load(inFilename);
		sys_.state.hardReset();
	} else {
		string tmp;
		double a,b,c,d;
        ifstream f(inFilename);

		while (getline(f,tmp))
		{
			if (tmp.length()<4 || tmp[0] == '#')
				continue;
			else {
				istringstream(tmp)>>a>>b>>c>>d;
				Part p;
				p.pos.setXYZ(a,b,0);
				p.m.setXYZ(c,d,0);
				sys_.insert(p);
			}
		}
	}
    sys_.setInteractionRange(irange);
	if (PBC) setPBCEnergies(sys_);
	unsigned N = sys_.size();
	
    emax = calcEmax(sys_);
    emax += emax*0.01; //add 1% to the value
	emin = -emax;

 

    const double dividerGlob=pow(10,precision);
    const int memsize = (emax-emin)*dividerGlob+1;
    unsigned long long * dosGlob = new unsigned long long [memsize];
    for (unsigned i=0;i<memsize;++i){
        dosGlob[i]=0;
    }
	
    
    cout << "# file: "<<inFilename<<endl;
    cout << "#    N: "<<sys_.size()<<endl;
    cout << "#    E: " << sys_.E()<<endl;
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
	
	
    #pragma omp parallel
    {
        int threadCount = omp_get_num_threads(); 
        int thread = omp_get_thread_num();

        /////////// import the system
        ifstream f(inFilename);
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


        if (sys.size()>64){
            cout<<"System size is more than 64. It is impossible task! Exiting."<<endl;
        }


        /////////// reserve the memory
        
        const double divider=pow(10,precision);
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


        unsigned long long totalStates=1;
        totalStates = totalStates << (sys.size()-1);
        unsigned long long statesPerThread = totalStates/threadCount;

        unsigned long long stateFrom = statesPerThread*thread;
        unsigned long long stateTo = stateFrom+statesPerThread;
        if (thread==(threadCount-1)) stateTo = totalStates;

        if (dbg){
            #pragma omp for ordered schedule(static,1)
            for (int t=0; t<threadCount; ++t)
            {
                if (thread == 0) {
                    cerr<<"# 2^N="<<totalStates<<endl;
                    cerr<<"# Distribution of states by threads:"<<endl;
                    cerr<<"# thrd\tfrom\tto\tdelta"<<endl;
                }
                assert( thread==t );
                #pragma omp ordered
                {
                    cerr<<"# "<<thread<<"\t"<<stateFrom<<"\t"<<stateTo<<"\t"<<stateTo-stateFrom<<endl;
                }
            }
        }

        #pragma omp barrier
        
        unsigned long long binState;
        unsigned partnum = 0;

        //variables to write out percents
        unsigned long long every = (stateTo-stateFrom)/100;
        unsigned long long everyCounter = every;
        unsigned percent = 0;

        string stateStr = std::bitset< 64 >( grayencode(stateFrom) ).to_string().substr(64-sys.size());
        reverse(stateStr.begin(),stateStr.end());
        sys.state.fromString(stateStr);
        /*#pragma omp critical
        {
            cout<<thread<<" "<<sys.state.toString()<<" "<<
                stateStr<<" "<<
                stateFrom<<endl;
            cout.flush();
        }*/

        if (thread==0) cerr<<"# ";
        for (binState=stateFrom+1;binState<=stateTo;++binState){
            if (sys.E()<emin)
                dos[0]+=2;
            else if (sys.E()>emax)
                dos[memsize-1]+=2;
            else
                dos[int((sys.E()-emin)*divider)]+=2;
            partnum = ffsll(binState)-1;// ffsll() - find first set bit in word
            sys.parts[partnum]->rotate();
            --everyCounter;
            if (thread==0 && everyCounter==0){
                everyCounter = every;
                cerr<<percent<<"%..";
                cerr.flush();
                ++percent;
            }
        }
        if (thread==0) cerr<<"\n";

        //gather dos here
        #pragma omp critical
        for (unsigned i=0;i<memsize;++i){
            if (dos[i]>0){
                dosGlob[i] += dos[i];
            }
        }
        
        delete[] dos;
    }
	printf("\n");
	unsigned int chek = 0;

    /////////// print the dos
    for (unsigned i=0;i<memsize;++i){
        if (dosGlob[i]>0){
            // fprintf(ofile,"%f\t%llu\n",(i/dividerGlob)+emin,dosGlob[i]);
            fprintf(ofile,"%.*f\t%llu\n",precision,(i/dividerGlob)+emin,dosGlob[i]);
			chek += dosGlob[i];
        }
    }
	
	if (chek != pow(2, N))
		cout << "ERROR sum(dos) != 2^N" << endl;

    delete[] dosGlob;

}