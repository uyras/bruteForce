#include "common.h"

Vect sizePBC;
bool PBC;
bool dbg;

string inFilename;
string outFilename;
int precision;
double irange;
double emin,emax;

Vect radiusPBC(const Vect& a, const Vect& b){
    Vect dist = a - b;
    if (sizePBC.x && fabs(dist.x) > sizePBC.x/2){
        if (a.x < b.x){
            dist.x += sizePBC.x;
        } else {
            dist.x -= sizePBC.x;
        }
    }

    if (sizePBC.y && fabs(dist.y) > sizePBC.y/2){
        if (a.y < b.y){
            dist.y += sizePBC.y;
        } else {
            dist.y -= sizePBC.y;
        }
    }

    if (sizePBC.z && fabs(dist.z) > sizePBC.z/2){
        if (a.z < b.z){
            dist.z += sizePBC.z;
        } else {
            dist.z -= sizePBC.z;
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

double calcEmax(PartArray &sys)
{
	double Emax = 0;
	for (Part* elem: sys.parts) {
		for (Part* neigh: sys.neighbours.at(elem->Id()))
		{
			if (neigh != elem)
			{
				Emax += fabs(hamiltonianDipolarPBC(neigh, elem));
			}
		}
	}
	return Emax;
}

void configParser(argument_parser & parser)
{
    auto params = parser.params();

    params.add_parameter(inFilename,  "input").nargs(1).required().metavar("INFILE{.mfsys|.txt}")
        .help("Path to text file with structure of the system. File format is described below.");
    params.add_parameter(precision,"-p","--precision").maxargs(1).default_value(0).metavar("INT")
        .help("Number of digits after decimal point for DOS storage. "
            "All numbers after this value will be rounded. "
            "Default is 0. The value may be below 0.");
    params.add_parameter(irange,"-r","--range").maxargs(1).default_value(999999).metavar("DOUBLE")
        .help("Maximal distance between particles which interaction counts "
        "in total energy. Partices located further are considered "
        "non-interacting.");
    params.add_parameter(outFilename,"-o","--out").maxargs(1).default_value("").metavar("file.txt")
        .help("Resulting file name where the program writes DOS."
        "By default program prints DOS to STDOUT (to screen).");
    params.add_parameter(sizePBC.x,"-x").maxargs(1).default_value(0)
        .help("Size of the system along X axis for periodical boundary conditions"
        "If set 0 then PBC along X-axis is not applied. Default 0.");
    params.add_parameter(sizePBC.y,"-y").maxargs(1).default_value(0)
        .help("Size of the system along Y axis for periodical boundary conditions"
        "If set 0 then PBC along Y-axis is not applied. Default 0.");
    params.add_parameter(sizePBC.z,"-z").maxargs(1).default_value(0)
        .help("Size of the system along Z axis for periodical boundary conditions"
        "If set 0 then PBC along Z-axis is not applied. Default 0.");
    params.add_parameter(dbg,"-d","--debug")
        .help("enable debug mode");

    parser.config().epilog("Parameters <E_min> <E_max> and <precision> mostly affects on the RAM consumption. "
        "<E_min> <E_max> are calculating automatically"
        "The total requred memory for each openMP process is calculated as:\n"
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
        "<Mx> and <My> are magnetic vector coordinates. This vector starts " 
        "at <X>,<Y> position, so it is relative.\n"
        "Feld delimiter is space or tab symbol.\n"
        "Empty lines and lines started from '#' are skipped.");
    return;
}
