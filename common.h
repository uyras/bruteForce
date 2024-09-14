
#include <PartArray.h>
#include <Part.h>
#include <Vect.h>
#include <optional>
#include <argumentum/argparse.h>

using namespace std;
using namespace argumentum;

extern Vect sizePBC;
extern bool PBC;
extern bool dbg;

extern string inFilename;
extern string outFilename;
extern int precision;
extern double irange;
extern double emin,emax;

Vect radiusPBC(const Vect& a, const Vect& b);

double hamiltonianDipolarPBC(Part *a, Part *b);

void setPBCEnergies(PartArray & sys);

double calcEmax(PartArray &sys);

void configParser(argument_parser & parser);