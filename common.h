
#include <PartArray.h>
#include <Part.h>
#include <Vect.h>
#include <optional>

using namespace std;

static Vect sizePBC;

Vect radiusPBC(const Vect& a, const Vect& b);

double hamiltonianDipolarPBC(Part *a, Part *b);

void setPBCEnergies(PartArray & sys);

double calcEmax(PartArray &sys);