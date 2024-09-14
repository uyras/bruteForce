#include "common.h"

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