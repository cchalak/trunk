/*************************************************************************
*  Copyright (C) 2006 by luc Scholtes                                    *
*  luc.scholtes@hmg.inpg.fr                                              *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#pragma once
#include<yade/pkg/dem/FrictPhys.hpp>
#include<yade/pkg/dem/DelaunayInterpolation.hpp>

class MeniscusPhysicalData {
public:
    double R;
    double volume;
    double distance;
    double surface;
    double energy;
    double force;
    double succion;
    double delta1;
    double delta2;
    double arcLength;
    //default ctor
    MeniscusPhysicalData() : R(0), volume(0), distance(0), surface(0), energy(0), force(0), succion(0), delta1(0), delta2(0), arcLength(0)  {}
    //ctor with input values
    MeniscusPhysicalData(const double& r, const double& v, const double& d, const double& s, const double& e, const double& f, const double& p, const double& a1, const double& a2, const double& arc) : R(r), volume(v), distance(d), surface(s), energy(e), force(f), succion(p), delta1(a1), delta2(a2), arcLength(arc){}

    //a minimal list of operators for the interpolation
    //these operators are requirements for the DataType template parameter of interpolate()
    //FIXME: algebra operations must include energy, perimeter, and any other new variable for including them in the interpolation
    MeniscusPhysicalData& operator+= (const MeniscusPhysicalData& m2) {
        R+=m2.R; volume+=m2.volume; distance+=m2.distance; surface+=m2.surface; energy+=m2.energy; force+=m2.force; succion+=m2.succion; delta1+=m2.delta1; delta2+=m2.delta2;
        return *this;}

    MeniscusPhysicalData operator* (const double& fact) const {
        return MeniscusPhysicalData(fact*R, fact*volume, fact*distance, fact*surface, fact*energy, fact*force, fact*succion, fact*delta1, fact*delta2, fact*arcLength);}

    const MeniscusPhysicalData& operator= (const MeniscusPhysicalData& m1) {
        R=m1.R; volume=m1.volume; distance=m1.distance; surface=m1.surface; energy=m1.energy; force=m1.force; succion=m1.succion; delta1=m1.delta1; delta2=m1.delta2; arcLength=m1.arcLength;
        return *this;}
};

//The structure for the meniscus: physical properties + cached values for fast interpolation
class Meniscus {
    public:
    typedef MeniscusPhysicalData Data;
        Data data;    //the variables of Laplace's problem
        DT::Cell_handle cell;        //pointer to the last location in the triangulation, for faster locate()
        std::vector<K::Vector_3> normals;// 4 normals relative to the current cell

        Meniscus() : data(), cell(DT::Cell_handle()), normals(std::vector<K::Vector_3>()) {}
};

class CapillaryPhys1 : public FrictPhys
{
	public :
		int currentIndexes [4]; // used for faster interpolation (stores previous positions in tables)
		Meniscus m;
		
		virtual ~CapillaryPhys1();

	YADE_CLASS_BASE_DOC_ATTRS_CTOR(CapillaryPhys1,FrictPhys,"Physics (of interaction) for Law2_ScGeom_CapillaryPhys_Capillarity.",
				 ((bool,meniscus,false,,"Presence of a meniscus if true"))
				 ((bool,isBroken,false,,"If true, capillary force is zero and liquid bridge is inactive."))
				 ((Real,capillaryPressure,0.,,"Value of the capillary pressure Uc defines as Ugas-Uliquid"))
				 ((Real,vMeniscus,0.,,"Volume of the menicus"))
				 ((Real,Delta1,0.,,"Defines the surface area wetted by the meniscus on the smallest grains of radius R1 (R1<R2)"))
				 ((Real,Delta2,0.,,"Defines the surface area wetted by the meniscus on the biggest grains of radius R2 (R1<R2)"))
				 ((Vector3r,fCap,Vector3r::Zero(),,"Capillary Force produces by the presence of the meniscus"))
				 ((Real,SInterface,0.,,"Fluid-Gaz Interfacial area"))
				 ((Real,arcLength,0.,,"Arc Length of the Fluid-Gaz Interface"))
				 ((short int,fusionNumber,0.,,"Indicates the number of meniscii that overlap with this one"))
				 ,createIndex();currentIndexes[0]=currentIndexes[1]=currentIndexes[2]=currentIndexes[3]=0;
				 );
	REGISTER_CLASS_INDEX(CapillaryPhys1,FrictPhys);
};
REGISTER_SERIALIZABLE(CapillaryPhys1);

// class EnergeticCapillaryPhys : public CapillaryPhys1
// {
// 	public :
// 		int currentIndexes [4]; // used for faster interpolation (stores previous positions in tables)
// 		
// 		virtual ~EnergeticCapillaryPhys();
// 
// 	YADE_CLASS_BASE_DOC_ATTRS_CTOR(EnergeticCapillaryPhys,CapillaryPhys1,"Physics (of interaction) for Law2_ScGeom_CapillaryPhys_Capillarity.",
// 				
// 				 ((Real,SInterface,0.,,"Fluid-Gaz Interfacial area"))
// 				 ((Real,arcLength,0.,,"Arc Length of the Fluid-Gaz Interface"))
// 				 ,createIndex();
// 				 );
// 	REGISTER_CLASS_INDEX(EnergeticCapillaryPhys,CapillaryPhys1);
// };
// REGISTER_SERIALIZABLE(EnergeticCapillaryPhys);
