/*************************************************************************
* Copyright (C) 2006 by luc Scholtes *
* luc.scholtes@hmg.inpg.fr *
* *
* This program is free software; it is licensed under the terms of the *
* GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

//Modifs : Parameters renamed as MeniscusParameters
//id1/id2 as id1 is the smallest grain, FIXME : wetting angle?
//FIXME : in triaxialStressController, change test about null force in updateStiffnessccc
//FIXME : needs "requestErase" somewhere

#include "Law2_ScGeom_CapillaryPhys_Capillarity1.hpp"
#include <yade/pkg/common/ElastMat.hpp>
#include <yade/pkg/dem/ScGeom.hpp>
#include <yade/pkg/dem/Ip2_FrictMat_FrictMat_CapillaryPhys1.hpp>
#include <yade/pkg/dem/Ip2_FrictMat_FrictMat_MindlinCapillaryPhys.hpp>
#include <yade/core/Omega.hpp>
#include <yade/core/Scene.hpp>
#include <yade/lib/base/Math.hpp>
// #include <yade/pkg/dem/CapillaryPhys1.hpp>

//#include "DelaunayInterpolation.hpp"

#include <iostream>
#include <fstream>

DT Law2_ScGeom_CapillaryPhys_Capillarity1::dtVbased;
DT Law2_ScGeom_CapillaryPhys_Capillarity1::dtPbased;

     Real Law2_ScGeom_CapillaryPhys_Capillarity1::intEnergy()
{
	Real energy=0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		CapillaryPhys1* phys = dynamic_cast<CapillaryPhys1*>(I->phys.get());
 		if(phys) {
		  Real liquidTension=0.073; //why to declare it twice? any other way to do it?
		  //BC: declare it like capillaryPressure in hpp, and make it 0.073 by default. Then it can be modified by the user.
			energy += liquidTension*phys->SInterface;}
// 			energy += 0.5*(phys->normalForce.squaredNorm()/phys->kn + phys->shearForce.squaredNorm()/phys->ks);}
 	}
	return energy;
} 

void Law2_ScGeom_CapillaryPhys_Capillarity1::triangulateData() {
    //test if R1>R2, Sinon changer
//   if (R1<R2)
//   {
//     v=R1;0
//     R1=R2;
//     R2=v;C
//   }

	if (solutions.size()>0) {
		LOG_WARN("Law2_ScGeom_CapillaryPhys_Capillarity1 asking triangulation for the second time. Ignored.");
		return;}
/// We get data from a file and input them in triangulations

    ifstream file (inputFilename.c_str());
    if (!file.is_open()) {
		LOG_ERROR("No data file found for capillary law. Check path and inputFilename.");
		return;}
//        ifstream file ("capillaryfile.txt");
       
    // convention R,v,d,s,e,f,p,a1,a2,dummy (just for the example, define your own,
    // dummy is because has too much values per line - with one useless extra colum,)
    MeniscusPhysicalData dat;
    double dummy;
    while ( file.good() ) {
//         file >>dat.R>>dat.volume>>dat.distance>>dat.surface>>dat.energy>>dat.force>>dat.succion>>dat.delta1>>dat.delta2>>dummy;
        file >>dat.succion>>dat.force>>dat.distance>>dat.volume>>dat.surface>>dat.arclength>>dat.delta1>>dat.delta2>>dat.R>>dummy;

        solutions.push_back(dat);
//         cout <<dat.R<<" "<<dat.volume<<" "<<dat.distance<<" "<<dat.surface<<" "<<dat.energy<<" "<<dat.force<<" "<<dat.succion<<" "<<dat.delta1<<" "<<dat.delta2<<endl;
//         cout <<dat.succion<<" "<<dat.force<<" "<<dat.distance<<" "<<dat.volume<<" "<<dat.surface<<" "<<dat.arcLength<<" "<<dat.delta1<<" "<<dat.delta2<<" "<<dat.R<<endl;

    }
    file.close();

    //We build two triangulations, one for imposed succion, the other for imposed volume
    for (unsigned int k=0; k<solutions.size(); k++) {
        DT::Vertex_handle vh = dtVbased.insert(K::Point_3(solutions[k].R, solutions[k].volume, solutions[k].distance));
        vh->id()=k;
        vh = dtPbased.insert(K::Point_3(solutions[k].R, solutions[k].succion, solutions[k].distance));
        vh->id()=k;
    }
// capillary = shared_ptr<capillarylaw>(new capillarylaw);
// capillary->fill("M(r=1)");
// capillary->fill("M(r=1.1)");
// capillary->fill("M(r=1.25)");
// capillary->fill("M(r=1.5)");
// capillary->fill("M(r=1.75)");
// capillary->fill("M(r=2)");
// capillary->fill("M(r=3)");
// capillary->fill("M(r=4)");
// capillary->fill("M(r=5)");
// capillary->fill("M(r=10)");
}

// int main()
// {
// // DT dtPbased;//here the coordinates are R,distance,P
// // DT dtVbased;//here they are R,distance,volume
// std::vector<MeniscusPhysicalData> solutions;
//
// /// We get data from a file and input them in triangulations
//
// ifstream file ("inputSolutions.dat");
// // convention R,v,d,s,e,f,p,a1,a2,dummy (just for the example, define your own,
// // dummy is because has too much values per line - with one useless extra colum,)
// MeniscusPhysicalData dat;
// double dummy;
// while ( file.good() ) {
// file >>dat.R>>dat.volume>>dat.distance>>dat.surface>>dat.energy>>dat.force>>dat.succion>>dat.delta1>>dat.delta2>>dummy;
// solutions.push_back(dat);
// cout <<dat.R<<" "<<dat.volume<<" "<<dat.distance<<" "<<dat.surface<<" "<<dat.energy<<" "<<dat.force<<" "<<dat.succion<<" "<<dat.delta1<<" "<<dat.delta2<<endl;
//
// }
// file.close();
//
// //We build two triangulations, one for imposed succion, the other for imposed volume
// for (int k=0; k<solutions.size(); k++) {
// DT::Vertex_handle vh = dtVbased.insert(K::Point_3(solutions[k].R, solutions[k].volume, solutions[k].distance));
// vh->id()=k;
// vh = dtPbased.insert(K::Point_3(solutions[k].R, solutions[k].succion, solutions[k].distance));
// vh->id()=k;
// }
//
// /// Now we will test a range of query points and write the results to a file
// /// The last points are out of range, hence should give zero.
// Meniscus m;
// m.data.distance=0.0; m.data.volume=0.0; m.data.R=0.0;
// ofstream ofile("out.dat");
//
// // for (double kk=0; kk<6; kk+=0.1) {
// // MeniscusPhysicalData res = interpolate(dtVbased,K::Point_3(m.data.R+kk, m.data.volume+kk, m.data.distance+kk), m, solutions);
// // MeniscusParameters solution(Pinterpol? capillary->Interpolate(R1,R2,Dinterpol, Pinterpol, currentIndexes) : MeniscusParameters());
// //
// //
// // ofile << res.succion<<" "<<func(m.data.R+kk,m.data.volume+kk,m.data.distance+kk)<<std::endl;
// // }
// // ofile.close();
//
// // /// A tricky way to perform complex tasks from a c++ program. We create a script for gnuplot and we use a system command to execute gnuplot on it.
// // ofstream script("plotscript");
// // script<< "plot './out.dat' using 0:1 title 'interpolated'\nreplot './out.dat' using 0:2 w l title 'analytical'\npause -1\n";
// // script.close();
// // return system( "gnuplot plotscript" ) ;
// }

YADE_PLUGIN((Law2_ScGeom_CapillaryPhys_Capillarity1));

using namespace std;


// void Law2_ScGeom_CapillaryPhys_Capillarity::postLoad(Law2_ScGeom_CapillaryPhys_Capillarity&){
//
// capillary = shared_ptr<capillarylaw>(new capillarylaw);
// capillary->fill("M(r=1)");
// capillary->fill("M(r=1.1)");
// capillary->fill("M(r=1.25)");
// capillary->fill("M(r=1.5)");
// capillary->fill("M(r=1.75)");
// capillary->fill("M(r=2)");
// capillary->fill("M(r=3)");
// capillary->fill("M(r=4)");
// capillary->fill("M(r=5)");
// capillary->fill("M(r=10)");
// }


// MeniscusParameters::MeniscusParameters()
// {
// V = 0;
// F = 0;
// delta1 = 0;
// delta2 = 0;
// };
//
// MeniscusParameters::MeniscusParameters(const MeniscusParameters &source)
// {
// V = source.V;
// F = source.F;
// delta1 = source.delta1;
// delta2 = source.delta2;
// }
//
// MeniscusParameters::~MeniscusParameters()
// {}

void Law2_ScGeom_CapillaryPhys_Capillarity1::action()
{
    if (!scene) cerr << "scene not defined!";
    shared_ptr<BodyContainer>& bodies = scene->bodies;
    if (dtPbased.number_of_vertices ()<1 ) triangulateData();
    if (fusionDetection && !bodiesMenisciiList.initialized) bodiesMenisciiList.prepare(scene);

    InteractionContainer::iterator ii = scene->interactions->begin();
    InteractionContainer::iterator iiEnd = scene->interactions->end();
    bool hertzInitialized = false;
    for (; ii!=iiEnd ; ++ii) {
/// interaction is real
        if ((*ii)->isReal()) {
            const shared_ptr<Interaction>& interaction = *ii;
            if (!hertzInitialized) {//NOTE: We are assuming that only one type is used in one simulation here
                if (CapillaryPhys1::getClassIndexStatic()==interaction->phys->getClassIndex()) hertzOn=false;
                else if (MindlinCapillaryPhys::getClassIndexStatic()==interaction->phys->getClassIndex()) hertzOn=true;
                else LOG_ERROR("The capillary law is not implemented for interactions using"<<interaction->phys->getClassName());
            }
            hertzInitialized = true;
            CapillaryPhys1* cundallContactPhysics=NULL;
            MindlinCapillaryPhys* mindlinContactPhysics=NULL;

/// contact physics depends on the contact law, that is used (either linear model or hertz model)
            if (!hertzOn) cundallContactPhysics = static_cast<CapillaryPhys1*>(interaction->phys.get());//use CapillaryPhys for linear model
            else mindlinContactPhysics = static_cast<MindlinCapillaryPhys*>(interaction->phys.get());//use MindlinCapillaryPhys for hertz model

            unsigned int id1 = interaction->getId1();
            unsigned int id2 = interaction->getId2();
            Body* b1 = (*bodies)[id1].get();
            Body* b2 = (*bodies)[id2].get();

/// interaction geometry search (this test is to compute capillarity only between spheres (probably a better way to do that)
            int geometryIndex1 = (*bodies)[id1]->shape->getClassIndex(); // !!!
            int geometryIndex2 = (*bodies)[id2]->shape->getClassIndex();
            if (!(geometryIndex1 == geometryIndex2)) continue;

/// definition of interacting objects (not necessarily in contact)
            ScGeom* currentContactGeometry = static_cast<ScGeom*>(interaction->geom.get());

/// Capillary components definition:
            Real liquidTension = 0.073; // superficial water tension at 20 Celsius degrees in N/m

/// Interacting Grains:
// If you want to define a ratio between YADE sphere size and real sphere size
            Real alpha=1;
            Real R1 = alpha*std::min(currentContactGeometry->radius2,currentContactGeometry->radius1) ;           
            Real R2 =alpha*std::max(currentContactGeometry->radius2,currentContactGeometry->radius1) ;

/// intergranular distance
            Real D = alpha*((b2->state->pos-b1->state->pos).norm()-(currentContactGeometry->radius1+ currentContactGeometry->radius2)); // scGeom->penetrationDepth could probably be used here?

            if ((currentContactGeometry->penetrationDepth>=0)|| D<=0 || createDistantMeniscii) { //||(scene->iter < 1) ) // a simplified way to define meniscii everywhere
                D=0; // defines fCap when spheres interpenetrate. D<0 leads to wrong interpolation has D<0 has no solution in the interpolation : this is not physically interpretable!! even if, interpenetration << grain radius.
                if (!hertzOn) {
                    if (fusionDetection && !cundallContactPhysics->meniscus) bodiesMenisciiList.insert((*ii));
                    cundallContactPhysics->meniscus=true;
                } else {
                    if (fusionDetection && !mindlinContactPhysics->meniscus) bodiesMenisciiList.insert((*ii));
                    mindlinContactPhysics->meniscus=true;
                }
            }
            Real Dinterpol = D/R1;

/// Suction (Capillary pressure):
            Real Pinterpol = 0;
//FIXME: why removing normalization?! (Bruno)
            if (!hertzOn) Pinterpol = cundallContactPhysics->isBroken ? 0 : capillaryPressure*R1/liquidTension;//??????//*(R2/liquidTension);
            else Pinterpol = mindlinContactPhysics->isBroken ? 0 : capillaryPressure*R1/liquidTension;//*(R2/liquidTension);
            if (!hertzOn) cundallContactPhysics->capillaryPressure = capillaryPressure;
            else mindlinContactPhysics->capillaryPressure = capillaryPressure;

/// Capillary solution finder:
            if ((Pinterpol>=0) && (hertzOn? mindlinContactPhysics->meniscus : cundallContactPhysics->meniscus)) {
//int* currentIndexes = hertzOn? mindlinContactPhysics->currentIndexes : cundallContactPhysics->currentIndexes;
//If P=0, we use null solution
//MeniscusParameters
// solution(Pinterpol? capillary->Interpolate(R1,R2,Dinterpol, Pinterpol, currentIndexes) : MeniscusParameters());
//FIXME: is it R1/R2 (less than 1) or R2/R1 (>1)?
                MeniscusPhysicalData solution = interpolate(dtPbased,K::Point_3(R2/R1, Pinterpol, Dinterpol), cundallContactPhysics->m, solutions);

/// capillary adhesion force
                Real Finterpol = solution.force;
                Vector3r fCap = Finterpol*R1*liquidTension*currentContactGeometry->normal;
                if (!hertzOn) cundallContactPhysics->fCap = fCap;
                else mindlinContactPhysics->fCap = fCap;
/// meniscus volume
//FIXME: hardcoding numerical constants is bad practice generaly, and it probably reveals a flaw in that case (Bruno)
                Real Vinterpol = solution.volume*pow(R1,3);
                Real SInterface = solution.surface*pow(R1,2);
                if (!hertzOn) { 
                    cundallContactPhysics->vMeniscus = Vinterpol;
                    cundallContactPhysics->SInterface = SInterface;
                    if (Vinterpol > 0) cundallContactPhysics->meniscus = true;
                    else cundallContactPhysics->meniscus = false;
                } else {
                    mindlinContactPhysics->vMeniscus = Vinterpol;
                    if (Vinterpol > 0) mindlinContactPhysics->meniscus = true;
                    else mindlinContactPhysics->meniscus = false;
                }
                if (!Vinterpol) {
                    if ((fusionDetection) || (hertzOn ? mindlinContactPhysics->isBroken : cundallContactPhysics->isBroken)) bodiesMenisciiList.remove((*ii));
                    if (D>0) scene->interactions->requestErase(interaction);
                }
///interfacial area

// Real SInterface = solution.surface;
// if (!hertzOn) {
// cundallContactPhysics->SInterface = SInterface;
// if (Vinterpol != 0) cundallContactPhysics->meniscus = true;
// else cundallContactPhysics->meniscus = false;
// } else {
// // mindlinContactPhysics->SInterface = SInterface;
// if (Vinterpol != 0) mindlinContactPhysics->meniscus = true;
// else mindlinContactPhysics->meniscus = false;
// }
// if (!Vinterpol) {
// if ((fusionDetection) || (hertzOn ? mindlinContactPhysics->isBroken : cundallContactPhysics->isBroken)) bodiesMenisciiList.remove((*ii));
// if (D>0) scene->interactions->requestErase(interaction);
// }
/// wetting angles
                if (!hertzOn) {
                    cundallContactPhysics->Delta1 = max(solution.delta1,solution.delta2);
                    cundallContactPhysics->Delta2 = min(solution.delta1,solution.delta2);
                } else {
                    mindlinContactPhysics->Delta1 = max(solution.delta1,solution.delta2);
                    mindlinContactPhysics->Delta2 = min(solution.delta1,solution.delta2);
                }
            }
///interaction is not real //If the interaction is not real, it should not be in the list
        } else if (fusionDetection) bodiesMenisciiList.remove((*ii));
    }
    if (fusionDetection) checkFusion();

    for (ii= scene->interactions->begin(); ii!=iiEnd ; ++ii) {
        if ((*ii)->isReal()) {
            CapillaryPhys1* cundallContactPhysics=NULL;
            MindlinCapillaryPhys* mindlinContactPhysics=NULL;
            if (!hertzOn) cundallContactPhysics = static_cast<CapillaryPhys1*>((*ii)->phys.get());//use CapillaryPhys for linear model
            else mindlinContactPhysics = static_cast<MindlinCapillaryPhys*>((*ii)->phys.get());//use MindlinCapillaryPhys for hertz model

            if ((hertzOn && mindlinContactPhysics->meniscus) || (!hertzOn && cundallContactPhysics->meniscus)) {
                if (fusionDetection) {//version with effect of fusion
//BINARY VERSION : if fusionNumber!=0 then no capillary force
                    short int& fusionNumber = hertzOn?mindlinContactPhysics->fusionNumber:cundallContactPhysics->fusionNumber;
                    if (binaryFusion) {
                        if (fusionNumber!=0) {	//cerr << "fusion" << endl;
                            hertzOn?mindlinContactPhysics->fCap:cundallContactPhysics->fCap = Vector3r::Zero();
                            continue;
                        }
                    }
//LINEAR VERSION : capillary force is divided by (fusionNumber + 1) - NOTE : any decreasing function of fusionNumber can be considered in fact
                    else if (fusionNumber !=0) hertzOn?mindlinContactPhysics->fCap:cundallContactPhysics->fCap /= (fusionNumber+1.);
                }
                scene->forces.addForce((*ii)->getId1(), hertzOn?mindlinContactPhysics->fCap:cundallContactPhysics->fCap);
                scene->forces.addForce((*ii)->getId2(),-(hertzOn?mindlinContactPhysics->fCap:cundallContactPhysics->fCap));
            }
        }
    }
//   Real Law2_ScGeom_CapillaryPhys_Capillarity1::capillaryEnergy()
// {
// 	Real energy=0;
// 	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
//
// 			energy += liquidTension*SInterface;
// 	}
// 	return energy;
// }
}

// capillarylaw::capillarylaw()
// {}

// void capillarylaw::fill(const char* filename)
// {
// data_complete.push_back(Tableau(filename));
//
// }

void Law2_ScGeom_CapillaryPhys_Capillarity1::checkFusion()
{
//Reset fusion numbers
    InteractionContainer::iterator ii = scene->interactions->begin();
    InteractionContainer::iterator iiEnd = scene->interactions->end();
    for( ; ii!=iiEnd ; ++ii ) {
        if ((*ii)->isReal()) {
            if (!hertzOn) static_cast<CapillaryPhys1*>((*ii)->phys.get())->fusionNumber=0;
            else static_cast<MindlinCapillaryPhys*>((*ii)->phys.get())->fusionNumber=0;
        }
    }

    list< shared_ptr<Interaction> >::iterator firstMeniscus, lastMeniscus, currentMeniscus;
    Real angle1 = -1.0;
    Real angle2 = -1.0;

    for ( int i=0; i< bodiesMenisciiList.size(); ++i ) { // i is the index (or id) of the body being tested
        CapillaryPhys1* cundallInteractionPhysics1=NULL;
        MindlinCapillaryPhys* mindlinInteractionPhysics1=NULL;
        CapillaryPhys1* cundallInteractionPhysics2=NULL;
        MindlinCapillaryPhys* mindlinInteractionPhysics2=NULL;
        if ( !bodiesMenisciiList[i].empty() ) {
            lastMeniscus = bodiesMenisciiList[i].end();
            for ( firstMeniscus=bodiesMenisciiList[i].begin(); firstMeniscus!=lastMeniscus; ++firstMeniscus ) { //FOR EACH MENISCUS ON THIS BODY...
                currentMeniscus = firstMeniscus;
                ++currentMeniscus;
                if (!hertzOn) {
                    cundallInteractionPhysics1 = YADE_CAST<CapillaryPhys1*>((*firstMeniscus)->phys.get());
                    if (i == (*firstMeniscus)->getId1()) angle1=cundallInteractionPhysics1->Delta1;//get angle of meniscus1 on body i
                    else angle1=cundallInteractionPhysics1->Delta2;
                }
                else {
                    mindlinInteractionPhysics1 = YADE_CAST<MindlinCapillaryPhys*>((*firstMeniscus)->phys.get());
                    if (i == (*firstMeniscus)->getId1()) angle1=mindlinInteractionPhysics1->Delta1;//get angle of meniscus1 on body i
                    else angle1=mindlinInteractionPhysics1->Delta2;
                }
                for ( ; currentMeniscus!= lastMeniscus; ++currentMeniscus) { //... CHECK FUSION WITH ALL OTHER MENISCII ON THE BODY
                    if (!hertzOn) {
                        cundallInteractionPhysics2 = YADE_CAST<CapillaryPhys1*>((*currentMeniscus)->phys.get());
                        if (i == (*currentMeniscus)->getId1()) angle2=cundallInteractionPhysics2->Delta1;//get angle of meniscus2 on body i
                        else angle2=cundallInteractionPhysics2->Delta2;
                    }
                    else {
                        mindlinInteractionPhysics2 = YADE_CAST<MindlinCapillaryPhys*>((*currentMeniscus)->phys.get());
                        if (i == (*currentMeniscus)->getId1()) angle2=mindlinInteractionPhysics2->Delta1;//get angle of meniscus2 on body i
                        else angle2=mindlinInteractionPhysics2->Delta2;
                    }
                    if (angle1==0 || angle2==0) cerr << "THIS SHOULD NOT HAPPEN!!"<< endl;

//cerr << "angle1 = " << angle1 << " | angle2 = " << angle2 << endl;

                    Vector3r normalFirstMeniscus = YADE_CAST<ScGeom*>((*firstMeniscus)->geom.get())->normal;
                    Vector3r normalCurrentMeniscus = YADE_CAST<ScGeom*>((*currentMeniscus)->geom.get())->normal;

                    Real normalDot = 0;
                    if ((*firstMeniscus)->getId1() == (*currentMeniscus)->getId1() || (*firstMeniscus)->getId2() == (*currentMeniscus)->getId2()) normalDot = normalFirstMeniscus.dot(normalCurrentMeniscus);
                    else normalDot = - (normalFirstMeniscus.dot(normalCurrentMeniscus));

                    Real normalAngle = 0;
                    if (normalDot >= 0 ) normalAngle = Mathr::FastInvCos0(normalDot);
                    else normalAngle = ((Mathr::PI) - Mathr::FastInvCos0(-(normalDot)));

                    if ((angle1+angle2)*Mathr::DEG_TO_RAD > normalAngle) {
                        if (!hertzOn) {
                            ++(cundallInteractionPhysics1->fusionNumber);    //count +1 if 2 meniscii are overlaping
                            ++(cundallInteractionPhysics2->fusionNumber);
                        }
                        else {
                            ++(mindlinInteractionPhysics1->fusionNumber);
                            ++(mindlinInteractionPhysics2->fusionNumber);
                        }
                    };
                }
            }
        }
    }
}
// Real Law2_ScGeom_CapillaryPhys_Capillarity1::intEnergy()
// {
// 	Real energy=0;
// 	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
// 		if(!I->isReal()) continue;
// // 		FrictPhys* phys = dynamic_cast<FrictPhys*>(I->phys.get());
// // 		if(phys) {
// 			energy += liquidTension*SInterface;
// // 			energy += 0.5*(phys->normalForce.squaredNorm()/phys->kn + phys->shearForce.squaredNorm()/phys->ks);}
//  	}
// 	return energy;
// }
// MeniscusParameters capillarylaw::Interpolate(Real R1, Real R2, Real D, Real P, int* index)
// { //cerr << "interpolate" << endl;
// if (R1 > R2) {
// Real R3 = R1;
// R1 = R2;
// R2 = R3;
// }
//
// Real R = R2/R1;
//cerr << "R = " << R << endl;
/*
MeniscusParameters result_inf;
MeniscusParameters result_sup;
MeniscusParameters result;
int i = 0;

for ( ; i < (NB_R_VALUES); i++)
{
Real data_R = data_complete[i].R;
//cerr << "i = " << i << endl;

if (data_R > R) // Attention a l'ordre ds lequel vont etre ranges les tableau R (croissant)
{
Tableau& tab_inf=data_complete[i-1];
Tableau& tab_sup=data_complete[i];

Real r=(R-tab_inf.R)/(tab_sup.R-tab_inf.R);

result_inf = tab_inf.Interpolate2(D,P,index[0], index[1]);
result_sup = tab_sup.Interpolate2(D,P,index[2], index[3]);

result.V = result_inf.V*(1-r) + r*result_sup.V;
result.F = result_inf.F*(1-r) + r*result_sup.F;
result.delta1 = result_inf.delta1*(1-r) + r*result_sup.delta1;
result.delta2 = result_inf.delta2*(1-r) + r*result_sup.delta2;

i=NB_R_VALUES;
//cerr << "i = " << i << endl;

}*/
// else if (data_complete[i].R == R)
// {
// result = data_complete[i].Interpolate2(D,P, index[0], index[1]);
// i=NB_R_VALUES;
// //cerr << "i = " << i << endl;
// }
// }
// return result;
// }
//
// Tableau::Tableau()
// {}

//Tableau::Tableau(const char* filename)

// {
// ifstream file (filename);
// file >> R;
// //cerr << "r = " << R << endl;
// int n_D; //number of D values
// file >> n_D;
//
// if (!file.is_open())
// {
// static bool first=true;
// if(first)
// {
// cout << "WARNING: cannot open files used for capillary law, all forces will be null. Instructions on how to download and install them is found here : https://yade-dem.org/wiki/CapillaryTriaxialTest." << endl;
// first=false;
// }
// return;
// }
// for (int i=0; i<n_D; i++)
// full_data.push_back(TableauD(file));
// file.close();
// //cerr << *this; // exemple d'utilisation de la fonction d'ecriture (this est le pointeur vers l'objet courant)
// }
//
// Tableau::~Tableau()
// {}

// MeniscusParameters Tableau::Interpolate2(Real D, Real P, int& index1, int& index2)
//
// { //cerr << "interpolate2" << endl;
// MeniscusParameters result;
// MeniscusParameters result_inf;
// MeniscusParameters result_sup;
//
// for ( unsigned int i=0; i < full_data.size(); ++i)
// {
// if (full_data[i].D > D ) // ok si D rang�s ds l'ordre croissant
//
// {
// Real rD = (D-full_data[i-1].D)/(full_data[i].D-full_data[i-1].D);
//
// result_inf = full_data[i-1].Interpolate3(P, index1);
// result_sup = full_data[i].Interpolate3(P, index2);
//
// result.V = result_inf.V*(1-rD) + rD*result_sup.V;
// result.F = result_inf.F*(1-rD) + rD*result_sup.F;
// result.delta1 = result_inf.delta1*(1-rD) + rD*result_sup.delta1;
// result.delta2 = result_inf.delta2*(1-rD) + rD*result_sup.delta2;
//
// i = full_data.size();
// }
// else if (full_data[i].D == D)
// {
// result=full_data[i].Interpolate3(P, index1);
//
// i=full_data.size();
// }
//
// }
// return result;
// }
//
// TableauD::TableauD()
// {}
//
// TableauD::TableauD(ifstream& file)
// {
// int i=0;
// Real x;
// int n_lines; //pb: n_lines is real!!!
// file >> n_lines;
// //cout << n_lines << endl;
//
// file.ignore(200, '\n'); // saute les caract�res (200 au maximum) jusque au caract�re \n (fin de ligne)*_
//
// if (n_lines!=0)
// for (; i<n_lines; ++i) {
// data.push_back(vector<Real> ());
// for (int j=0; j < 6; ++j) // [D,P,V,F,delta1,delta2]
// {
// file >> x;
// data[i].push_back(x);
// }
// }
// D = data[i-1][0];
// }
//
// MeniscusParameters TableauD::Interpolate3(Real P, int& index)
// { //cerr << "interpolate3" << endl;
// MeniscusParameters result;
// int dataSize = data.size();
//
// if (index < dataSize && index>0)
// {
// if (data[index][1] >= P && data[index-1][1] < P)
// {
// //compteur1+=1;
// Real Pinf=data[index-1][1];
// Real Finf=data[index-1][3];
// Real Vinf=data[index-1][2];
// Real Delta1inf=data[index-1][4];
// Real Delta2inf=data[index-1][5];
//
// Real Psup=data[index][1];
// Real Fsup=data[index][3];
// Real Vsup=data[index][2];
// Real Delta1sup=data[index][4];
// Real Delta2sup=data[index][5];
//
// result.V = Vinf+((Vsup-Vinf)/(Psup-Pinf))*(P-Pinf);
// result.F = Finf+((Fsup-Finf)/(Psup-Pinf))*(P-Pinf);
// result.delta1 = Delta1inf+((Delta1sup-Delta1inf)/(Psup-Pinf))*(P-Pinf);
// result.delta2 = Delta2inf+((Delta2sup-Delta2inf)/(Psup-Pinf))*(P-Pinf);
// return result;
//
// }
// }
// //compteur2+=1;
// for (int k=1; k < dataSize; ++k) // Length(data) ??
//
// { //cerr << "k = " << k << endl;
// if ( data[k][1] > P) // OK si P rangés ds l'ordre croissant
//
// { //cerr << "if" << endl;
// Real Pinf=data[k-1][1];
// Real Finf=data[k-1][3];
// Real Vinf=data[k-1][2];
// Real Delta1inf=data[k-1][4];
// Real Delta2inf=data[k-1][5];
//
// Real Psup=data[k][1];
// Real Fsup=data[k][3];
// Real Vsup=data[k][2];
// Real Delta1sup=data[k][4];
// Real Delta2sup=data[k][5];
//
// result.V = Vinf+((Vsup-Vinf)/(Psup-Pinf))*(P-Pinf);
// result.F = Finf+((Fsup-Finf)/(Psup-Pinf))*(P-Pinf);
// result.delta1 = Delta1inf+((Delta1sup-Delta1inf)/(Psup-Pinf))*(P-Pinf);
// result.delta2 = Delta2inf+((Delta2sup-Delta2inf)/(Psup-Pinf))*(P-Pinf);
// index = k;
//
// k=dataSize;
// }
// else if (data[k][1] == P)
//
// { //cerr << "elseif" << endl;
// result.V = data[k][2];
// result.F = data[k][3];
// result.delta1 = data[k][4];
// result.delta2 = data[k][5];
// index = k;
//
// k=dataSize;
// }
//
// }
// return result;
// }
//
// TableauD::~TableauD()
// {}
//
// std::ostream& operator<<(std::ostream& os, Tableau& T)
// {
// os << "Tableau : R=" << T.R << endl;
// for (unsigned int i=0; i<T.full_data.size(); i++) {
// os << "TableauD : D=" << T.full_data[i].D << endl;
// for (unsigned int j=0; j<T.full_data[i].data.size();j++) {
// for (unsigned int k=0; k<T.full_data[i].data[j].size(); k++)
// os << T.full_data[i].data[j][k] << " ";
// os << endl;
// }
// }
// os << endl;
// return os;
// }

BodiesMenisciiList1::BodiesMenisciiList1(Scene * scene)
{
    initialized=false;
    prepare(scene);
}

bool BodiesMenisciiList1::prepare(Scene * scene)
{
//cerr << "preparing bodiesInteractionsList" << endl;
    interactionsOnBody.clear();
    shared_ptr<BodyContainer>& bodies = scene->bodies;

    Body::id_t MaxId = -1;
    BodyContainer::iterator bi = bodies->begin();
    BodyContainer::iterator biEnd = bodies->end();
    for( ; bi!=biEnd ; ++bi )
    {
        MaxId=max(MaxId, (*bi)->getId());
    }
    interactionsOnBody.resize(MaxId+1);
    for ( unsigned int i=0; i<interactionsOnBody.size(); ++i )
    {
        interactionsOnBody[i].clear();
    }

    InteractionContainer::iterator ii = scene->interactions->begin();
    InteractionContainer::iterator iiEnd = scene->interactions->end();
    for( ; ii!=iiEnd ; ++ii ) {
        if ((*ii)->isReal()) {
            if (static_cast<CapillaryPhys1*>((*ii)->phys.get())->meniscus) insert(*ii);
        }
    }

    return initialized=true;
}

bool BodiesMenisciiList1::insert(const shared_ptr< Interaction >& interaction)
{
    interactionsOnBody[interaction->getId1()].push_back(interaction);
    interactionsOnBody[interaction->getId2()].push_back(interaction);
    return true;
}


bool BodiesMenisciiList1::remove(const shared_ptr< Interaction >& interaction)
{
    interactionsOnBody[interaction->getId1()].remove(interaction);
    interactionsOnBody[interaction->getId2()].remove(interaction);
    return true;
}

list< shared_ptr<Interaction> >& BodiesMenisciiList1::operator[] (int index)
{
    return interactionsOnBody[index];
}

int BodiesMenisciiList1::size()
{
    return interactionsOnBody.size();
}

void BodiesMenisciiList1::display()
{
    list< shared_ptr<Interaction> >::iterator firstMeniscus;
    list< shared_ptr<Interaction> >::iterator lastMeniscus;
    for ( unsigned int i=0; i<interactionsOnBody.size(); ++i )
    {
        if ( !interactionsOnBody[i].empty() )
        {
            lastMeniscus = interactionsOnBody[i].end();
//cerr << "size = "<<interactionsOnBody[i].size() << " empty="<<interactionsOnBody[i].empty() <<endl;
            for ( firstMeniscus=interactionsOnBody[i].begin(); firstMeniscus!=lastMeniscus; ++firstMeniscus )
            {
                if ( *firstMeniscus ) {
                    if ( firstMeniscus->get() )
                        cerr << "(" << ( *firstMeniscus )->getId1() << ", " << ( *firstMeniscus )->getId2() <<") ";
                    else cerr << "(void)";
                }
            }
            cerr << endl;
        }
        else cerr << "empty" << endl;
    }
}

BodiesMenisciiList1::BodiesMenisciiList1()
{
    initialized=false;
}

