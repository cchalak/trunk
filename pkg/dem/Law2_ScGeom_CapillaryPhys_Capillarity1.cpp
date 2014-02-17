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

#include <yade/pkg/dem/Law2_ScGeom_CapillaryPhys_Capillarity1.hpp>
#include <yade/pkg/common/ElastMat.hpp>

#include <yade/pkg/dem/ScGeom.hpp>
#include <yade/pkg/dem/Ip2_FrictMat_FrictMat_CapillaryPhys1.hpp>
#include <yade/pkg/dem/Ip2_FrictMat_FrictMat_MindlinCapillaryPhys.hpp>
#include <yade/core/Omega.hpp>
#include <yade/core/Scene.hpp>
#include <yade/lib/base/Math.hpp>

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

     Real Law2_ScGeom_CapillaryPhys_Capillarity1::waterVolume()
{
	Real volume=0;
	FOREACH(const shared_ptr<Interaction>& I, *scene->interactions){
		if(!I->isReal()) continue;
		CapillaryPhys1* phys = dynamic_cast<CapillaryPhys1*>(I->phys.get());
 		if(phys) {
		  volume += phys->vMeniscus;}
 	}
	return volume;
} 

void Law2_ScGeom_CapillaryPhys_Capillarity1::triangulateData() {
    /// We get data from a file and input them in triangulations
    if (solutions.size()>0) {LOG_WARN("Law2_ScGeom_CapillaryPhys_Capillarity1 asking triangulation for the second time. Ignored."); return;}
    ifstream file (inputFilename.c_str());
    if (!file.is_open()) { LOG_ERROR("No data file found for capillary law. Check path and inputFilename."); return;}

    // convention R,v,d,s,e,f,p,a1,a2,dummy (just for the example, define your own,
    // dummy is because has too much values per line - with one useless extra colum,)
    MeniscusPhysicalData dat;
    double dummy;
    while ( file.good() ) {
        file >>dat.succion>>dat.force>>dat.distance>>dat.volume>>dat.surface>>dat.arcLength>>dat.delta1>>dat.delta2>>dat.R>>dummy;
        solutions.push_back(dat);
    }
    file.close();
    // Make lists of points with index, so we can use range insertion, more efficient
    // see http://doc.cgal.org/latest/Triangulation_3/index.html#Triangulation_3SettingInformationWhileInserting
    std::vector< std::pair<K::Point_3,unsigned> > pointsP, pointsV;
    for (unsigned int k=0; k<solutions.size(); k++) {
        pointsP.push_back(std::make_pair(K::Point_3(solutions[k].R, solutions[k].succion, solutions[k].distance),k));
        pointsV.push_back(std::make_pair(K::Point_3(solutions[k].R, solutions[k].volume, solutions[k].distance),k));
    }
    // and now range insertion
    dtPbased.insert(pointsP.begin(), pointsP.end());
    dtVbased.insert(pointsV.begin(), pointsV.end());
}


 
YADE_PLUGIN((Law2_ScGeom_CapillaryPhys_Capillarity1));

using namespace std;




void Law2_ScGeom_CapillaryPhys_Capillarity1::action()
{

    if (IsPressureImposed ==false) {
    InteractionContainer::iterator ii = scene->interactions->begin();
    InteractionContainer::iterator iiEnd = scene->interactions->end();  
    Real p0=capillaryPressure;
    Real pente;
    Real eps=0.01;
    GetVolumeforGivenSuction(p0);
    Real V0=waterVolume();

    Real p1=capillaryPressure+0.1;
    GetVolumeforGivenSuction(p1);
    Real V1=waterVolume();
    cout<<"V0="<< V0 << endl;
    cout<<"V1="<< V1 << endl;
    while (abs(Real (VolumeofWater-V1))>eps){

      pente= (p1-p0)/(V1-V0);
      p0=p1;
      V0=V1;
      p1=p1-pente*(V1-VolumeofWater);
      GetVolumeforGivenSuction(p1);
      V1=waterVolume();

      capillaryPressure=p1;
      cout<<"eps="<< (VolumeofWater-V1) << endl;
      cout<<"p0="<< p0 << endl;
      cout<<"V0="<< V0 << endl;
      cout<<"p1="<< p1 << endl;
      cout<<"V1="<< V1 << endl;
     }
   
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
}
    if (IsPressureImposed ==true) {
           
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
           
/// the parameter that takes into account the rugosity of the particles. 
            Real epsilon = 0;
/// Interacting Grains:
// If you want to define a ratio between YADE sphere size and real sphere size
            Real alpha=1;
            Real R1 = alpha*std::max(currentContactGeometry->radius2,currentContactGeometry->radius1);
            Real R2 =alpha*std::min(currentContactGeometry->radius2,currentContactGeometry->radius1);
            Real factor = std::max(R2/R1,1-R2/R1);
            R1 = R1-epsilon*factor;           
            R2 =R2-epsilon*(1-factor);

/// intergranular distance
            Real D = alpha*((b2->state->pos-b1->state->pos).norm()-(currentContactGeometry->radius1+ currentContactGeometry->radius2))+epsilon; // scGeom->penetrationDepth could probably be used here?

            if ((currentContactGeometry->penetrationDepth>=0)|| D<=0 || createDistantMeniscii) { //||(scene->iter < 1) ) // a simplified way to define meniscii everywhere
//                 D=0; // defines fCap when spheres interpenetrate. D<0 leads to wrong interpolation has D<0 has no solution in the interpolation : this is not physically interpretable!! even if, interpenetration << grain radius.
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
                MeniscusPhysicalData solution = interpolate1(dtPbased,K::Point_3(R2/R1, Pinterpol, Dinterpol), cundallContactPhysics->m, solutions);
//                 MeniscusPhysicalData solution = interpolate2(dtVbased,K::Point_3(R2/R1, Vinterpol, Dinterpol), cundallContactPhysics->m, solutions);
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

}

}
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

void Law2_ScGeom_CapillaryPhys_Capillarity1::GetVolumeforGivenSuction(Real suction)
{
     
    if (!scene) cerr << "scene not defined!";
    shared_ptr<BodyContainer>& bodies = scene->bodies;
    if (dtPbased.number_of_vertices ()<1 ) triangulateData();
    if (fusionDetection && !bodiesMenisciiList.initialized) bodiesMenisciiList.prepare(scene);
//     Real V=0;
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
           
/// the parameter that takes into account the rugosity of the particles. 
            Real epsilon = 0;
/// Interacting Grains:
// If you want to define a ratio between YADE sphere size and real sphere size
            Real alpha=1;
            Real R1 = alpha*std::max(currentContactGeometry->radius2,currentContactGeometry->radius1);
            Real R2 =alpha*std::min(currentContactGeometry->radius2,currentContactGeometry->radius1);
            Real factor = std::max(R2/R1,1-R2/R1);
            R1 = R1-epsilon*factor;           
            R2 =R2-epsilon*(1-factor);

/// intergranular distance
            Real D = alpha*((b2->state->pos-b1->state->pos).norm()-(currentContactGeometry->radius1+ currentContactGeometry->radius2))+epsilon; // scGeom->penetrationDepth could probably be used here?

            if ((currentContactGeometry->penetrationDepth>=0)|| D<=0 || createDistantMeniscii) { //||(scene->iter < 1) ) // a simplified way to define meniscii everywhere
//                 D=0; // defines fCap when spheres interpenetrate. D<0 leads to wrong interpolation has D<0 has no solution in the interpolation : this is not physically interpretable!! even if, interpenetration << grain radius.
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
            if (!hertzOn) Pinterpol = cundallContactPhysics->isBroken ? 0 : suction*R1/liquidTension;//??????//*(R2/liquidTension);
            else Pinterpol = mindlinContactPhysics->isBroken ? 0 : suction*R1/liquidTension;//*(R2/liquidTension);
            if (!hertzOn) cundallContactPhysics->capillaryPressure = suction;
            else mindlinContactPhysics->capillaryPressure = suction;

/// Capillary solution finder:
            if ((Pinterpol>=0) && (hertzOn? mindlinContactPhysics->meniscus : cundallContactPhysics->meniscus)) {
//int* currentIndexes = hertzOn? mindlinContactPhysics->currentIndexes : cundallContactPhysics->currentIndexes;
//If P=0, we use null solution
//MeniscusParameters
// solution(Pinterpol? capillary->Interpolate(R1,R2,Dinterpol, Pinterpol, currentIndexes) : MeniscusParameters());
//FIXME: is it R1/R2 (less than 1) or R2/R1 (>1)?
                MeniscusPhysicalData solution = interpolate1(dtPbased,K::Point_3(R2/R1, Pinterpol, Dinterpol), cundallContactPhysics->m, solutions);
//                 MeniscusPhysicalData solution = interpolate2(dtVbased,K::Point_3(R2/R1, Vinterpol, Dinterpol), cundallContactPhysics->m, solutions);
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
// 	V+= Vinterpol;
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
//     return V;

//     for (ii= scene->interactions->begin(); ii!=iiEnd ; ++ii) {
//         if ((*ii)->isReal()) {
//             CapillaryPhys1* cundallContactPhysics=NULL;
//             MindlinCapillaryPhys* mindlinContactPhysics=NULL;
//             if (!hertzOn) cundallContactPhysics = static_cast<CapillaryPhys1*>((*ii)->phys.get());//use CapillaryPhys for linear model
//             else mindlinContactPhysics = static_cast<MindlinCapillaryPhys*>((*ii)->phys.get());//use MindlinCapillaryPhys for hertz model
// 
//             if ((hertzOn && mindlinContactPhysics->meniscus) || (!hertzOn && cundallContactPhysics->meniscus)) {
//                 if (fusionDetection) {//version with effect of fusion
// //BINARY VERSION : if fusionNumber!=0 then no capillary force
//                     short int& fusionNumber = hertzOn?mindlinContactPhysics->fusionNumber:cundallContactPhysics->fusionNumber;
//                     if (binaryFusion) {
//                         if (fusionNumber!=0) {	//cerr << "fusion" << endl;
//                             hertzOn?mindlinContactPhysics->fCap:cundallContactPhysics->fCap = Vector3r::Zero();
//                             continue;
//                         }
//                     }
// //LINEAR VERSION : capillary force is divided by (fusionNumber + 1) - NOTE : any decreasing function of fusionNumber can be considered in fact
//                     else if (fusionNumber !=0) hertzOn?mindlinContactPhysics->fCap:cundallContactPhysics->fCap /= (fusionNumber+1.);
//                 }
//                 scene->forces.addForce((*ii)->getId1(), hertzOn?mindlinContactPhysics->fCap:cundallContactPhysics->fCap);
//                 scene->forces.addForce((*ii)->getId2(),-(hertzOn?mindlinContactPhysics->fCap:cundallContactPhysics->fCap));
//             }
//         }
//     }

}
