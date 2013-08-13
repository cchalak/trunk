// 2010 © Václav Šmilauer <eudoxos@arcig.cz>

#include "Scene.hpp"
#include "Body.hpp"
#include "BodyContainer.hpp"
#include "Clump.hpp"
#ifdef YADE_OPENMP
	#include<omp.h>
#endif


CREATE_LOGGER(BodyContainer);
 
BodyContainer::BodyContainer(): lowestFree(0)
{}


BodyContainer::~BodyContainer(){}

void BodyContainer::clear(){
	body.clear(); lowestFree=0;
}

Body::id_t BodyContainer::findFreeId(){
	Body::id_t max=body.size();
	for(; lowestFree<max; lowestFree++){
		if(!(bool)body[lowestFree]) return lowestFree;
	}
	return body.size();
}

Body::id_t BodyContainer::insert(shared_ptr<Body>& b){
	Body::id_t newId=findFreeId();
	return insert(b,newId);
}

Body::id_t BodyContainer::insert(shared_ptr<Body>& b, Body::id_t id){
	assert(id>=0);
	if((size_t)id>=body.size()) body.resize(id+1);
	
	const shared_ptr<Scene>& scene=Omega::instance().getScene(); 
	b->iterBorn=scene->iter;
	b->timeBorn=scene->time;
	b->id=id;
	
	scene->doSort = true;
	
	body[id]=b;
	return id;
}

bool BodyContainer::erase(Body::id_t id){
	if(!exists(id)) return false;
	lowestFree=min(lowestFree,id);
	
	const shared_ptr<Body>& b=Body::byId(id);   //If the body is the last member of clump, the clump should be removed as well
	if ((b) and (b->isClumpMember())) {
		const shared_ptr<Body>& clumpBody=Body::byId(b->clumpId);
		const shared_ptr<Clump> clump=YADE_PTR_CAST<Clump>(clumpBody->shape);
		Clump::del(clumpBody, b);
		
		if (clump->members.size()==0) {            //Clump has no members any more. Remove it
			this->erase(b->clumpId);
		}
	}
	
	const shared_ptr<Scene>& scene=Omega::instance().getScene();
	for(Body::MapId2IntrT::iterator it=b->intrs.begin(),end=b->intrs.end(); it!=end; ++it) {  //Iterate over all bodie's interactions
		scene->interactions->requestErase((*it).second);
	}
	body[id]=shared_ptr<Body>();
	
	return true;
}

