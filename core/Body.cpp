
#include<yade/core/Body.hpp>
#include<limits>
#include<yade/core/Scene.hpp>
#include<yade/core/Omega.hpp>
#include<yade/core/InteractionContainer.hpp>

//! This could be -1 if id_t is re-typedef'ed as `int'
const Body::id_t Body::ID_NONE=Body::id_t(-1);

const shared_ptr<Body>& Body::byId(Body::id_t _id, Scene* rb){return (*((rb?rb:Omega::instance().getScene().get())->bodies))[_id];}
const shared_ptr<Body>& Body::byId(Body::id_t _id, shared_ptr<Scene> rb){return (*(rb->bodies))[_id];}

// return list of interactions of this particle
python::list Body::py_intrs(){
	python::list ret;
	for(Body::MapId2IntrT::iterator it=this->intrs.begin(),end=this->intrs.end(); it!=end; ++it) {  //Iterate over all bodie's interactions
		if(!(*it).second->isReal()) continue;
		ret.append((*it).second);
	}
	return ret;
}

