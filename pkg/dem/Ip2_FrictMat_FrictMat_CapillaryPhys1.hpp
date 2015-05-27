/*************************************************************************
*  Copyright (C) 2007 by Bruno CHAREYRE                                  *
*  bruno.chareyre@hmg.inpg.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg/common/Dispatching.hpp>
#include<yade/pkg/common/ElastMat.hpp>


class Ip2_FrictMat_FrictMat_CapillaryPhys1 : public IPhysFunctor
{
	public :
		virtual void go(	const shared_ptr<Material>& b1,
					const shared_ptr<Material>& b2,
					const shared_ptr<Interaction>& interaction);

	FUNCTOR2D(FrictMat,FrictMat);
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(Ip2_FrictMat_FrictMat_CapillaryPhys1,IPhysFunctor, "RelationShips to use with Law2_ScGeom_CapillaryPhys_Capillarity1\n\n In these RelationShips all the interaction attributes are computed. \n\n.. warning::\n\tas in the others :yref:`Ip2 functors<IPhysFunctor>`, most of the attributes are computed only once, when the interaction is new.",
				       ((bool,computeDefault,true,,"bool to assign the default value of computeBridge.")),;
	  
	);
	
};
REGISTER_SERIALIZABLE(Ip2_FrictMat_FrictMat_CapillaryPhys1);



