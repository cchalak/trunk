//
// C++ Interface: Law2_ScGeom_CapillaryPhys_Capillarity
/*************************************************************************
* Copyright (C) 2006 by luc Scholtes *
* luc.scholtes@hmg.inpg.fr *
* *
* This program is free software; it is licensed under the terms of the *
* GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include <yade/core/GlobalEngine.hpp>
#include <set>
#include <boost/tuple/tuple.hpp>

#include <vector>
#include <list>
#include <utility>
#include <yade/pkg/dem/CapillaryPhys1.hpp>
// #include "DelaunayInterpolation.hpp"

#include <string>
// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Delaunay_triangulation_3.h>
// #include <CGAL/Cartesian.h>
#include <iostream>
#include <fstream>

/**
This law allows one to take into account capillary forces/effects between spheres coming from the presence of interparticular liquid bridges (menisci).
refs:
- (french, lot of documentation) L. Scholtes, PhD thesis -> http://tel.archives-ouvertes.fr/tel-00363961/en/
- (english, less...) L. Scholtes et al. Micromechanics of granular materials with capillary effects. International Journal of Engineering Science 2009,(47)1, 64-75

The law needs ascii files M(r=i) with i=R1/R2 to work (downloaded from https://yade-dem.org/wiki/CapillaryTriaxialTest). They contain a set of results from the resolution of the Laplace-Young equation for different configurations of the interacting geometry and must be placed in the bin directory (where yade exec file is situated) to be taken into account.
The control parameter is the capillary pressure (or suction) Delta_u, defined as the difference between gas and liquid pressure: Delta_u = u_gas - u_liquid
Liquid bridges properties (volume V, extent over interacting grains delta1 and delta2) are computed as a result of Delta_u and the interacting geometry (spheres radii and interparticular distance)

Rk: - the formulation is valid only for pendular menisci involving two grains (pendular regime).
- an algorithm was developed by B. Chareyre to identify menisci overlaps on each spheres (menisci fusion).
- some assumptions can be made to reduce capillary forces when menisci overlap (binary->F_cap=0 if at least 1 overlap, linear->F_cap=F_cap/numberOfOverlaps)
*/

/// !!! This version is deprecated. It should be updated to the new formalism -> ToDo !!!

/// a class to store meniscus parameters -> Rk: is it really needed since CapillaryPhys exist?
// class MeniscusParameters1
// {
// public :
// Real V;
// Real F;
// Real S;
// Real L;
// Real delta1;
// Real delta2;
// int index1;
// int index2;
//
// MeniscusParameters1();
// MeniscusParameters1(const MeniscusParameters1 &source);
// ~MeniscusParameters1();
// };

/// R = ratio(RadiusParticle1 on RadiusParticle2). Here, 10 R values from interpolation files (yade/extra/capillaryFiles), R = 1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
//const int NB_R_VALUES = 10;

// class capillarylaw1; // fait appel a la classe def plus bas
class Interaction;


///This container class is used to check if meniscii overlap. Wet interactions are put in a series of lists, with one list per body.
class BodiesMenisciiList1
{
private:
vector< list< shared_ptr<Interaction> > > interactionsOnBody;

//shared_ptr<Interaction> empty;

public:
BodiesMenisciiList1();
BodiesMenisciiList1(Scene* body);
bool prepare(Scene* scene);
bool insert(const shared_ptr<Interaction>& interaction);
bool remove(const shared_ptr<Interaction>& interaction);
list< shared_ptr<Interaction> >& operator[] (int index);
int size();
void display();


bool initialized;
};

/// This is the constitutive law
class Law2_ScGeom_CapillaryPhys_Capillarity1 : public GlobalEngine
{
public :
    void checkFusion();
// shared_ptr<capillarylaw1> capillary;

    static DT dtVbased;
    static DT dtPbased;
    std::vector<MeniscusPhysicalData> solutions;


    BodiesMenisciiList1 bodiesMenisciiList;

    void action();
     Real intEnergy();
     Real waterVolume();
//     {return 0;}
//     void postLoad(Law2_ScGeom_CapillaryPhys_Capillarity1&);
    void triangulateData();

    YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(Law2_ScGeom_CapillaryPhys_Capillarity1,GlobalEngine,"This law allows one to take into account capillary forces/effects between spheres coming from the presence of interparticular liquid bridges (menisci).\n\nThe control parameter is the capillary pressure (or suction) Uc = ugas - Uliquid. Liquid bridges properties (volume V, extent over interacting grains delta1 and delta2) are computed as a result of the defined capillary pressure and of the interacting geometry (spheres radii and interparticular distance).\n\nReferences: in english [Scholtes2009b]_; more detailed, but in french [Scholtes2009d]_.\n\nThe law needs ascii files M(r=i) with i=R1/R2 to work (see https://yade-dem.org/index.php/CapillaryTriaxialTest). These ASCII files contain a set of results from the resolution of the Laplace-Young equation for different configurations of the interacting geometry."
    "\n\nIn order to allow capillary forces between distant spheres, it is necessary to enlarge the bounding boxes using :yref:`Bo1_Sphere_Aabb::aabbEnlargeFactor` and make the Ig2 define define distant interactions via:yref:`interactionDetectionFactor<Ig2_Sphere_Sphere_ScGeom::interactionDetectionFactor>`. It is also necessary to disable interactions removal by the constitutive law (:yref:`Law2<Law2_ScGeom_FrictPhys_CundallStrack::neverErase>=True`). The only combinations of laws supported are currently capillary law + :yref:`Law2_ScGeom_FrictPhys_CundallStrack` and capillary law + :yref:`Law2_ScGeom_MindlinPhys_Mindlin` (and the other variants of Hertz-Mindlin).\n\nSee CapillaryPhys-example.py for an example script.",
                   ((Real,capillaryPressure,0.,,"Value of the capillary pressure Uc defines as Uc=Ugas-Uliquid"))
                   ((bool,fusionDetection,false,,"If true potential menisci overlaps are checked"))
                   ((bool,initialized,false,," "))
                   ((bool,binaryFusion,true,,"If true, capillary forces are set to zero as soon as, at least, 1 overlap (menisci fusion) is detected"))
                   ((bool,hertzOn,false,,"|yupdate| true if hertz model is used"))
                   ((string,inputFilename,string("capillaryfile.txt"),,"the file with meniscus solutions, used for interpolation."))
                   ((bool,createDistantMeniscii,false,,"Generate meniscii between distant spheres? Else only maintain the existing one. For modeling a wetting path this flag should always be false. For a drying path it should be true for one step (initialization) then false, as in the logic of [Scholtes2009c]_"))
	          ,,
		  .def("intEnergy",&Law2_ScGeom_CapillaryPhys_Capillarity1::intEnergy,"define the energy of interfaces in unsaturated pendular state")
		  .def("waterVolume",&Law2_ScGeom_CapillaryPhys_Capillarity1::waterVolume,"return the total value of water in the sample")
		   
     );
};

// class TableauD
// {
// public:
// Real D;
// std::vector<std::vector<Real> > data;
// MeniscusParameters Interpolate3(Real P, int& index);
// TableauD();
// TableauD(std::ifstream& file);
// ~TableauD();
// };
//
// // Fonction d'ecriture de tableau, utilisee dans le constructeur pour test
// class Tableau;
// std::ostream& operator<<(std::ostream& os, Tableau& T);
//
// class Tableau
// {
// public:
// Real R;
// std::vector<TableauD> full_data;
// MeniscusParameters Interpolate2(Real D, Real P, int& index1, int& index2);
// std::ifstream& operator<< (std::ifstream& file);
// Tableau();
// Tableau(const char* filename);
// ~Tableau();
// };
//
// class capillarylaw
// {
// public:
// capillarylaw();
// std::vector<Tableau> data_complete;
// MeniscusParameters Interpolate(Real R1, Real R2, Real D, Real P, int* index);
// void fill (const char* filename);
// };

REGISTER_SERIALIZABLE(Law2_ScGeom_CapillaryPhys_Capillarity1);

/////////////////////////////////////////////////////////////////////
/*
TYPES:
Triangulation_vertex_base_with_id_3: we redefine a vertex base including an index for each vertex (available in CGAL in 2D but not in 3D),
MeniscusPhysicalData: the physical variables describing a capillary bridge, with a few algebraic operators
Meniscus: a structure combining MeniscusPhysicalData with some cached data allowing faster operations in multiple queries (pointer to the last cell found and its normals)

FUNCTIONS:
getIncidentVtxWeights: an interpolation algorithm which is 20x faster than the natural neighbor interpolation of CGAL.
returns a list of vertices incident to a query point and their respective weights
interpolate: uses the results of getIncidentVtxWeights combined with a data array to return a weighted average,
may be used with arbitrary data types provided they have the required algebraic operators
main: example usage
*/
/////////////////////////////////////////////////////////////////////

// namespace CGAL {
//
// helpful array for permutations
// int comb [6] = {1, 2, 3, 0, 1, 2};
//
// //Vertex base including an index for each vertex, adapted from CGAL::Triangulation_vertex_base_with_id_2
// template < typename GT, typename Vb = Triangulation_vertex_base_3<GT> >
// class Triangulation_vertex_base_with_id_3 : public Vb
// {
// int _id;
// public:
// typedef typename Vb::Cell_handle Cell_handle;
// typedef typename Vb::Point Point;
// template < typename TDS3 >
// struct Rebind_TDS {
// typedef typename Vb::template Rebind_TDS<TDS3>::Other Vb3;
// typedef Triangulation_vertex_base_with_id_3<GT, Vb3> Other;
// };
//
// Triangulation_vertex_base_with_id_3(): Vb() {}
// Triangulation_vertex_base_with_id_3(const Point & p): Vb(p) {}
// Triangulation_vertex_base_with_id_3(const Point & p, Cell_handle c): Vb(p, c) {}
// Triangulation_vertex_base_with_id_3(Cell_handle c): Vb(c) {}
// int id() const { return _id; }
// int& id() { return _id; }
// };
//
// The function returning vertices and their weights for an arbitrary point in R3 space.
// The returned triplet contains:
// - the output iterator pointing to the filled vector of vertices
// - the sum of weights for normalisation
// - a bool telling if the query point was located inside the convex hull of the data points
// template <class Dt, class OutputIterator>
// Triple< OutputIterator, // iterator with value type std::pair<Dt::Vertex_handle, Dt::Geom_traits::FT>
// typename Dt::Geom_traits::FT, // Should provide 0 and 1
// bool >
// getIncidentVtxWeights(const Dt& dt,
// const typename Dt::Geom_traits::Point_3& Q,
// OutputIterator nn_out, typename Dt::Geom_traits::FT & norm_coeff,
// std::vector<typename Dt::Geom_traits::Vector_3>& normals,
// typename Dt::Cell_handle& start = CGAL_TYPENAME_DEFAULT_ARG Dt::Cell_handle())
// {
// typedef typename Dt::Geom_traits Gt;
// typedef typename Gt::Point_3 Point;
// typedef typename Dt::Cell_handle Cell_handle;
// typedef typename Dt::Locate_type Locate_type;
// typedef typename Gt::FT Coord_type;
// CGAL_triangulation_precondition (dt.dimension()== 3);
// Locate_type lt; int li, lj;
// Cell_handle c = dt.locate( Q, lt, li, lj, start);
// bool updateNormals = (c!=start || normals.size()<4);
// if (updateNormals) normals.clear();
// if ( lt == Dt::VERTEX )
// {
// *nn_out++= std::make_pair(c->vertex(li),Coord_type(1));
// return make_triple(nn_out,norm_coeff=Coord_type(1),true);
// }
// else if (dt.is_infinite(c))
// return make_triple(nn_out, Coord_type(1), false);//point outside the convex-hull
// norm_coeff=0;
// for ( int k=0;k<4;k++ )
// {
// if (updateNormals) {
// normals.push_back(cross_product(c->vertex(comb[k])->point()-c->vertex(comb[k+1])->point(),c->vertex(comb[k])->point()-c->vertex(comb[k+2])->point()));
// normals[k] = normals[k]/
// ((c->vertex(k)->point()-c->vertex(comb[k])->point())*normals[k]);
// }
// Coord_type closeness = ((Q-c->vertex(comb[k])->point())*normals[k]);
// Coord_type w = closeness;
// *nn_out++= std::make_pair(c->vertex(k),w);
// norm_coeff += w;
// }
// start = c;
// return make_triple(nn_out,norm_coeff,true);
// }
//
//
//
// } //END NAMESPACE CGAL
//
// typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
// typedef CGAL::Simple_cartesian<double> K;
// typedef CGAL::Delaunay_triangulation_3<K>::Geom_traits Traits;
// typedef CGAL::Triangulation_vertex_base_with_id_3<Traits> Vb;
// typedef CGAL::Triangulation_cell_base_3<Traits> Cb;
// typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
// typedef CGAL::Delaunay_triangulation_3<Traits,Tds> DT;
// typedef std::vector< std::pair< DT::Vertex_handle, K::FT> > Vertex_weight_vector;
//
// template <class Dt, class DataOwner>
// typename DataOwner::Data interpolate (const Dt& dt, const typename Dt::Geom_traits::Point_3& Q, DataOwner& owner, const std::vector<typename DataOwner::Data>& rawData)
// {
// K::FT norm;
// Vertex_weight_vector coords;
// CGAL::Triple<std::back_insert_iterator<Vertex_weight_vector>,K::FT, bool> result = CGAL::getIncidentVtxWeights(dt, Q,std::back_inserter(coords), norm, owner.normals , owner.cell);
//
// typename DataOwner::Data data = typename DataOwner::Data();//initialize null solution
// if (!result.third) return data;// out of the convex hull, we return the null solution
// //else, we compute the weighted sum
// for (int k=0; k<coords.size(); k++) data += (rawData[coords[k].first->id()]*coords[k].second);
// return data*(1./result.second);
// }
//
//
// 