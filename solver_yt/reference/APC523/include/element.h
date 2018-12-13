#ifndef ELEMENT_H
#define ELEMENT_H

#include "common_definitions.h"

class Element
{	
	public:
		/* One element has four vertices. {vertices} is their coordinates. */
		matrix vertices;

		/* {id} is the identification number for each element.
		   The way we enumerate element can be shown in the following diagram ({ELE_NUMBER} = 4)
		
			|-----|-----|
			|  2  |  3  |
			|-----|-----|
			|  0  |  1  |
			|-----|-----|

		*/
		int id;

		/* {x_list} contains coordinates at each point where there is a degree of freedom */
		matrix x_list;

		/* Constructor */
		Element(int id_);

		/* Initialize {vertices} and {x_list} */
		void get_coordinates();		

};

/* {ElementFace} is inherited from {Element} */
class ElementFace: public Element
{
	public:
		/* {face_id} is the identification number for each face.
		   The way we enumerate face can be shown in the following diagram
	
			|---3---|
			|		|
			0		1
			|		|
			|---2---|

		*/
		int face_id;


		/* If the face has a neighbour face, {neigh_id} is the id for the the neighbour element. */
		int neigh_id;
		/* If the face has a neighbour face, {neigh_face_id} is the id for this this neighbour face. */
		int neigh_face_id;

		/* Whether the face is a boundary face*/
		bool is_boundary;

		/* The normal vector pointing outward */
		vec normal_vector;

		/* Constructor */		
    	ElementFace(int id_, int face_id_);		
};

#endif