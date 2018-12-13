#ifndef CELL_H
#define CELL_H

#include "common_definitions.h"

class Cell
{	
	public:

		/* {id_} is the identification number for each control volume cell.
		   The way we enumerate cell can be shown in the following diagram ({CELL_NUMBER} = 4)
		
			|-----|-----|
			|  2  |  3  |
			|-----|-----|
			|  0  |  1  |
			|-----|-----|
	
		*/
		int id_;

    // To do: Build specific structures for things like "points" and "tensors"

		/* One cell has four vertices. {vertices} are their coordinates. */
		matrix vertices_;

    /* cell centor coordinates */
    vec cell_center_;

    /* Store face center coordinates */
    map<DIRECTION, vec> face_centers_;

		/* Store neighbouring cell ids. If the neighbour is already outside the boundary, assign {OUTSIDE_CELL_ID} */
    map<DIRECTION, int> neighbour_ids_;

		/* Constructor */
		Cell(int id);

		/* Initialize {vertices_} */
		void ini_coordinates();		

		/* Initialize {neighbour_ids_} */
		void ini_neighbours();	


};


#endif