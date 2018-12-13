#include "cell.h"

Cell::Cell(int id)
{

  // To do: Pass flags to selectively initiate class variables
  // This is important to enhance performance!
  id_ = id;
  ini_coordinates();
  ini_neighbours();

  Delta_x_; //??
  Delta_y_;
  Nabla_x_;
  Nable_y_;
  delta_x_;
  delta_y_;
}

void Cell::ini_coordinates()
{
  int index_1 = id_%LINE_CELL_NUMBER;
  int index_2 = id_/LINE_CELL_NUMBER;
  double step_x = (X_2[0] - X_1[0])/LINE_CELL_NUMBER;
  double step_y = (X_2[1] - X_1[1])/LINE_CELL_NUMBER;

  vec x_lower_left{X_1[0]+step_x*index_1, X_1[1]+step_y*index_2};

  for (int i = 0; i < VERTICES_PER_LINE; ++i)
  {
    for (int j = 0; j < VERTICES_PER_LINE; ++j)
    {
      vec v;
      v.push_back(x_lower_left[0] + step_x*j);
      v.push_back(x_lower_left[1] + step_y*i);
      vertices_.push_back(v);
    }
  }
  cell_center_.push_back(x_lower_left[0] + step_x/2);
  cell_center_.push_back(x_lower_left[1] + step_y/2);


  // Rubbish code, needed to be changed
  vec v_west;
  v_west.push_back(cell_center_[0] - step_x/2);
  v_west.push_back(cell_center_[1]);  
  face_centers_[WEST] = v_west;

  vec v_east;
  v_east.push_back(cell_center_[0] + step_x/2);
  v_east.push_back(cell_center_[1]);  
  face_centers_[EAST] = v_east;

  vec v_south;
  v_south.push_back(cell_center_[0]);
  v_south.push_back(cell_center_[1] - step_y/2);  
  face_centers_[SOUTH] = v_south;

  vec v_north;
  v_north.push_back(cell_center_[0]);
  v_north.push_back(cell_center_[1] + step_y/2);  
  face_centers_[NORTH] = v_north;    
}

void Cell::ini_neighbours()
{

  // Western neighbouring cell
  if (id_%LINE_CELL_NUMBER == 0)
    neighbour_ids_[WEST] = OUTSIDE_CELL_ID; 
  else
    neighbour_ids_[WEST] = id_ - 1;

  // Eastern neighbouring cell
  if ((id_+1)%LINE_CELL_NUMBER == 0)
    neighbour_ids_[EAST] = OUTSIDE_CELL_ID;
  else
    neighbour_ids_[EAST] = id_ + 1;

  // Southern neighbouring cell
  if (id_/LINE_CELL_NUMBER == 0)
    neighbour_ids_[SOUTH] = OUTSIDE_CELL_ID;
  else
    neighbour_ids_[SOUTH] = id_ - LINE_CELL_NUMBER;

  // Northern neighbouring cell
  if (id_/LINE_CELL_NUMBER + 1 == LINE_CELL_NUMBER)
    neighbour_ids_[NORTH] = OUTSIDE_CELL_ID;
  else
    neighbour_ids_[NORTH] = id_ + LINE_CELL_NUMBER;


}

