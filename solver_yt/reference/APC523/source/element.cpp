 #include "element.h"


//**************************************************************************************
// {public} functions of {Element} start from here
Element::Element(int id_)
{
    id = id_;
    get_coordinates();

}

void Element::get_coordinates()
{
    int index_1 = id%LINE_DIVISION;
    int index_2 = id/LINE_DIVISION;
    double step_x = (X_2[0] - X_1[0])/LINE_DIVISION;
    double step_y = (X_2[1] - X_1[1])/LINE_DIVISION;
    double sub_step_x = step_x/ORDER;
    double sub_step_y = step_y/ORDER;

    double x_lower_left[DIM] = {X_1[0]+step_x*index_1, X_1[1]+step_y*index_2};

    for (int i = 0; i < VERTICES_PER_LINE; ++i)
    {
        for (int j = 0; j < VERTICES_PER_LINE; ++j)
        {
            vec v;
            v.push_back(x_lower_left[0] + step_x*j);
            v.push_back(x_lower_left[1] + step_y*i);
            vertices.push_back(v);
        }
    }
    for (int i = 0; i < DOF_PER_ELE; ++i)
    {
        int sub_index_1 = i%(ORDER + 1);
        int sub_index_2 = i/(ORDER + 1);
        vec x;
        x.push_back(x_lower_left[0] + sub_step_x*sub_index_1);
        x.push_back(x_lower_left[1] + sub_step_y*sub_index_2);
        x_list.push_back(x);
    }    
}
// {public} functions of {Element} end here
//**************************************************************************************




//**************************************************************************************
// {public} functions of {ElementFace} start from here

ElementFace::ElementFace(int id_, int face_id_):Element(id_)
{
    id = id_;
    face_id = face_id_;
    is_boundary = false;

    switch(face_id)
    {
        case 0:
        {
            normal_vector.push_back(-1);
            normal_vector.push_back(0);
            break;              
        }
        case 1:
        {
            normal_vector.push_back(1);
            normal_vector.push_back(0);
            break;              
        }
        case 2:
        {
            normal_vector.push_back(0);
            normal_vector.push_back(-1);
            break;              
        }
        case 3:
        {
            normal_vector.push_back(0);
            normal_vector.push_back(1);
            break;              
        }      
    }

    if (id < LINE_DIVISION && face_id == 2)
    {
       is_boundary = true;
    }
    else if((id > LINE_DIVISION*(LINE_DIVISION - 1) - 1)  && face_id == 3)
    {
       is_boundary = true;
    }
    else if(id%LINE_DIVISION == 0 && face_id == 0)  
    {
        is_boundary = true;
    }
    else if(((id - LINE_DIVISION + 1)%LINE_DIVISION == 0) && face_id == 1)
    {
        is_boundary = true;
    }

    if (is_boundary == false)
    {
        switch(face_id)
        {
            case 0:
            {
                neigh_id = id - 1;
                neigh_face_id = 1;
                break;              
            }
            case 1:
            {
                neigh_id = id + 1;
                neigh_face_id = 0;
                break;              
            }
            case 2:
            {
                neigh_id = id - LINE_DIVISION;
                neigh_face_id = 3;
                break;              
            }
            case 3:
            {
                neigh_id = id + LINE_DIVISION;
                neigh_face_id = 2;
                break;              
            } 
        }
    }
}

// {public} functions of {ElementFace} end here
//**************************************************************************************