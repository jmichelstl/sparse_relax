/*
 * Author: Jonathan Michel, Rochester Institute of Technology
 * 
 * This code relaxes a two-dimensional eleastic network with periodic boundary
 * conditions in the horizontal network. A compression is applied, followed
 * by a shear deformation. The program reads in a network specified as a set 
 * of points, followed by a set of bonds connecting those points. Bonds are 
 * specified as pairs of integer indices into the list of points, with markers
 * indicating whether bonds should wrap around from one edge of the network to
 * the other. Next, a compression is applied, folowed by a traction on the
 * top face of the network. A sparse QR solver is employed to solve for a zero-
 * force configuaration of the network, given the applied boundary conditions.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <map>
#include <string>
#include <stdio.h>
#include <cmath>
#include "SuiteSparseQR.hpp"
#include <sstream>
#include <cfloat>
#include <stdlib.h>
#include <unistd.h>
#include <algorithm>
#include <set>
#include "cholmod.h"
#include "file_lists.hpp"
#include <time.h>
#include <cstring>

using namespace std;

//Data structures for representing vectors, bonds in the network, and 
//interacting triples.
typedef tuple<double, double> vec2d;
typedef tuple<int, int, double, char> edge_datum;
typedef tuple<int, int, int, char> triple_datum;

#define FLT_TOL 1e-10
#define TRUE 1

//Find the energy due to stretching and compression of bonds
void get_stretch_shear_energy(vector<vec2d> points, vector<edge_datum> edges, double offset, cholmod_dense *u, double ks, double mu, double &e_ss){

    double dx, dy, length, dux, duy, de;
    int iter, idx1, idx2;
    edge_datum curr_edge;

    e_ss = 0;

    for(iter = 0; iter < edges.size(); iter++){
	curr_edge = edges[iter];
        idx1 = get<0>(curr_edge);
        idx2 = get<1>(curr_edge);
	length = get<2>(curr_edge);

        dx = get<0>(points[idx2])+get<3>(curr_edge)*offset-get<0>(points[idx1]);
        dy = get<1>(points[idx2]) - get<1>(points[idx1]);
        dux =((double *)u->x)[idx2*2] - ((double *)u->x)[idx1*2];
        duy = ((double *)u->x)[idx2*2 + 1] - ((double *) u->x)[idx1*2 + 1];
	de = (dux*dx + duy*dy) / length;
	e_ss += .5 * ks * de*de;
	e_ss += .5 * mu * (dux*dux + duy*duy - de*de);
    }
}

//Find the energy due to bending of bonds
void get_bending_energy(vector<vec2d> points, vector<triple_datum> tdata, double offset, cholmod_dense *u, double kb, double &e_bend){

    double dx, dy, length_sq, dux, duy, du_c_r;
    int iter, idx1, idx2, idx3;
    triple_datum tdat;

    e_bend = 0;

    for(iter = 0; iter < tdata.size(); iter ++){
        tdat = tdata[iter];
        idx1 = get<0>(tdat);
        idx2 = get<1>(tdat);
        idx3 = get<2>(tdat);

        dx = get<0>(points[idx2])+offset*get<3>(tdat)-get<0>(points[idx1]);
        dy = get<1>(points[idx2]) - get<1>(points[idx1]);
        length_sq = dx*dx + dy*dy;
        dux = ((double *)u->x)[idx2*2] + ((double *)u->x)[idx3*2] - 2*((double *)u->x)[idx1*2];
        duy = ((double *) u->x)[idx2*2 + 1] + ((double *)u->x)[idx3*2 + 1] - 2*((double *)u->x)[idx1*2 + 1];
        du_c_r = dux * dy - duy * dx;
        e_bend += .5 * kb * du_c_r * du_c_r / length_sq;
    }
}

//Read the network from a file
//
//In:
//	-input, an open stream for a network description file, pointing to the
//	beginning of the file
//	-offset, a quantity indicating the horizontal shift that should be added
//	to wrap a bond from one edge of the network to the opposite edge. This
//	will be set if a valid header can be read.
//	-points, a list of vertices in the network, which will be filled if
//	reading from the provided file is successful
//	edges, a list of bonds in the network, which will be filled if reading
//	from the provided file is successful
//
//Out:
//	-A boolean flag indicating whether a network was successfully read from
//	the provided file
bool read_network(ifstream &input, double &offset, vector<vec2d> &points, vector<edge_datum> &edges,double &min_x, double &max_x, double &min_y, double &max_y){

    string nextline;
    int num_read, num_points, num_edges, idx1, idx2;
    double x, y, dx, dy;
    char mult;

    min_x = FLT_MAX;
    max_x = FLT_MIN;
    min_y = FLT_MAX;
    max_y = FLT_MIN;

    //Attempt to read the header for the network description file. If too little
    //information or invalid data are found, close the input stream and report
    //a failure to read a network.
    getline(input, nextline);
    num_read = sscanf(nextline.c_str(), "%d %d %lf", &num_points, &num_edges, &offset);
    if(num_read < 3 || offset < 0 || num_points <= 0 || num_edges < 0){
	cerr << "Reading failed. The provided file has an invalid header.\n";
	input.close();
	return false;
    }

    //If a valid header was read, attempt to read the specified number of points
    while(! input.eof() && points.size() < num_points){
	getline(input, nextline);
	num_read = sscanf(nextline.c_str(), "%lf %lf", &x, &y);
	if(num_read == 2){
	    points.push_back(make_tuple(x, y));
	    if(x < min_x) min_x = x;
	    if(x > max_x) max_x = x;
	    if(y < min_y) min_y = y;
	    if(y > max_y) max_y = y;
	}
    }

    //If the program has reached the end of the input stream, too little
    //information has been provided, and the file is invalid.
    if(input.eof()){
	cerr << "Reading failed. Too few points were read.\n";
	input.close();
	return false;
    }


    while(!input.eof() && edges.size() < num_edges){
	getline(input, nextline);
	num_read = sscanf(nextline.c_str(),"%d %d %hhd\n", &idx1, &idx2, &mult);

	//Ensure two point indices and an offset are given. If so, and the point
	//indices are in bounds, find the length of the edge, and add it to
	//the list.
	if(num_read == 3){
            if(min(idx1, idx2) > -1 && max(idx1, idx2) < num_points){
		dx = get<0>(points[idx2])+mult*offset - get<0>(points[idx1]);
		dy = get<1>(points[idx2]) - get<1>(points[idx1]);
		edges.push_back(make_tuple(idx1, idx2, sqrt(dx*dx+dy*dy), mult));
	    }
	    else{
	        cerr << "Reading failed due to an out-of-bounds point index.\n";
		input.close();
		return false;
	    }
	}
    }

    input.close();

    if(edges.size() < num_edges) cerr << "Too few edges were found.\n";

    return edges.size() == num_edges;
}

void assign_index(void *index_list, int itype, int index, int value){

    if(itype == CHOLMOD_LONG){
	((SuiteSparse_long *) index_list)[index] = value;
    }
    else ((int *) index_list)[index] = value;
}

//Find all triples in the network
vector<triple_datum> get_triples(vector<vec2d> points, vector<edge_datum> edges, double offset){

    //Map from each point to the points with which it shares a bond
    map<int, vector<tuple<int, char>>> adjacency_lists;
    vector<triple_datum> triples;
    int idx1, idx2, idx3, i, j;
    double x, y, dx2, dy2, dx3, dy3, dot, mag_sq;
    tuple<int, char> pair1, pair2;
    vector<tuple<int, char>> neighbors;

    for(edge_datum next_edge : edges){
	idx1 = get<0>(next_edge);
	idx2 = get<1>(next_edge);

	if(adjacency_lists.find(idx1) == adjacency_lists.end()){
	    adjacency_lists.insert(make_pair(idx1, vector<tuple<int, char>>()));
	}
	if(adjacency_lists.find(idx2) == adjacency_lists.end()){
	    adjacency_lists.insert(make_pair(idx2, vector<tuple<int, char>>()));
	}

	adjacency_lists[idx1].push_back(make_tuple(idx2, get<3>(next_edge)));
	adjacency_lists[idx2].push_back(make_tuple(idx1, -get<3>(next_edge)));
    }

    //Once the adjacency lists are built, identify all antiparallel bonds
    //sharing a common vertex. Whenever such a pair of bonds is identified,
    //create a triple, represented by the integer index of the common point,
    //or "hinge", followed by the integer indices of the two outside points,
    //followed by the offset that should be added for the bond from the hinge
    //to the first outside point.
    for(auto iter=adjacency_lists.begin();iter != adjacency_lists.end();iter++){
        
	idx1 = iter->first;
	neighbors = iter->second;
	x = get<0>(points[idx1]);
	y = get<1>(points[idx1]);

	for(i = 0; i < neighbors.size(); i ++){
	    pair1 = neighbors[i];
            idx2 = get<0>(pair1);
	    dx2 = get<0>(points[idx2]) + offset*get<1>(pair1) - x;
	    dy2 = get<1>(points[idx2]) - y;

	    for(j = i+1; j < neighbors.size(); j++){
                pair2 = neighbors[j];
		idx3 = get<0>(pair2);

		dx3 = get<0>(points[idx3]) + offset*get<1>(pair2) - x;
		dy3 = get<1>(points[idx3]) - y;

		dot = dx2*dx3 + dy2*dy3;
		mag_sq = (dx2*dx2 + dy2*dy2)*(dx3*dx3 + dy3*dy3);

		//If the dot product of the vectors pointing from the hinge
		//to either end point is negative, and the square of this
		//dot product is equal to the product of the square of the
		//magnitudes of the two displacement vectors, the two 
		//displacement vectors are antiparallel, and a valid triple has
		//been identified.
		if(dot < 0 && abs(dot*dot  - mag_sq) < FLT_TOL){
		    triples.push_back(make_tuple(idx1,idx2,idx3,get<1>(pair1)));
		}
	    }
	}
    }

    return triples;
}

//Auxilliary function to update an entry in a map from an index pair to a
//sparse matrix element corresponding to that index pair.
void update_matrix(int idx1, int idx2, double value, map<tuple<int, int>, double> &sp_mat){

    tuple<int, int> idx_pair = idx1<idx2 ?  make_tuple(idx1,idx2) : make_tuple(idx2,idx1);

    if(sp_mat.find(idx_pair) == sp_mat.end()){
	sp_mat.insert(make_pair(idx_pair, value));
    }
    else sp_mat[idx_pair] += value;
}

//Calculate contributions to a stiffness matrix due to a pair of vertices 
//joined by a bond.
void calc_stretching_elements(vector<vec2d> points, vector<edge_datum> edges, double offset, double stiffness, map<tuple<int, int>, double> &sp_mat){

    double dx, dy, l, inv_l_cubed, mat_elem;
    vec2d p1, p2;
    int xIndex1, xIndex2, yIndex1, yIndex2;

    for(edge_datum edat : edges){
	p1 = points[get<0>(edat)];
	p2 = points[get<1>(edat)];
	dx = get<0>(p2) + offset*get<3>(edat) - get<0>(p1);
	dy = get<1>(p2) - get<1>(p1);

	l = sqrt(dx*dx + dy*dy);
        inv_l_cubed = 1 / (l * l * l);

        xIndex1 = 2*get<0>(edat);
        xIndex2 = 2*get<1>(edat);
        yIndex1 = 2*get<0>(edat) + 1;
        yIndex2 = 2*get<1>(edat) + 1;

	//Add matrix elements corresponding to partial derivatives with respect
	//to x coordinates
	mat_elem = stiffness * inv_l_cubed * dx * dx;
        update_matrix(xIndex1, xIndex2, -mat_elem, sp_mat);
        update_matrix(xIndex1, xIndex1, mat_elem, sp_mat);
        update_matrix(xIndex2, xIndex2, mat_elem, sp_mat);

	//Add matrix elements corresponding to partial derivatives with respect
	//to one x coordinate and one y coordinate.
	mat_elem = stiffness * inv_l_cubed * dx * dy;
	update_matrix(xIndex1, yIndex2, -mat_elem, sp_mat);
	update_matrix(yIndex1, xIndex2, -mat_elem, sp_mat);
	update_matrix(xIndex1, yIndex1, mat_elem, sp_mat);
	update_matrix(xIndex2, yIndex2, mat_elem, sp_mat);

	//Add matrix elements corresponding to derivatives with respect to y
	//coordinates
	mat_elem = stiffness * inv_l_cubed * dy * dy;
	update_matrix(yIndex1, yIndex2, -mat_elem, sp_mat);
	update_matrix(yIndex1, yIndex1, mat_elem, sp_mat);
	update_matrix(yIndex2, yIndex2, mat_elem, sp_mat);
    }

}	

//Calculate the contributions to a stiffness matrix due to a triple of vertices
void calc_bending_elements(vector<vec2d> points, vector<triple_datum> triples, double offset, double kb, map<tuple<int, int>, double> &sp_mat){

    double dx, dy, xx, yy, xy, len_sq;
    int idx1, idx2, idx3;

    for(triple_datum tdat : triples){
        //Find the indices of the three vertices concerned in a given bending
        //interaction
	idx1 = get<0>(tdat);
	idx2 = get<1>(tdat);
	idx3 = get<2>(tdat);

        //To evaluate derivatives of the bending energy with respect to the
        //hinge and terminal coordinates, first obtain the components of the
        //displacement vector from the hinge to one of the end points.
        dx = get<0>(points[idx2]) + offset*get<3>(tdat) - get<0>(points[idx1]);
        dy = get<1>(points[idx2]) - get<1>(points[idx1]);
	len_sq = dx*dx + dy*dy;
	xx = kb*dx*dx / len_sq;
	yy = kb*dy*dy / len_sq;
	xy = kb*dx*dy / len_sq;

	//Add diagonal elements of the stiffness matrix
	update_matrix(2*idx1, 2*idx1, 4*yy, sp_mat);
	update_matrix(2*idx2, 2*idx2, yy, sp_mat);
	update_matrix(2*idx3, 2*idx3, yy, sp_mat);
	update_matrix(2*idx1 + 1, 2*idx1 + 1, 4*xx, sp_mat);
	update_matrix(2*idx2 + 1, 2*idx2 + 1, xx, sp_mat);
	update_matrix(2*idx3 + 1, 2*idx3 + 1, xx, sp_mat);
	
	//Add elements of the stiffness matrix involving the x component of
	//the hinge coordinate and another coordinate
	update_matrix(2*idx1, 2*idx1 + 1, -4*xy, sp_mat);
	update_matrix(2*idx1, 2*idx2, -2*yy, sp_mat);
	update_matrix(2*idx1, 2*idx3, -2*yy, sp_mat);
	update_matrix(2*idx1, 2*idx2 + 1, 2*xy, sp_mat);
	update_matrix(2*idx1, 2*idx3 + 1, 2*xy, sp_mat);

	//Add elements of the stiffness matrix inolving the y component of the
	//the hinge coordinate and another coordinate.
	update_matrix(2*idx1 + 1, 2*idx2, 2*xy, sp_mat);
	update_matrix(2*idx1 + 1, 2*idx3, 2*xy, sp_mat);
	update_matrix(2*idx1 + 1, 2*idx2 + 1, -2*xx, sp_mat);
	update_matrix(2*idx1 + 1, 2*idx3 + 1, -2*xx, sp_mat);

	//Add elements of the stiffness matrix involoving the x component of
	//the first end point and one other coordinate
	update_matrix(idx2*2, idx2*2 + 1, -xy, sp_mat);
	update_matrix(idx2*2, idx3*2, yy, sp_mat);
	update_matrix(idx2*2, idx3*2 + 1, -xy, sp_mat);

	//Add elements of the stiffness matrix involving the y component of
	//the first end point and another coordinate
	update_matrix(idx2*2 + 1, idx3*2, -xy, sp_mat);
	update_matrix(idx2*2 + 1, idx3*2 + 1, xx, sp_mat);

	//Add elements of the stiffness matrix involving the x component of the
	//second end point and the y component of the second end point
	update_matrix(idx3*2, idx3*2 + 1, -xy, sp_mat);
    }
}

//Calculate the contributions to the stiffness matrix due to shearing of the
//background gel to which the network is coupled.
void calc_shearing_elements(vector<vec2d> points, vector<edge_datum> edges, double offset, double shear_mod, map<tuple<int, int>, double> &sp_mat){

     int idx1, idx2;
     vec2d p1, p2;

     //Components of displacement vectors, components of unit vectors, the
     //length of a bond, and the finalized matrix element to be used to
     //update a network
     double dx, dy, rh_x_sq, rh_y_sq, inv_len_sq, mat_elem;

     for(edge_datum edat : edges){
	 idx1 = get<0>(edat);
	 idx2 = get<1>(edat);
	 p1 = points[idx1];
	 p2 = points[idx2];

         //Find components of displacement vectors, components of unit vectors,
	 //and the length of a bond.
	 dx = get<0>(p2) + get<3>(edat) * offset - get<0>(p1);
	 dy = get<1>(p1) - get<1>(p2);
	 inv_len_sq = get<2>(edat);
	 inv_len_sq = 1 / inv_len_sq / inv_len_sq;
	 rh_x_sq = dx*dx * inv_len_sq;
	 rh_y_sq = dy*dy * inv_len_sq;

         //Compute derivatives of the strain energy with respect to the
	 //components of the displacement field, and update a sparse matrix
	 //representation of the stiffness matrix.
	 mat_elem = shear_mod * (1 - rh_x_sq);
	 update_matrix(idx1*2, idx1*2, mat_elem, sp_mat);
	 update_matrix(idx2*2, idx2*2, mat_elem, sp_mat);
	 update_matrix(idx1*2, idx2*2, -mat_elem, sp_mat);

	 mat_elem = shear_mod * (1 - rh_y_sq);
	 update_matrix(idx1*2 + 1, idx1*2 + 1, mat_elem, sp_mat);
	 update_matrix(idx2*2 + 1, idx2*2 + 1, mat_elem, sp_mat);
	 update_matrix(idx1*2 + 1, idx2*2 + 1, -mat_elem, sp_mat);
     }

}

//Given data for a stiffness matrix, specified as ordered pairs of indices
//mapped to matrix elements for those index pairs, construct a sparse 
//representation of the stiffness matrix for a network.
cholmod_sparse *get_sparse_stiffness_mat(vector<vec2d> points, vector<edge_datum> edges, vector<triple_datum> &tdata, double offset, double ks, double kb, double mu, cholmod_common *cc){

    map<tuple<int, int>, double> sparse_map;
    cholmod_triplet *kmat_triplet_form;
    cholmod_sparse *kmat_compressed;
    int elem_num, kitype, iter;

    if(kb > 0){
        tdata = get_triples(points, edges, offset);
    }

    //Obtain stiffness matrix elements related to stretching and bending of
    //filaments, and shearing of the background gel
    calc_stretching_elements(points, edges, offset, ks, sparse_map);
    if(kb > 0){
        calc_bending_elements(points, tdata, offset, kb, sparse_map);
    }

    if(mu > 0){
        calc_shearing_elements(points, edges, offset, mu, sparse_map);
    }

    //Create a sparse matrix representation of the stiffness matrix of the
    //form row index, column index, value
    kmat_triplet_form = cholmod_l_allocate_triplet(2*points.size(), 2*points.size(), sparse_map.size(), -1, CHOLMOD_REAL, cc);
    kitype = kmat_triplet_form->itype;

    elem_num = 0;
    for(auto iter = sparse_map.begin(); iter != sparse_map.end(); iter ++){
	if(abs(iter->second) > FLT_TOL){
	    assign_index(kmat_triplet_form->i, kitype, elem_num, get<0>(iter->first));
	    assign_index(kmat_triplet_form->j, kitype, elem_num, get<1>(iter->first));
	    ((double *) kmat_triplet_form->x)[elem_num] = iter->second;
	    elem_num++;
	}
    }

    kmat_triplet_form->nnz = elem_num;

    //Convert the sparse matrix representation of the stiffness matrix to
    //compressed row format
    kmat_compressed = cholmod_l_triplet_to_sparse(kmat_triplet_form,elem_num,cc);

    cholmod_l_free_triplet(&kmat_triplet_form, cc);

    return kmat_compressed;
}

//Utility function to find the number of digits in the base-ten representation
//of an integer, for the purpose of appending a zero-padded integer to the
//base name for a report file
int get_dec_digits(int value){
    int num_digits = 0;

    do{
        num_digits++;
        value /= 10;
    }while(value > 0);

    return num_digits;
}

//Utility function to create a file name, given a base, the number of decimal
//digits to use in an appended count, and the position of the report file in
//a series.
string report_name(string base, int num_digits, int count, string ext){
    int cdigits, padding, iter;
    ostringstream oss;

    padding = num_digits - get_dec_digits(count);
    oss << base << "_";
    for(iter = 0; iter < padding; iter++) oss << 0;
    oss << count << "." << ext;

    return oss.str();
}

//Displace the top vertices of a network to apply either a compression or
//a shear
void add_displacement(cholmod_dense *rhs, vector<int> top, double disp, char mode){

    //If the mode of displacement is compression, the displacement should
    //be downward. Otherwise, it should be rightward.
    int offset = mode == 'c' ? 1 : 0;
    double factor = mode == 'c' ? -1 : 1;

    for(int next_index : top){
	((double *) rhs->x)[next_index * 2 + offset] = factor * disp;
    }
}

//Create an affine displacement field, with some compressive strain, c and a 
//shear displacement, s, such that the x and y displacements of each point are
//s * (y - ymin) and c * (y-ymin), respectively.
void affine_displacement(cholmod_dense *rhs, vector<vec2d> points, double c_strain, double s_strain, double ymin){

    int iter;
    double y;

    for(iter = 0; iter < points.size(); iter ++){
	y = get<1>(points[iter]);
	((double *) rhs->x)[iter * 2] = s_strain * (y - ymin);
	((double *) rhs->x)[iter * 2 + 1] = c_strain * (ymin - y);
    }
}

//Print the components of a displacement field to a file
void report_displacements(string file_name, cholmod_dense *u, int npoints){

    FILE *report_file;
    int iter;

    report_file = fopen(file_name.c_str(), "w");
    for(iter = 0; iter < npoints; iter++){
	fprintf(report_file, "%2.10le\t%2.10le\n", ((double *) u->x)[iter*2], ((double *) u->x)[iter*2 + 1]);
    }

    fclose(report_file);
}
	
//Prepare to pose the problem of finding a zero-force configuration for non-
//constrained nodes in a network as a matrix problem. Define operators
//projecting a vector from a 2N-dimensional space to a 2R-dimensional space
//and vice-versa, where N and R are the total number and the number of relaxed
//vertices in the network, respectively. Use these projection operators to
//define the left-hand and right-hand side of the matrix equation.
void get_lhs_rhs(int num_points, set<int> bottom, set<int> top, cholmod_sparse *kmat, cholmod_sparse **lhs, cholmod_sparse **rhs, cholmod_sparse **sp_prn, cholmod_common *cc){

    cholmod_triplet *prn;
    cholmod_sparse *sp_pnr;
    int nr, iter, pos, kitype;

    //Find the number of vertices to be relaxed, and construct a sparse
    //matrix representation of the projection from the reduced to the full-
    //dimensional space in triplet form
    nr = 2*(num_points - bottom.size() - top.size());
    prn = cholmod_l_allocate_triplet(2*num_points, nr, nr, 0, CHOLMOD_REAL, cc);
    prn->nnz = prn->nzmax;
    kitype = prn->itype;

    pos = 0;
    for(iter = 0; iter < num_points && pos < nr/2; iter++){
	//If a point is neither on the bottom nor on the top of the network,
	//it should be accounted for in the mapping
        if(bottom.find(iter) == bottom.end() && top.find(iter) == top.end()){
	    assign_index(prn->i, kitype, 2*pos, 2*iter);
	    assign_index(prn->j, kitype, 2*pos, 2*pos);
	    ((double *) prn->x)[2*pos] = 1;
 
	    assign_index(prn->i, kitype, 2*pos + 1, 2*iter + 1);
	    assign_index(prn->j, kitype, 2*pos + 1, 2*pos + 1);
	    ((double *) prn->x)[2*pos + 1] = 1;

	    pos ++;
	}
    }

    //Convert the previously defined projection matrix to compressed form, and
    //calculated its transpose.
    *sp_prn = cholmod_l_triplet_to_sparse(prn, prn->nnz, cc);
    sp_pnr = cholmod_l_transpose(*sp_prn, 1, cc);
    //cholmod_l_print_sparse(sp_pnr, "PN->R", cc);

    //Define the final lhs and rhs matrices
    *rhs = cholmod_l_ssmult(sp_pnr, kmat, 0, TRUE, TRUE, cc);
    //cholmod_l_print_sparse(*rhs, "RHS", cc);
    *lhs = cholmod_l_ssmult(*rhs, *sp_prn, 0, TRUE, TRUE, cc);
    //cholmod_l_print_sparse(*lhs, "LHS", cc);

    cholmod_l_free_triplet(&prn, cc);
    cholmod_l_free_sparse(&sp_pnr, cc);
}

//Given a network to compress and shear, the mechanical, attributes of the 
//network, and the program of compression and shear that should be executed,
//subject the network to the desired strain protocol, perform structural
//relaxation, and report the results.
void do_relaxation_run(string dat_name, string log_name, string disp_base, string e_name, double ks, double kb, double mu, double min_c, double max_c, double c_inc, double min_s, double max_s, double s_inc, bool strain_mode, bool affine, cholmod_common *cc){

    //Basic attributes of the network
    double offset, min_x, max_x, min_y, max_y, pe, e_ss, e_bend, compression;
    double shear, sqrt_dof;
    int num_read, iter, elems, num_digits, disp_count, mitype, sitype;
    time_t start;
    //Structures for obtaining a description of the network
    ifstream input;
    vector<vec2d> points;
    vector<edge_datum> edges;
    vector<int> bottom, top;
    vector<triple_datum> tdata;

    //Arrays to keep track of the displacements, and specify the structural
    //relaxation problem to be solved.
    cholmod_dense *u, *u_final, *rhs, *bc_vec;

    //Sparse matrices for computing the psuedo-forces on the degrees of
    //freedom to be relaxed
    cholmod_sparse *kmat, *cond_mat, *mult, *prn;
    double unity[2] = {1, 0}, neg1[2] = {-1, 0}, null[2] = {0, 0};
    bool write_disps = false;
    FILE *log = NULL, *efile = NULL;

    //If a non-empty log file name was passed, open a log file for reading
    if(! log_name.compare("") == 0){
        log = fopen(log_name.c_str(), "w");
    }

    //If a log file is to be written, make a record of the elapsed time
    if(log != NULL) time(&start);

    input.open(dat_name);

    //If no file was read, the program cannot proceed
    if(! input.is_open()) return;

    //Ensure network data can be read
    if(! read_network(input, offset, points, edges, min_x,max_x, min_y, max_y)){
        cerr << "A properly formatted network description could not be read.\n";
	return;
    }

    //Identify the top and bottom nodes in the network
    for(iter = 0; iter < points.size(); iter ++){
        if(get<1>(points[iter]) < min_y + FLT_TOL) bottom.push_back(iter);
	if(get<1>(points[iter]) > max_y - FLT_TOL) top.push_back(iter);
    }

    set<int> bset(bottom.begin(), bottom.end());
    set<int> tset(top.begin(), top.end());

    //If displacements are specified as strains, convert to absolute 
    //displacement
    if(strain_mode){
        min_c *= (max_y - min_y);
        max_c *= (max_y - min_y);
        c_inc *= (max_y - min_y);
	if(c_inc != c_inc) c_inc = FLT_MAX - max_c;

	min_s *= (max_x - min_x);
	max_s *= (max_x - min_x);
	s_inc *= (max_x - min_x);
	if(s_inc != s_inc) s_inc = FLT_MAX - max_c;
    }

    //Build the stiffness matrix
    kmat = get_sparse_stiffness_mat(points,edges,tdata, offset, ks, kb, mu, cc);

    //If displacement files are to be written, prepare to generate file names
    if(disp_base.compare("") != 0){
	write_disps = true;
	disp_count = 1;
	num_digits = get_dec_digits(2+(int)((max_c-min_c)/c_inc)+(int)((max_s-min_s)/s_inc));
    }

    //If a non-empty energy file name was passed, open an energy file for
    //reading
    if(e_name.compare("") != 0){
        efile = fopen(e_name.c_str(), "w");
    }

    //Obtain matrices for constructing the LHS and RHS of the linear system
    //of equations to be solved for u.
    get_lhs_rhs(points.size(), bset, tset, kmat, &cond_mat, &mult, &prn, cc);

    //Initialize the RHS of the linear system of equations
    bc_vec = cholmod_l_zeros(2*points.size(), 1, CHOLMOD_REAL, cc);
    u_final = cholmod_l_zeros(2*points.size(), 1, CHOLMOD_REAL, cc);
    rhs = cholmod_l_zeros(2*points.size() - 2*bottom.size() - 2*top.size(), 1, CHOLMOD_REAL, cc);

    sqrt_dof = sqrt(rhs->nrow);

    //Apply compression and perform relaxation
    compression = min_c;
    while(compression < max_c + FLT_TOL){
	shear = min_s;
	while(shear < max_s + FLT_TOL){
            //Impose compression and shear
	    if(affine){
                affine_displacement(bc_vec, points, compression/(max_y-min_y),  shear/(max_y - min_y), min_y);
            }
	    else{
                add_displacement(bc_vec, top, compression, 'c');
                add_displacement(bc_vec, top, shear, 's');
	    }
            cholmod_l_sdmult(mult, 0, neg1, null, bc_vec, rhs, cc);

            u = SuiteSparseQR <double>(cond_mat, rhs, cc);

	    //Project the displacement vector for the reduced-dimension
	    //subspace consisting of just internal vertices' displacements on
	    //to the full 2*N-dimensional displacement space, and sum this
	    //result with the boundary condition vector.
	    cholmod_l_copy_dense2(bc_vec, u_final, cc);
            cholmod_l_sdmult(prn, 0, unity, unity, u, u_final, cc);

	    //If desired, report the displacement field after relaxation
	    if(write_disps){
                report_displacements(report_name(disp_base, num_digits, disp_count, "disp"), u_final, points.size());
	        disp_count ++;
	    }

	    //Find the magnitude of the difference between the LHS and RHS
	    //of the linear system, by multiplying the 2Rx1 displacement vector
	    //u by the LHS, and subtracting the result from the RHS
	    cholmod_l_sdmult(cond_mat, 0, neg1, unity, u, rhs, cc);
	    
	    if(log != NULL){
	        fprintf(log, "Residual RMS force: %2.10le\n", cholmod_l_norm_dense(rhs, 2, cc));
	    }
            else{
	        printf("Residual RMS force: %2.10le\n", cholmod_l_norm_dense(rhs, 2, cc));
	    }

	    if(ks > 0 || mu > 0){
                get_stretch_shear_energy(points,edges,offset,u_final,ks,mu,e_ss);
	    }
	    if(kb > 0){
                get_bending_energy(points, tdata, offset, u_final, kb, e_bend);
	    }
	    pe = e_ss + e_bend;
            if(efile != NULL){
	        fprintf(efile, "%12.8lf\t%12.8lf\t%12.10le\n", compression, shear, pe);
	    }
	    else{
	        printf("%12.8lf\t%12.8lf\t%12.10le\n", compression, shear, pe);
	    }
            shear += s_inc;
	
	    cholmod_l_free_dense(&u, cc);
	}

	compression += c_inc;
    }

    //Clean up allocated memory and exit
    cholmod_l_free_sparse(&kmat, cc);
    cholmod_l_free_sparse(&cond_mat, cc);
    cholmod_l_free_sparse(&mult, cc);
    cholmod_l_free_sparse(&prn, cc);
    cholmod_l_free_dense(&rhs, cc);
    cholmod_l_free_dense(&bc_vec, cc);
    cholmod_l_free_dense(&u_final, cc);
    cholmod_l_finish(cc);
    if(efile != NULL) fclose(efile);
    if(log != NULL){
	fprintf(log, "Elapsed time: %d seconds.\n", (int) difftime(time(NULL), start));
	fclose(log);
    }
}

//Obtain the range for either shear or compression
void get_loading_range(double &min, double &max, double &inc, string ltype){

    int num_read;
    string response;

    do{
        cout << "Enter the minimum value, maximum value and increment for " << ltype << ": ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf %lf %lf", &min, &max, &inc);
        if(num_read < 0 || min < 0 || max < 0 || max < min || inc <0){
            cerr << "Enter three numbers, with a max greater than or equal to the min.\n";
        }
        else break;
    }while(true);
}

//Split a line into pieces, using whitespace characters as delimiters
vector<string> split_line(string line){
    vector<string> result;
    char *full_line, *token;

    full_line = (char *) malloc((1 + line.size()) * sizeof(char));
    strcpy(full_line, line.c_str());

    token = strtok(full_line, " \t");
    while(token != NULL){
	result.push_back(string(token));
	token = strtok(NULL, " \t");
    }

    free(full_line);

    return result;
}

//Parse a file containing instructions for a batch of simulations. The input
//file should countain the following:
//	-base name for network description (.dat) files
//	-mechanical attributes of network fibers and the background gel
//	-the range of compressions to apply
//	-the range of shears to apply
//Information should be stated in the order given above. Blank lines, and lines
//beginning with a pound sign, will be ignored.
//
//Return value:
//	-true, if all necessary information could be obtained
//	-false, otherwise
bool parse_batch_file(string file_name, string &base_name, string &tag, double &ks, double &kb, double &mu, double &cmin, double &cmax, double &cinc, double &smin, double &smax, double &sinc){

    string nextline;
    ifstream input;
    vector<string> tokens;
    int num_read;
    bool all_found;

    //Make sure the batch file can be opened for reading
    input.open(file_name);
    if(! input.is_open()){
	cerr << "The file " << file_name << " could not be read.\n";    
	return false;
    }

    //Read the base name and an optional tag
    while(! input.eof()){
        getline(input, nextline);
	if(nextline.compare("") == 0 || nextline[0] == '#') continue;

	tokens = split_line(nextline);

	//If at least one string of non-whitespace characters is found,
	//initialize the base name, and, optionally, a tag to append to the
	//base name.
        if(tokens.size() >= 1){
	    base_name = tokens[0];

	    if(tokens.size() > 1) tag = tokens[1];

	    break;
	}
    }

    //If the end of the file has been reached, not enough information could
    //be read, and reading has failed.
    if(input.eof()){
	cerr << "Not enough information was found.\n";
	input.close();
	return false;
    }

    //Obtain the mechanical attributes of the network and background gel
    while(! input.eof()){
	getline(input, nextline);
	if(nextline.compare("") == 0 || nextline[0] == '#') continue;

	num_read = sscanf(nextline.c_str(), "%lf %lf %lf", &ks, &kb, &mu);
	if(num_read < 3 || ks < 0 || kb < 0 || mu < 0){
	    cerr << "Error: invalid format for mechanical attributes.";
	    input.close();
	    return false;
	}

	else break;
    }

    //Make sure there is still something to read
    if(input.eof()){
	cerr << "Not enough information was provided in the batch file.\n";
	input.close();
	return false;
    }

    //Obtain the parameters for compressing the network
    while(! input.eof()){
	getline(input, nextline);
	if(nextline.compare("") == 0 || nextline[0] == '#') continue;

	num_read = sscanf(nextline.c_str(), "%lf %lf %lf", &cmin, &cmax, &cinc);
	if(num_read < 3 || cmin < 0 || cmax < 0 || cinc < 0 || cmax < cmin){
	    cerr << "Error: invalid specification of compression.";
	    input.close();
	    return false;
	}

	else break;
    }

    //Make sure there is still something to read
    if(input.eof()){
	cerr << "Not enough information was provided in the batch file.\n";
	input.close();
	return false;
    }

    //Obtain the parameters for shearing the network
    while(! input.eof()){
	getline(input, nextline);
	if(nextline.compare("") == 0 || nextline[0] == '#') continue;

	num_read = sscanf(nextline.c_str(), "%lf %lf %lf", &smin, &smax, &sinc);
	if(num_read < 3 || smin < 0 || smax < 0 || sinc < 0 || smax < smin){
	    cerr << "Error: invalid format for mechanical attributes.";
	    input.close();
	    return false;
	}

	else{
	    all_found = true;
	    break;
	}
    }

    input.close();
    return all_found;
}

//Prompt the user for a base name for a batch of dat files, mechanical 
//information about the network, and a testing protocol, and then perform
//a batch of simulations using these settings.
void run_relaxation_batch(string file_name, bool strain_mode, bool affine){

    string base, tag;
    //Lists of file names for i/o
    vector<string> dat_files, log_files, disp_files, edat_files;

    //Network mechanical attributes and loading parameters
    double ks, kb, mu, min_c, max_c, c_inc, min_s, max_s, s_inc;

    //Data for allocating and manipulating vectors and matrices
    size_t total_mem, available_mem;
    cholmod_common common;
    int iter;
    bool code;

    //Attempt to parse the specified batch file to obtain information about the
    //simulation run
    if(! parse_batch_file(file_name, base, tag, ks, kb, mu, min_c, max_c, c_inc, min_s, max_s, s_inc)){
	cerr << "File parsing failed.\n";
	return;
    }

    //Generate lists of dat file and output file names
    code = get_file_lists(base, tag, dat_files, log_files, disp_files, edat_files);

    if(! code){
	cerr << "File lists could not be read.\n";
	return;
    }

    //Prepare data structures for creating matrices
    cholmod_l_start(&common);

    //If the starting and ending values are the same for either compression
    //or shear, just make the increment arbitrarily large
    if(max_s - min_s < FLT_TOL) s_inc = 1e12;
    if(max_c - min_c < FLT_TOL) c_inc = 1e12;

    //Perform network relaxation
    for(iter = 0; iter < dat_files.size(); iter++){
        do_relaxation_run(dat_files[iter], log_files[iter], disp_files[iter], edat_files[iter], ks, kb, mu, min_c,  max_c, c_inc, min_s, max_s, s_inc, strain_mode, affine, &common);
    }

    //Clean up the sparse matrix workspace
    cholmod_l_finish(&common);
}

//Prompt for a report file name
string get_report_name(string prompt){

    string response;

    cout << prompt << ", or enter to decline: ";
    getline(cin, response);
    return response;
}

//Perform an interactive session, in which the user is prompted for a specific
//dat file, and asked about options for reporting results.
void run_relaxation_interactive(bool strain_mode, bool affine){

    //Strings for file i/o
    string response, dat_name, log_name, disp_base, e_name;

    //Network mechanical attributes and loading parameters
    double ks, kb, mu, min_c, max_c, c_inc, min_s, max_s, s_inc;

    //Data for allocating and manipulating matrices and vectors
    size_t total_mem, available_mem;
    cholmod_common common;
    int num_read;

    cholmod_l_start(&common);

    //Prompt for the file name containing the network description
    cout << "Enter the name of the network description file: ";
    getline(cin, dat_name);
   
    //Prompt for the mechanical attributes of the network
    do{
        cout << "Enter the stretching modulus, bending modulus, and background shear modulus: ";
        getline(cin, response);
        num_read = sscanf(response.c_str(), "%lf %lf %lf", &ks, &kb, &mu);
        if(num_read < 3 || ks < 0 || kb < 0 || mu < 0){
            cerr << "Enter three non-negative numbers.\n";
        }
        else break;
    }while(true);

    //Find the range of compression and shear
    get_loading_range(min_c, max_c, c_inc, "compression");
    get_loading_range(min_s, max_s, s_inc, "shear");

    //If the starting and ending values are the same for either compression
    //or shear, just make the increment arbitrarily large
    if(max_s - min_s < FLT_TOL) s_inc = 1e12;
    if(max_c - min_c < FLT_TOL) c_inc = 1e12;

    //Offer the opportunity to write a log, and report displacements and
    //energies to files
    log_name = get_report_name("Type a log name");
    disp_base = get_report_name("Type a base name for displacement files");
    e_name = get_report_name("Type an energy file");

    //Perform network relaxation
    do_relaxation_run(dat_name, log_name, disp_base, e_name, ks, kb, mu, min_c,  max_c, c_inc, min_s, max_s, s_inc, strain_mode, affine, &common);
    
    //Free data associated with linear algebra computations
    cholmod_l_finish(&common);
}

int main(int argc, char **argv){

    char c;
    char *batch_file;
    bool batch_mode = false, strain_mode = false,affine = true;

    //Check for a command line flag. If the "-b" option is specified, along
    //with a non-empty string, attempt to read information from a file with
    //the name of the provided argument, and run the program in batch mode.
    //If "-b" is specified with no argument, advise the user of incorrect
    //usage and exit.
    while((c = getopt(argc, argv, "b:sz")) != -1){
	switch(c) {
	    case 'b':
		batch_file = (char *) malloc((1 + strlen(optarg))*sizeof(char));
                strcpy(batch_file, optarg);
		batch_mode = true;
		break;
	    case 's':
		//Specify displacements as relative strains, rather than 
		//absolute quantities
                strain_mode = true;
		break;
	    //Relaxation should be performed about a displacement field of
	    //zero, rather than an affine displacement field
	    case 'z':
                affine = false;
		break;
	    case '?' :
	        if(optopt == 'b'){
		    cerr << "Option \"b\" requires a file name.\n";
		    return -1;
	        }
	        else if(isprint(optopt)){
		    fprintf(stderr, "Unrecognized option: %c\n", optopt);
	        }
	        else{
		    cerr << stderr, "Unknown option character.\n";
	        }
	        break;
    	    default:
	        break;
	}
    }

    //Perform energy relaxation of a strained elastic, filamentous lattice
    //with horizontal periodic boundary conditions.
    if(batch_mode){
	string bstring = string(batch_file);
	free(batch_file);
	run_relaxation_batch(bstring, strain_mode, affine);
    }
    else{
        run_relaxation_interactive(strain_mode, affine);
    }

    return 0;
}
