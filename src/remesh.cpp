#include "remesh.h"

#include <igl/edge_lengths.h>
#include <igl/edge_flaps.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/collapse_edge.h>
#include <igl/adjacency_list.h>
#include <igl/flip_edge.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/principal_curvature.h>
#include <set>
#include <cstdlib>
#include <stdlib.h>     /* abs */

Remesher::Remesher(MatrixXd V, MatrixXi F)
{
	OV = V;
	OF = F;

	CV = OV;
	CF = OF;

	counter = 0;
	CL.resize(1, 1);
	mark_feature();

}

Remesher::~Remesher()
{
}

void Remesher::reset() {
	CL.resize(1, 1);
	CV = OV;
	CF = OF;
	counter = 0;
	CV_feature.resize(OV.rows(), 1);
	CV_feature.setZero();

}

void Remesher::mark_feature() {
	//for edge topology
	MatrixXi E, EF, EI, EN;
	Eigen::VectorXi EMAP;
	MatrixXd N;

	igl::edge_flaps(OF, E, EMAP, EF, EI);//re-compute edge topology

	OE = E;
	OE_feature.resize(E.rows(), 1);
	CV_feature.resize(OV.rows(), 1);
	OE_feature.setZero();
	CV_feature.setZero();
	igl::per_vertex_normals(OV, OF, N);
	
	double angle = 3.1415926*30. / 180.;

	for (int e = 0; e < E.rows(); e++)
	{
		int f_0 = EF(e, 0), f_1 = EF(e, 1);
		int i_0 = EI(e, 0), i_1 = EI(e, 1);
		int v_2 = OF(f_0, i_0), v_3 = OF(f_1, i_1);

		Vector3d n_0 = N.row(v_2);
		Vector3d n_1 = N.row(v_3);
		double d = n_0.dot(n_1);
		if (d < -1.0) d = -1.0;
		else if (d >  1.0) d = 1.0;
		if (abs(acos(d))>angle)//feature
		{
			OE_feature(e, 0) = 1;

			int v_0 = E(e, 0), v_1 = E(e, 1);

			CV_feature(v_0, 0) = 1;
			CV_feature(v_1, 0) = 1;
		}
	}
}

void Remesher::edge_neighbor(VectorXi EMAP, MatrixXi E, MatrixXi EF, MatrixXi EI, MatrixXi &EN) {
	
	EN.resize(E.rows(), 4);

	for (int e = 0; e < E.rows(); e++)
	{
		int f_0 = EF(e, 0), f_1 = EF(e, 1);
		int i_0 = EI(e, 0), i_1 = EI(e, 1);

		int e_1 = EMAP(((i_0 + 1) % 3)*CF.rows() + f_0);
		int e_2 = EMAP(((i_0 + 2) % 3)*CF.rows() + f_0);
		int e_3 = EMAP(((i_1 + 1) % 3)*CF.rows() + f_1);
		int e_4 = EMAP(((i_1 + 2) % 3)*CF.rows() + f_1);

		EN(e, 0) = e_1;
		EN(e, 1) = e_2;
		EN(e, 2) = e_3;
		EN(e, 3) = e_4;

	}

}
void Remesher::init_target_length(double length) {
	
	MatrixXd PD1, PD2;
	VectorXd PV1, PV2;
	
	igl::principal_curvature(OV, OF, PD1, PD2, PV1, PV2);

	CL.resize(OV.rows(), 1);

	for (int i = 0; i < OV.rows(); i++)
	{
		if (PV1(i)==0)
		{
			CL(i, 0) = length;
		}
		else if (abs(PV2(i)/(PV1(i))) < 1./2.)
		{
			CL(i, 0) = length/2.;
		}
		else
		{
			CL(i, 0) = length;
		}

	}
}

void Remesher::step_by_step(double length, MatrixXd &NV, MatrixXi &NF) {

	MatrixXd L;//edge length
			   // for edge information
	MatrixXi E, EF, EI, EN;
	Eigen::VectorXi EMAP;

	double low = length * 4. / 5.;
	double high = length * 4. / 3.;

	switch (counter)
	{
	case 0://split
		   // compute edge length
		igl::edge_flaps(CF, E, EMAP, EF, EI);//compute edge topology
		igl::edge_lengths(CV, E, L);//compute edge length

		// compute edge neighbor
		edge_neighbor(EMAP, E, EF, EI, EN);

		//split edge
		split_edge(L, E, EF, EI, EN, high);
		
		counter++;
		break;

	case 1:
		//collapse short edge
		igl::edge_flaps(CF, E, EMAP, EF, EI);//re-compute edge topology
		igl::edge_lengths(CV, E, L);//compute edge length

		collapse_edge(L, E, EF, EI, EMAP, low, high);
		
		counter++;
		break;

	case 2:

		//equalize valences or flip edges
		igl::edge_flaps(CF, E, EMAP, EF, EI);//re-compute edge topology
		edge_neighbor(EMAP, E, EF, EI, EN);

		equalize_valence(E, EF, EI, EN, high);

		counter++;
		break;

	case 3:
		//tangential relaxation
		tangential_relaxation();
		counter++;
		break;
	case 4://project to surface
		project_surface();
		counter = 0;
		break;
	default:
		break;
	}

	//set current results
	NV = CV;
	NF = CF;
}

void Remesher::incremental_remeshing(double length, MatrixXd &NV, MatrixXi &NF, bool adaptive) {

	//CV = OV;
	//CF = OF;

	MatrixXd L;//edge length
	// for edge information
	MatrixXi E, EF, EI, EN;
	Eigen::VectorXi EMAP;

	double low = length * 4. / 5.;
	double high = length * 4. / 3.;

	int iter = 1;

	if (adaptive&&CL.rows()==1)
	{
		init_target_length(length);
	}

	for (int i = 0; i < iter; i++)
	{

		// compute edge length
		igl::edge_flaps(CF, E, EMAP, EF, EI);//compute edge topology
		igl::edge_lengths(CV, E, L);//compute edge length

		// compute edge neighbor
		edge_neighbor(EMAP, E, EF, EI, EN);

		//split edge
		split_edge(L, E, EF, EI, EN, high, adaptive);
		
		//collapse short edge
		igl::edge_flaps(CF, E, EMAP, EF, EI);//re-compute edge topology
		igl::edge_lengths(CV, E, L);//compute edge length

		collapse_edge(L, E, EF, EI, EMAP, low, high, adaptive);

		////equalize valences or flip edges
		igl::edge_flaps(CF, E, EMAP, EF, EI);//re-compute edge topology
		edge_neighbor(EMAP, E, EF, EI, EN);

		equalize_valence(E, EF, EI, EN, high, adaptive);

		//tangential relaxation
		tangential_relaxation();

		//project to surface
		project_surface();
		//set current results
		NV = CV;
		NF = CF;
	}


}

void Remesher::split_edge(MatrixXd L, MatrixXi E, MatrixXi EF, MatrixXi EI, MatrixXi EN, double high, bool adaptive) {
	
	assert(CV_feature.rows() == CV.rows());

	int edge_num = L.rows();
	int init_f_num = CF.rows();
	//traverse all the edges
	int e = 0;
	//bool ex
	while(e < edge_num){

		assert(L.rows() == edge_num);
		bool determin = false;
		int v_0 = E(e, 0), v_1 = E(e, 1);

		if (adaptive)
		{
			L(e, 0) > (4. / 3.)*((CL(v_0, 0)+CL(v_1, 0))/2.);
		}
		else
		{
			determin = L(e, 0) > high;
		}

		if (determin)
		{
			//split this edge
			
			/* update V */
			//add new vertex to n+1
			int f_0 = EF(e, 0), f_1 = EF(e, 1);
			int i_0 = EI(e, 0), i_1 = EI(e, 1);
			int v_2 = CF(f_0, i_0), v_3 = CF(f_1, i_1);
			Vector3d new_p = (CV.row(v_0) + CV.row(v_1)) / 2.;
			if (
				((CV.row(v_0) + CV.row(v_1)) / 2. -CV.row(v_2)).norm()==0
				||
				((CV.row(v_0) + CV.row(v_1)) / 2. - CV.row(v_3)).norm() == 0
				)
			{
				e++;
				continue;
			}

			CV.conservativeResize(CV.rows() + 1, CV.cols());
			int new_v = CV.rows() - 1;
			CV.row(new_v) = new_p;
			// determine if new vertex is also feature
			CV_feature.conservativeResize(CV_feature.rows() + 1, CV_feature.cols());
			if (CV_feature(v_0, 0) == 1 && CV_feature(v_1, 0) == 1)
			{
				CV_feature(new_v, 0) = 1;
			}
			else
			{
				CV_feature(new_v, 0) = 0;
			}

			if (adaptive)
			{
				CL.conservativeResize(CL.rows() + 1, CL.cols());
				CL(new_v, 0) = (CL(v_0, 0) + CL(v_1, 0)) / 2.;
			}

			/* updata F */ 
			//modify original 2 faces
			

			CF(f_0, (i_0 + 2) % 3) = new_v;
			CF(f_1, (i_1 + 1) % 3) = new_v;
			//add 2 new faces
			CF.conservativeResize(CF.rows() + 2, CF.cols());
			
			int new_f0 = CF.rows() - 2, new_f1 = CF.rows() - 1;

			CF.row(new_f0) = Vector3i(new_v, v_1, v_2);//f01
			CF.row(new_f1) = Vector3i(v_1, new_v, v_3);//f11


			/* update E, EF, EI, EMAP */
			E(e, 1) = new_v;
			E.conservativeResize(E.rows() + 3, E.cols());
			int e_0 = e, new_e1 = E.rows() - 3, new_e2 = E.rows() - 2, new_e3 = E.rows() - 1;
			E.row(new_e1) = Vector2i(new_v, v_1);
			E.row(new_e2) = Vector2i(new_v, v_2);
			E.row(new_e3) = Vector2i(v_3, new_v);

			

			// update L
			L(e, 0) = (CV.row(new_v) - CV.row(v_0)).norm();
			L.conservativeResize(L.rows() + 3, L.cols());
			L(new_e1, 0) = (CV.row(new_v) - CV.row(v_1)).norm();
			L(new_e2, 0) = (CV.row(new_v) - CV.row(v_2)).norm();
			L(new_e3, 0) = (CV.row(v_3) - CV.row(new_v)).norm();

			// update EF and EI, modify 4, add 3, stay 1
			int e_1 = EN(e_0, 0), e_2 = EN(e_0, 1), e_3 = EN(e_0, 2), e_4 = EN(e_0, 3);//only need take care of e1 and e4
			
			if (EF(e_1, 0) == f_0 )
			{
				EF(e_1, 0) = new_f0;
				EI(e_1, 0) = 0;
			}
			else
			{
				EF(e_1, 1) = new_f0;
				EI(e_1, 1) = 0;
			}
			
			if (EF(e_4, 0) == f_1)
			{
				EF(e_4, 0) = new_f1;
				EI(e_4, 0) = 1;
			}
			else
			{
				EF(e_4, 1) = new_f1;
				EI(e_4, 1) = 1;
			}

			EF.conservativeResize(EF.rows() + 3, EF.cols());
			EI.conservativeResize(EI.rows() + 3, EI.cols());

			EF.row(new_e1) = Vector2i(new_f0, new_f1);
			EF.row(new_e2) = Vector2i(f_0, new_f0);
			EF.row(new_e3) = Vector2i(f_1, new_f1);

			EI.row(new_e1) = Vector2i(2, 2);
			EI.row(new_e2) = Vector2i((i_0 + 1) % 3, 1);
			EI.row(new_e3) = Vector2i((i_1 + 2) % 3, 0);

			//finally, update EN, update 8
			EN(e_0, 0) = new_e2;
			EN(e_0, 3) = new_e3;

			if (EF(e_1, 0) == new_f0)
			{
				EN(e_1, 0) = new_e2;
				EN(e_1, 1) = new_e1;
			}
			else
			{
				EN(e_1, 2) = new_e2;
				EN(e_1, 3) = new_e1;
			}
			
			if (EF(e_2, 0) == f_0)
			{
				EN(e_2, 1) = new_e2;
			}
			else
			{
				EN(e_2, 3) = new_e2;
			}

			if (EF(e_3, 0) == f_1)
			{
				EN(e_3, 0) = new_e3;
			}
			else
			{
				EN(e_3, 2) = new_e3;
			}

			if (EF(e_4,0) == new_f1)
			{
				EN(e_4, 0) = new_e1;
				EN(e_4, 1) = new_e3;
			}
			else
			{
				EN(e_4, 2) = new_e1;
				EN(e_4, 3) = new_e3;
			}

			EN.conservativeResize(EN.rows() + 3, EN.cols());

			EN(new_e1, 0) = e_1;
			EN(new_e1, 1) = new_e2;
			EN(new_e1, 2) = new_e3;
			EN(new_e1, 3) = e_4;

			EN(new_e2, 0) = new_e1;
			EN(new_e2, 1) = e_1;
			EN(new_e2, 2) = e_2;
			EN(new_e2, 3) = e_0;

			EN(new_e3, 0) = e_4;
			EN(new_e3, 1) = new_e1;
			EN(new_e3, 2) = e_0;
			EN(new_e3, 3) = e_3;

			//finally, restart with this edge and add 3 more candidates
			//edge_num += 3;
			//e--;
		}
		e++;
	}
}

void Remesher::collapse_edge(MatrixXd L, MatrixXi E, MatrixXi EF, MatrixXi EI, VectorXi EMAP, double low, double high, bool adaptive) {

	typedef std::set<std::pair<double, int> > PriorityQueue;
	PriorityQueue Q;
	std::vector<PriorityQueue::iterator > Qit;
	// If an edge were collapsed, we'd collapse it to these points:
	MatrixXd C;
	int num_collapsed;

	MatrixXi F;
	MatrixXd V;

	F = CF;
	V = CV;
	igl::edge_flaps(F, E, EMAP, EF, EI);
	Qit.resize(E.rows());

	C.resize(E.rows(), V.cols());
	VectorXd costs(E.rows());
	Q.clear();
	for (int e = 0; e<E.rows(); e++)
	{
		double cost = e;
		RowVectorXd p(1, 3);
		igl::shortest_edge_and_midpoint(e, V, F, E, EMAP, EF, EI, cost, p);
		C.row(e) = p;
		Qit[e] = Q.insert(std::pair<double, int>(cost, e)).first;
	}
	num_collapsed = 0;

	bool something_collapsed = false;
	
	// collapse edge
	const int max_iter = std::ceil(0.01*Q.size());
	bool determin = false;

	if (adaptive)
	{
		int e_local = (*Q.begin()).second;
		int v_0 = E(e_local, 0), v_1 = E(e_local, 1);
		determin = (*Q.begin()).first < (4./5.) * (CL(v_0, 0) + CL(v_1, 0)) / 2.;
	}
	else
	{
		determin = (*Q.begin()).first < low;
	}

	while (determin)
	{
		if (adaptive)
		{
			int e_local = (*Q.begin()).second;
			int v_0 = E(e_local, 0), v_1 = E(e_local, 1);
			determin = (*Q.begin()).first < (4. / 5.) * (CL(v_0, 0) + CL(v_1, 0)) / 2.;
		}
		else
		{
			determin = (*Q.begin()).first < low;
		}
		int ce = (*Q.begin()).second;
		// check if collapse, larger than 4/3 length

		if (!igl::collapse_edge(
			igl::shortest_edge_and_midpoint, V, F, E, EMAP, EF, EI, Q, Qit, C))
		{
			break;
		}
		something_collapsed = true;
		num_collapsed++;
	}

	//remove 0 0 0 row
	int num_row = F.rows();
	for (int i = 0; i < F.rows(); i++)
	{
		if (i == num_row)
		{
			break;
		}
		if (F(i, 0) == 0 && F(i, 1) == 0&&F(i, 2) == 0)
		{
			num_row--;
			F.block(i, 0, F.rows() - i - 1, F.cols()) = F.bottomRows(F.rows() - i - 1);
			F.conservativeResize(num_row, F.cols());
			i--;
		}
	}

	CV = V;
	CF = F;
}

void Remesher::equalize_valence(MatrixXi E, MatrixXi EF, MatrixXi EI, MatrixXi EN, double high, bool adaptive) {

	int deviation_pre, deviation_post;
	std::vector<std::vector<int> > A;
	igl::adjacency_list(CF, A);

	for (int e = 0; e < E.rows(); e++)
	{
		int v_0 = E(e, 0), v_1 = E(e, 1);
		int f_0 = EF(e, 0), f_1 = EF(e, 1);
		int i_0 = EI(e, 0), i_1 = EI(e, 1);
		int v_2 = CF(f_0, i_0), v_3 = CF(f_1, i_1);
		int e_1 = EN(e, 0), e_2 = EN(e, 1), e_3 = EN(e, 2), e_4 = EN(e, 3);
		int e_5, e_6;
		// determine if valid
		if (EF(e_1, 0) == f_0)
		{
			e_5 = EN(e_1, 3);
		}
		else
		{
			e_5 = EN(e_1, 1);
		}

		if (EF(e_4, 0) == f_1)
		{
			e_6 = EN(e_4, 2);
		}
		else
		{
			e_6 = EN(e_4, 0);
		}
		if (e_5 == e_6)
		{
			continue;
		}

		// determine if needed
		deviation_pre = abs((int)(A[v_0].size() - 6))
			+ abs((int)(A[v_1].size() - 6))
			+ abs((int)(A[v_2].size() - 6))
			+ abs((int)(A[v_3].size() - 6));

		deviation_post = abs((int)(A[v_0].size() - 1 - 6))
			+ abs((int)(A[v_1].size() - 1 - 6))
			+ abs((int)(A[v_2].size() + 1 - 6))
			+ abs((int)(A[v_3].size() + 1 - 6));
		bool determine = false;
		if (adaptive)
		{
			determine = deviation_post < deviation_pre && (CV.row(v_2) - CV.row(v_3)).norm() < (4./3.)*(CL(v_2, 0) + CL(v_3, 0)) / 2.;

		}
		else
		{
			determine = deviation_post < deviation_pre && (CV.row(v_2) - CV.row(v_3)).norm() < high;
		}
		if (determine) // do flipping
		{

			//update A
			// update v_0
			for (std::vector<int>::iterator it = A[v_0].begin(); it != A[v_0].end(); it++)
			{
				if (*it == v_1)
				{
					it = A[v_0].erase(it);
					break;
				}
			}
			// update v_1
			for (std::vector<int>::iterator it = A[v_1].begin(); it != A[v_1].end(); it++)
			{
				if (*it == v_0)
				{
					it = A[v_1].erase(it);
					break;
				}
			}
			//update v_2 and v_3
			A[v_2].push_back(v_3);
			A[v_3].push_back(v_2);

			// update E
			E(e, 0) = v_2;
			E(e, 1) = v_3;

			// update F
			CF.row(f_0) = Vector3i(v_1, v_2, v_3);
			CF.row(f_1) = Vector3i(v_0, v_3, v_2);
			
			// ** update edge topology: EF, EI, EN, totatl 2+2+4 = 8 values to be modified, for boundary edge only 4 to be modified
			//EF(e) not change
			
			EI(e, 0) = 0;
			EI(e, 1) = 0;
			EN(e, 0) = e_4;
			EN(e, 1) = e_1;
			EN(e, 2) = e_2;
			EN(e, 3) = e_3;

			//e_1 is always near f_0, no need change, similary to e_3
			if (EF(e_2, 0) == f_0)
			{
				EF(e_2, 0) = f_1;
				EI(e_2, 0) = 1;

				EN(e_2, 0) = e_3;
				EN(e_2, 1) = e;
			}
			else
			{
				EF(e_2, 1) = f_1;
				EI(e_2, 1) = 1;

				EN(e_2, 2) = e_3;
				EN(e_2, 3) = e;
			}

			if (EF(e_4, 0) == f_1)
			{
				EF(e_4, 0) = f_0;
				EI(e_4, 0) = 1;

				EN(e_4, 0) = e_1;
				EN(e_4, 1) = e;
			}
			else
			{
				EF(e_4, 1) = f_0;
				EF(e_4, 1) = 1;

				EN(e_4, 2) = e_1;
				EN(e_4, 3) = e;
			}

			if (EF(e_1, 0) == f_0)
			{
				EI(e_1, 0) = 2;

				EN(e_1, 0) = e;
				EN(e_1, 1) = e_4;
			}
			else
			{
				EI(e_1, 1) = 2;
			
				EN(e_1, 2) = e;
				EN(e_1, 3) = e_4;
			}

			if (EF(e_3, 0) == f_1)
			{
				EI(e_3, 0) = 2;

				EN(e_3, 0) = e;
				EN(e_3, 1) = e_2;
			}
			else
			{
				EI(e_3, 1) = 2;

				EN(e_3, 2) = e;
				EN(e_3, 3) = e_2;
			}
		}
	}
}

void Remesher::tangential_relaxation() {

	assert(CV_feature.rows() == CV.rows());
	
	std::vector<std::vector<int>> A;
	igl::adjacency_list(CF, A);
	
	assert(CV.rows() == A.size());

	MatrixXd N, Q;
	Q = CV;
	igl::per_vertex_normals(CV, CF, N);

	// set q
	for (int i = 0; i < A.size(); i++)
	{
		if (A[i].size() == 0)
		{
			continue;
		}
		Vector3d p = CV.row(i);
		Vector3d q = Vector3d(0., 0., 0.);
		for (int j = 0; j < A[i].size(); j++)
		{
			q += CV.row(A[i][j]);
		}
		q /= A[i].size();
		Q.row(i) = q;
	}

	// get p'
	for (int i = 0; i < A.size(); i++)
	{
		if (A[i].size() == 0)
		{
			continue;
		}
		Vector3d p = CV.row(i);
		Vector3d q = Q.row(i);
		Vector3d n = N.row(i);
		Vector3d p_p = q + (n.dot(p - q))*n;

		if (CV_feature(i, 0) == 1)//this is a feature vertex
		{
			Vector3d proj;
			project_to_feature(p_p, proj);
			p_p = proj;
		}

		CV.row(i) = p_p;
	}
}

void Remesher::project_to_feature(Vector3d p, Vector3d &proj) {

	double min_dist = FLT_MAX, dist;
	int v_0, v_1;
	Vector3d p_0, p_1;
	for (int e = 0; e < OE.rows(); e++)
	{
		if (OE_feature(e,0) ==0)
		{
			continue;
		}

		v_0 = OE(e, 0);
		v_1 = OE(e, 1);

		p_0 = OV.row(v_0);
		p_1 = OV.row(v_1);

		//compute dist to the line
		dist = ((p - p_0).cross(p - p_1)).norm() / (p_0 - p_1).norm();

		if (dist < min_dist)
		{
			min_dist = dist;
			proj = p_0 + (p_1 - p_0) * (p_1 - p_0).dot(p - p_0) / ((p_1 - p_0).dot(p_1 - p_0));
		}
	}
}

void Remesher::project_surface() {

	MatrixXd C;
	VectorXd sqrD, I;
	igl::point_mesh_squared_distance(CV, OV, OF, sqrD, I, C);

	CV = C;
}
