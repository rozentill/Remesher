#pragma once
#include <Eigen\Core>

using namespace Eigen;

class Remesher
{
public:

	Remesher(MatrixXd V, MatrixXi F);
	~Remesher();

	void reset();

	void edge_neighbor(VectorXi EMAP, MatrixXi E, MatrixXi EF, MatrixXi EI, MatrixXi &EN);
	
	void incremental_remeshing(double length, MatrixXd &NV, MatrixXi &NF, bool adaptive = false);
	void step_by_step(double length, MatrixXd &NV, MatrixXi &NF);

	void init_target_length(double length);//for adaptive
	void mark_feature();
	void project_to_feature(Vector3d p, Vector3d &proj);
	void split_edge(MatrixXd L, MatrixXi E, MatrixXi EF, MatrixXi EI, MatrixXi EN, double high, bool adaptive=false);
	void collapse_edge(MatrixXd L, MatrixXi E, MatrixXi EF, MatrixXi EI, VectorXi EMAP, double low, double high, bool adaptive = false);
	void equalize_valence(MatrixXi E, MatrixXi EF, MatrixXi EI, MatrixXi EN, double high, bool adaptive = false);
	void tangential_relaxation();
	void project_surface();

private:

	MatrixXd CV, OV;
	MatrixXi CF, OF, CV_feature, OE, OE_feature;

	MatrixXd CL;//length for adaptive, per vertex

	int counter = 0;
};
