// assignment.cpp : Defines the entry point for the console application.
#include <iostream>
#include <cstring>
#include <utility> 
#include <set>
#include <stack>

#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/loop.h>
#include <igl/readOBJ.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/collapse_edge.h>
#include <igl/edge_flaps.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/snap_points.h>

#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>

#include "remesh.h"

using namespace Eigen;

void main_remesh(int argc, char**argv) {

	using namespace std;
	using namespace igl;
	using namespace Eigen;

	Eigen::MatrixXi OF, F, NF;
	Eigen::MatrixXd OV, V, NV;

	//read input
	if (argc>1)
	{
		string filename = argv[1];
		readOBJ(filename, V, F);
	}
	else
	{
		readOBJ("../data/inputs/closed/cow_head.obj", V, F);
	}

	Remesher* remesher;
	remesher = new Remesher(V, F);

	// Init the viewer
	igl::opengl::glfw::Viewer viewer;

	///*** Initialize Menu ***/
	//
	//// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	// Customize the menu
	double doubleVariable = 0.1f; // Shared between two menus

	// Draw additional windows
	menu.callback_draw_custom_window = [&]()
	{
		// Define next window position + size
		ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
		ImGui::Begin(
			"Remeshing", nullptr,
			ImGuiWindowFlags_NoSavedSettings
		);

		static double len = 0.0;

		if (ImGui::InputDouble("length", &len))
		{
			len = std::max(0.0, std::min(100.0, len));
		}

		if (ImGui::Button("Remesh", ImVec2(-1, 0)))
		{
			remesher->incremental_remeshing(len, NV, NF);
			viewer.data().clear();
			viewer.data().set_mesh(NV, NF);
		}

		if (ImGui::Button("Step", ImVec2(-1, 0)))
		{
			remesher->step_by_step(len, NV, NF);
			viewer.data().clear();
			viewer.data().set_mesh(NV, NF);
		}

		if (ImGui::Button("Adaptive", ImVec2(-1, 0)))
		{
			remesher->incremental_remeshing(len, NV, NF, true);
			viewer.data().clear();
			viewer.data().set_mesh(NV, NF);
		}

		if (ImGui::Button("Reset", ImVec2(-1, 0)))
		{
			remesher->reset();
			viewer.data().clear();
			viewer.data().set_mesh(V, F);
		}


		ImGui::End();
	};
	

	// Plot the mesh
	viewer.data().set_mesh(V, F);
	viewer.launch();
}

int main(int argc, char**argv)
{

	main_remesh(argc, argv);

	return 0;
}
