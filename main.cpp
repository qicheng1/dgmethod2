////////////////////////////////////////////////////////////////////////////////////////////
// main1.cpp :
//

#include <iostream>
#include <fstream>

#include <AMGSolver.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <Geometry.h>
#include <TemplateElement.h>
#include <FEMSpace.h>
#include <Operator.h>
#include <Functional.h>
#include <EasyMesh.h>
#include "Matrix.h"
#include "RightVec.h"

#define PI (4.0*atan(1.0))

double u(const double *);
double f(const double *);
double g(const double *);

int main(int argc, char ** argv)
{
  EasyMesh mesh;
  mesh.readData(argv[1]);

  //普通单元
  TemplateGeometry<2>	triangle_template_geometry;
  triangle_template_geometry.readData("triangle.tmp_geo");
  CoordTransform<2,2>	triangle_coord_transform;
  triangle_coord_transform.readData("triangle.crd_trs");
  TemplateDOF<2>	triangle_template_dof(triangle_template_geometry);
  triangle_template_dof.readData("triangle.DG.1.tmp_dof");
  BasisFunctionAdmin<double,2,2> triangle_basis_function(triangle_template_dof);
  triangle_basis_function.readData("triangle.DG.1.bas_fun");
  UnitOutNormal<2> triangle_unit_normal;
  triangle_unit_normal.readData("triangle.out_nrm");

  std::vector<TemplateElement<double,2,2> > template_element(1);
  template_element[0].reinit(triangle_template_geometry,
			     triangle_template_dof,
			     triangle_coord_transform,
			     triangle_basis_function,
                 triangle_unit_normal);
  //DG单元
  TemplateGeometry<1> interval_template_geometry;
  interval_template_geometry.readData("interval.tmp_geo");
  CoordTransform<1,2> interval_to2d;
  interval_to2d.readData("interval.to2d.crd_trs");

  std::vector<TemplateDGElement<1>> template_dgelement(1);
  template_dgelement[0].reinit(interval_template_geometry,interval_to2d);
  
  DGFEMSpace<double,2> dgfem_space(mesh, template_element,template_dgelement);
	
  int n_element = mesh.n_geometry(2);
  dgfem_space.element().resize(n_element);
  for (int i = 0;i < n_element;i ++)
    dgfem_space.element(i).reinit(dgfem_space,i,0);
  
  int n_boundary = mesh.n_geometry(1);
  dgfem_space.dgElement().resize(n_boundary);
  for(int i = 0; i<n_boundary; i++)
      dgfem_space.dgElement(i).reinit(dgfem_space,i,0);

  dgfem_space.buildElement();
  dgfem_space.buildDGElement();
  dgfem_space.buildDof();
  dgfem_space.buildDofBoundaryMark();

  AMatrix stiff_matrix(dgfem_space);
  stiff_matrix.set_algebric_accuracy(10);
  stiff_matrix.build_sparsity_pattern();
  stiff_matrix.build_element_matrix();
  stiff_matrix.build_dgelement_matrix();
  
  std::ofstream output("mat.txt");
  stiff_matrix.print(output);
  output.close();
   
  FEMFunction<double,2> solution(dgfem_space);
  
  RightVec right_hand_side(dgfem_space);
  right_hand_side.init();
  right_hand_side.set_algebric_accuracy(10);
  right_hand_side.build_element_rhs(&f);
  right_hand_side.build_dgelement_rhs(&g);

  std::ofstream output2("vec.txt");
  right_hand_side.print(output2,3,true,false);
  output2.close();

  /*BoundaryFunction<double,2> boundary(BoundaryConditionInfo::DIRICHLET, 1, &u);
  BoundaryConditionAdmin<double,2> boundary_admin(fem_space);
  boundary_admin.add(boundary);
  boundary_admin.apply(stiff_matrix, solution, right_hand_side);*/

  //AMGSolver solver(stiff_matrix);
  //solver.solve(solution, right_hand_side, 1.0e-08, 100000);	
  
  SolverControl solver_control(10000, 1e-08);
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(stiff_matrix,solution,right_hand_side,PreconditionIdentity());

  solution.writeOpenDXData("u.dx");
  double error = Functional::L2Error(solution, FunctionFunction<double>(&u), 10);
  std::cerr << "\nL2 error = " << error << std::endl;
  return 0;
}

double u(const double * p)
{
  return sin(PI*p[0]) * sin (2* PI*p[1])+p[0];
};

double f(const double * p)
{
  return 5*PI*PI*sin(PI*p[0]) * sin (2* PI*p[1]);
};

double g(const double * p)
{
    return u(p);
};

//
// end of file
////////////////////////////////////////////////////////////////////////////////////////////
