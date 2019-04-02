#include "RightVec.h"

void RightVec::init()
{
    this->reinit(space->n_dof());
}
void RightVec::build_element_rhs(double (*f)(const double *))
{
    auto the_ele = space->beginElement();
    const auto end_ele = space->endElement();
    for(; the_ele != end_ele; the_ele++)
    {
        std::vector<int>& ele_dof = the_ele->dof();
        double volume = the_ele->templateElement().volume();
        const QuadratureInfo<2>& quad_info = the_ele->findQuadratureInfo(algebric_accuracy);
        std::vector<point_t> q_point = the_ele->local_to_global(quad_info.quadraturePoint());
        std::vector<double> jacobian = the_ele->local_to_global_jacobian(quad_info.quadraturePoint());
        int n_quadrature_point = q_point.size();
        std::vector<std::vector<double>> basis_value = the_ele->basis_function_value(q_point);
        for(int l = 0; l < n_quadrature_point; l++)
        {
            double Jxw = quad_info.weight(l)*jacobian[l]*volume;
            double f_value = (*f)(q_point[l]);
            for(int j = 0; j<ele_dof.size(); j++){
                this->operator[](ele_dof[j]) += Jxw*f_value*basis_value[j][l];
            }
        }
    }
}
void RightVec::build_dgelement_rhs(double (*g)(const double *))
{
   auto the_dgele = space->beginDGElement();
   const auto end_dgele = space->endDGElement();
   for(;the_dgele != end_dgele; the_dgele++)
    {
        if((the_dgele->p_neighbourElement(1))!=NULL) continue;
        element_t& nei = the_dgele->neighbourElement(0);
        double volume = the_dgele->templateElement().volume();
        const QuadratureInfo<1>& quad_info = the_dgele->findQuadratureInfo(algebric_accuracy);
        std::vector<double> jacobian = the_dgele->local_to_global_jacobian(quad_info.quadraturePoint());
        int n_quadrature_point = quad_info.n_quadraturePoint();
        std::vector<point_t> q_point = the_dgele->local_to_global(quad_info.quadraturePoint());
        std::vector<std::vector<double>> basis_value = nei.basis_function_value(q_point);
        std::vector<std::vector<std::vector<double>>> basis_gradient = nei.basis_function_gradient(q_point);
        std::vector<std::vector<double>> out_normal = unitOutNormal(q_point,nei,*the_dgele);
        const std::vector<int>& element_dof = nei.dof();
        int n_element_dof = element_dof.size();
        for (int l = 0;l < n_quadrature_point;l ++) {
            double g_value = (*g)(q_point[l]);
            double Jxw = quad_info.weight(l)*jacobian[l]*volume;
            double h = Jxw;
            for (int j = 0;j < n_element_dof;j ++) {
	            this->operator[](element_dof[j]) += Jxw*100*g_value*basis_value[j][l]/h;
                this->operator[](element_dof[j]) -= Jxw*g_value*innerProduct(basis_gradient[j][l],out_normal[l]);
            }
        }
    }
}

void RightVec::set_algebric_accuracy(int i)
{
    algebric_accuracy = i;
}
