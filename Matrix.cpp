#include "Matrix.h"
#include <vector>
#include <Geometry.h>
#include <iostream>

void AMatrix::build_sparsity_pattern()
{
    auto the_dg_ele = space->beginDGElement();
    const auto end_dg_ele = space->endDGElement();
    int ndof = space->n_dof();
    int max_dof_couple = 0;
    for(;the_dg_ele != end_dg_ele; the_dg_ele++)
    {
        element_t& nei0 = the_dg_ele->neighbourElement(0);
        int n = nei0.n_dof();
        if(2*n>max_dof_couple) max_dof_couple = 2*n;
        if(the_dg_ele->p_neighbourElement(1)==NULL) continue;
        element_t& nei1 = the_dg_ele->neighbourElement(1);
        n += nei1.n_dof();
        if(2*n>max_dof_couple) max_dof_couple = 2*n;

    }
    //std::cout<<max_dof_couple<<std::endl;
    sp.reinit(ndof, ndof, std::min(max_dof_couple, ndof));
    
    auto the_ele = space->beginElement();
    const auto end_ele = space->endElement();
    for(; the_ele != end_ele; the_ele++)
    {
        std::vector<int>& dof0 = the_ele->dof();
        for(int i = 0; i<dof0.size(); i++)
            for(int j = 0; j<dof0.size(); j++)
                sp.add(dof0[i],dof0[j]);
    }
    auto the_dg_ele0 = space->beginDGElement();
    const auto end_dg_ele0 = space->endDGElement();
    for(; the_dg_ele0 != end_dg_ele0; the_dg_ele0++)
    {
        if(the_dg_ele0->p_neighbourElement(1)==NULL) continue;
        element_t& ele0 = the_dg_ele0->neighbourElement(0);
        element_t& ele1 = the_dg_ele0->neighbourElement(1);

        std::vector<int>& ele_dof0 = ele0.dof();
        std::vector<int>& ele_dof1 = ele1.dof();
        for(int i = 0; i<ele_dof0.size(); i++)
            for(int j = 0; j<ele_dof1.size(); j++)
                {
                    sp.add(ele_dof0[i],ele_dof1[j]);
                    sp.add(ele_dof1[j],ele_dof0[i]);
                }
    }
    sp.compress();
}
void AMatrix::build_element_matrix()
{
    this->reinit(sp);
    auto the_ele = space->beginElement();
    const auto end_ele = space->endElement();
    for(; the_ele != end_ele; the_ele++)
    {
        std::vector<int>& the_ele_dof = the_ele->dof();
        double volume = the_ele->templateElement().volume();
        const QuadratureInfo<2>& quad_info = the_ele->findQuadratureInfo(algebric_accuracy);
        const std::vector<point_t> q_point = the_ele->local_to_global(quad_info.quadraturePoint());
        int n_quadrature_point = quad_info.n_quadraturePoint();
        const std::vector<double> jacobian = the_ele->local_to_global_jacobian(quad_info.quadraturePoint());
        std::vector<std::vector<std::vector<double>>> basis_gradient = the_ele->basis_function_gradient(q_point);
        for(int l = 0; l < n_quadrature_point; l++)
        {
            double Jxw = quad_info.weight(l)*jacobian[l]*volume;
            for(int j = 0; j<the_ele_dof.size(); j++){
                for(int k = 0; k<the_ele_dof.size(); k++){
                    double val = Jxw*innerProduct(basis_gradient[j][l],basis_gradient[k][l]);
                    this->add(the_ele_dof[j],the_ele_dof[k],val);
                }
            }
        }
    }

}

void AMatrix::build_dgelement_matrix()
{
    auto the_dgele = space->beginDGElement();
    const auto end_dgele = space->endDGElement();
    for(; the_dgele != end_dgele; the_dgele++)
    {
        element_t& nei0 = the_dgele->neighbourElement(0);
        std::vector<int> nei_dof = nei0.dof();
        int temp_n = nei_dof.size();
        double volume = the_dgele->templateElement().volume();
        const QuadratureInfo<1>& quad_info = the_dgele->findQuadratureInfo(algebric_accuracy);
        std::vector<point_t> q_point = the_dgele->local_to_global(quad_info.quadraturePoint());
        std::vector<double> jacobian = the_dgele->local_to_global_jacobian(quad_info.quadraturePoint());
        int n_quadrature_point = q_point.size();
        std::vector<std::vector<double>> unit_normal0 = unitOutNormal(q_point,nei0,*the_dgele);
        std::vector<std::vector<double>> basis_value = nei0.basis_function_value(q_point);
        std::vector<std::vector<std::vector<double>>> basis_gradient = nei0.basis_function_gradient(q_point);
        double average = 1;
        std::vector<std::vector<std::vector<double>>> unit_normal;
        unit_normal.push_back(unit_normal0); 
        if(the_dgele->p_neighbourElement(1) != NULL)
        {
            element_t& nei1 = the_dgele->neighbourElement(1);
            std::vector<int>& nei_dof1 = nei1.dof();
            std::vector<std::vector<double>> unit_normal1 = unitOutNormal(q_point,nei1,*the_dgele);
            std::vector<std::vector<double>> basis_value1 = nei1.basis_function_value(q_point);
            std::vector<std::vector<std::vector<double>>> basis_gradient1 = nei1.basis_function_gradient(q_point);
            nei_dof.insert(nei_dof.end(),nei_dof1.begin(),nei_dof1.end());
            unit_normal.push_back(unit_normal1);
            //unit_normal.insert(unit_normal.end(),unit_normal1.begin(),unit_normal1.end());
            basis_value.insert(basis_value.end(),basis_value1.begin(),basis_value1.end());
            basis_gradient.insert(basis_gradient.end(),basis_gradient1.begin(),basis_gradient1.end());
            average = 0.5;
        }    
        for(int l = 0; l<n_quadrature_point; l++)
        {
            double Jxw = quad_info.weight(l)*jacobian[l]*volume;
            double h = Jxw;
            for(int j = 0; j<nei_dof.size(); j++){
                for(int k = 0; k<nei_dof.size(); k++){
                    int flag0 = 0,flag1 = 0;
                    if(j>temp_n-1) flag0 = 1; 
                    if(k>temp_n-1) flag1 = 1;

                    double val1 = Jxw*100*(basis_value[j][l]*basis_value[k][l])/h;
                    val1 *= innerProduct(unit_normal[flag0][l],unit_normal[flag1][l]);
                    this->add(nei_dof[j],nei_dof[k],val1);
                    double val2 = Jxw*average*basis_value[j][l]*innerProduct(basis_gradient[k][l],unit_normal[flag0][l]);
                    val2 = -val2;
                    this->add(nei_dof[j],nei_dof[k],val2);
                    this->add(nei_dof[k],nei_dof[j],val2);
                }
            }
        }
    }
}


void AMatrix::set_algebric_accuracy(int i)
{
    algebric_accuracy = i;
}


