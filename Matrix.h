#ifndef _AMATRIX_H_
#define _AMATRIX_H_
#include <lac/sparse_matrix.h>
#include <lac/vector.h>
#include <vector>
#include <DGFEMSpace.h>

class AMatrix : public dealii::SparseMatrix<double>
{
    public:
    typedef typename dealii::SparseMatrix<double> spmat_t;
    typedef DGFEMSpace<double,2> dgspace_t;
    typedef typename dealii::Vector<double> vec_t;
    typedef typename dgspace_t::dg_element_t dg_element_t;
    typedef typename dgspace_t::element_t element_t;
    typedef typename dgspace_t::template_t template_t;
    typedef typename dgspace_t::dg_template_t dg_template_t;
    typedef typename dealii::SparsityPattern sparsitypattern;
    typedef typename AFEPack::Point<2> point_t;
    private:
    dgspace_t * space;
    sparsitypattern sp;
    int algebric_accuracy;
    double lambda = 1.0;
    public:
    AMatrix() = default;
    AMatrix(dgspace_t& sp):space(&sp){};
    ~AMatrix() = default;
    public:
    void build_sparsity_pattern();
    void build_element_matrix();
    void build_dgelement_matrix();
    void set_algebric_accuracy(int i);
    
};
#endif
