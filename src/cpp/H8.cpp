/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "H8.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

namespace
{
    struct GaussPoint
    {
        double xi;
        double eta;
        double zeta;
    };

    const double H8Gauss = 1.0 / sqrt(3.0);

    const double LocalCoordinates[8][3] =
    {
        {-1.0, -1.0, -1.0},
        { 1.0, -1.0, -1.0},
        { 1.0,  1.0, -1.0},
        {-1.0,  1.0, -1.0},
        {-1.0, -1.0,  1.0},
        { 1.0, -1.0,  1.0},
        { 1.0,  1.0,  1.0},
        {-1.0,  1.0,  1.0}
    };

    const GaussPoint GaussPoints[8] =
    {
        {-H8Gauss, -H8Gauss, -H8Gauss},
        { H8Gauss, -H8Gauss, -H8Gauss},
        { H8Gauss,  H8Gauss, -H8Gauss},
        {-H8Gauss,  H8Gauss, -H8Gauss},
        {-H8Gauss, -H8Gauss,  H8Gauss},
        { H8Gauss, -H8Gauss,  H8Gauss},
        { H8Gauss,  H8Gauss,  H8Gauss},
        {-H8Gauss,  H8Gauss,  H8Gauss}
    };

    void ComputeNaturalGradients(double xi, double eta, double zeta, double gradients[8][3])
    {
        for (unsigned int node = 0; node < 8; node++)
        {
            double xi_i = LocalCoordinates[node][0];
            double eta_i = LocalCoordinates[node][1];
            double zeta_i = LocalCoordinates[node][2];

            gradients[node][0] = 0.125 * xi_i * (1.0 + eta * eta_i) * (1.0 + zeta * zeta_i);
            gradients[node][1] = 0.125 * eta_i * (1.0 + xi * xi_i) * (1.0 + zeta * zeta_i);
            gradients[node][2] = 0.125 * zeta_i * (1.0 + xi * xi_i) * (1.0 + eta * eta_i);
        }
    }

    void BuildConstitutiveMatrix(const CH8Material& material, double constitutive[6][6])
    {
        double lambda = material.E * material.Nu / ((1.0 + material.Nu) * (1.0 - 2.0 * material.Nu));
        double mu = material.E / (2.0 * (1.0 + material.Nu));

        for (unsigned int i = 0; i < 6; i++)
            for (unsigned int j = 0; j < 6; j++)
                constitutive[i][j] = 0.0;

        for (unsigned int i = 0; i < 3; i++)
            for (unsigned int j = 0; j < 3; j++)
                constitutive[i][j] = lambda;

        constitutive[0][0] += 2.0 * mu;
        constitutive[1][1] += 2.0 * mu;
        constitutive[2][2] += 2.0 * mu;
        constitutive[3][3] = mu;
        constitutive[4][4] = mu;
        constitutive[5][5] = mu;
    }

    bool Invert3x3(const double matrix[3][3], double inverse[3][3], double& determinant)
    {
        determinant = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
                    - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
                    + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);

        if (determinant <= 0.0)
            return false;

        double inverse_det = 1.0 / determinant;

        inverse[0][0] = inverse_det * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]);
        inverse[0][1] = inverse_det * (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]);
        inverse[0][2] = inverse_det * (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]);
        inverse[1][0] = inverse_det * (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]);
        inverse[1][1] = inverse_det * (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]);
        inverse[1][2] = inverse_det * (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]);
        inverse[2][0] = inverse_det * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
        inverse[2][1] = inverse_det * (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]);
        inverse[2][2] = inverse_det * (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
        return true;
    }

    bool BuildBMatrix(const CNode* const* nodes, double xi, double eta, double zeta, double B[6][24], double& detJ)
    {
        double natural_gradients[8][3];
        ComputeNaturalGradients(xi, eta, zeta, natural_gradients);

        double jacobian[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        for (unsigned int a = 0; a < 8; a++)
            for (unsigned int i = 0; i < 3; i++)
                for (unsigned int j = 0; j < 3; j++)
                    jacobian[i][j] += nodes[a]->XYZ[i] * natural_gradients[a][j];

        double inverse_jacobian[3][3];
        if (!Invert3x3(jacobian, inverse_jacobian, detJ))
            return false;

        double cartesian_gradients[8][3];
        for (unsigned int a = 0; a < 8; a++)
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                cartesian_gradients[a][i] = 0.0;
                for (unsigned int j = 0; j < 3; j++)
                    cartesian_gradients[a][i] += natural_gradients[a][j] * inverse_jacobian[j][i];
            }
        }

        for (unsigned int i = 0; i < 6; i++)
            for (unsigned int j = 0; j < 24; j++)
                B[i][j] = 0.0;

        for (unsigned int a = 0; a < 8; a++)
        {
            unsigned int col = 3 * a;
            double dndx = cartesian_gradients[a][0];
            double dndy = cartesian_gradients[a][1];
            double dndz = cartesian_gradients[a][2];

            B[0][col] = dndx;
            B[1][col + 1] = dndy;
            B[2][col + 2] = dndz;
            B[3][col] = dndy;
            B[3][col + 1] = dndx;
            B[4][col + 1] = dndz;
            B[4][col + 2] = dndy;
            B[5][col] = dndz;
            B[5][col + 2] = dndx;
        }

        return true;
    }
}

CH8::CH8()
{
    NEN_ = 8;
    nodes_ = new CNode*[NEN_];

    ND_ = 24;
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

CH8::~CH8()
{
}

bool CH8::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int node_numbers[8];
    unsigned int MSet;

    for (unsigned int i = 0; i < 8; i++)
        Input >> node_numbers[i];
    Input >> MSet;

    ElementMaterial_ = dynamic_cast<CH8Material*>(MaterialSets) + MSet - 1;
    for (unsigned int i = 0; i < 8; i++)
        nodes_[i] = &NodeList[node_numbers[i] - 1];

    return true;
}

void CH8::Write(COutputter& output)
{
    for (unsigned int i = 0; i < 8; i++)
        output << setw(9) << nodes_[i]->NodeNumber;
    output << setw(12) << ElementMaterial_->nset << endl;
}

void CH8::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    CH8Material* material = dynamic_cast<CH8Material*>(ElementMaterial_);
    double constitutive[6][6];
    BuildConstitutiveMatrix(*material, constitutive);

    double full[24][24];
    for (unsigned int i = 0; i < 24; i++)
        for (unsigned int j = 0; j < 24; j++)
            full[i][j] = 0.0;

    for (unsigned int gp = 0; gp < 8; gp++)
    {
        double B[6][24];
        double detJ;
        if (!BuildBMatrix(nodes_, GaussPoints[gp].xi, GaussPoints[gp].eta, GaussPoints[gp].zeta, B, detJ))
        {
            cerr << "*** Error *** H8 element has a non-positive Jacobian determinant." << endl;
            exit(6);
        }

        double DB[6][24];
        for (unsigned int i = 0; i < 6; i++)
            for (unsigned int j = 0; j < 24; j++)
            {
                DB[i][j] = 0.0;
                for (unsigned int k = 0; k < 6; k++)
                    DB[i][j] += constitutive[i][k] * B[k][j];
            }

        for (unsigned int i = 0; i < 24; i++)
            for (unsigned int j = 0; j < 24; j++)
                for (unsigned int k = 0; k < 6; k++)
                    full[i][j] += B[k][i] * DB[k][j] * detJ;
    }

    unsigned int index = 0;
    for (unsigned int col = 0; col < 24; col++)
        for (unsigned int row = col; row < 24; row++)
            Matrix[index++] = full[row][col];
}

void CH8::ElementStress(double* stress, double* Displacement)
{
    clear(stress, 6);

    CH8Material* material = dynamic_cast<CH8Material*>(ElementMaterial_);
    double constitutive[6][6];
    BuildConstitutiveMatrix(*material, constitutive);

    double element_displacement[24];
    for (unsigned int i = 0; i < 24; i++)
    {
        if (LocationMatrix_[i])
            element_displacement[i] = Displacement[LocationMatrix_[i] - 1];
        else
            element_displacement[i] = 0.0;
    }

    for (unsigned int gp = 0; gp < 8; gp++)
    {
        double B[6][24];
        double detJ;
        if (!BuildBMatrix(nodes_, GaussPoints[gp].xi, GaussPoints[gp].eta, GaussPoints[gp].zeta, B, detJ))
        {
            cerr << "*** Error *** H8 element has a non-positive Jacobian determinant." << endl;
            exit(6);
        }

        double strain[6];
        for (unsigned int i = 0; i < 6; i++)
        {
            strain[i] = 0.0;
            for (unsigned int j = 0; j < 24; j++)
                strain[i] += B[i][j] * element_displacement[j];
        }

        for (unsigned int i = 0; i < 6; i++)
            for (unsigned int j = 0; j < 6; j++)
                stress[i] += constitutive[i][j] * strain[j] / 8.0;
    }
}