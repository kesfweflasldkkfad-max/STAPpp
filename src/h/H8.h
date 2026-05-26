/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

//! H8 solid element class
class CH8 : public CElement
{
public:
	CH8();
	~CH8();

	virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);
	virtual void Write(COutputter& output);
	virtual void ElementStiffness(double* Matrix);
	virtual void ElementStress(double* stress, double* Displacement);
};