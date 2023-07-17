/**
This code is copied from part of:

The Arc-flow Vector Packing Solver (VPSolver).
Copyright (C) 2013-2016, Filipe Brandao
Faculdade de Ciencias, Universidade do Porto
Porto, Portugal. All rights reserved. E-mail: <fdabrandao@dcc.fc.up.pt>.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/
#ifndef VPS_UTILS_HPP
#define VPS_UTILS_HPP

#include "instance.hpp"
#include "arcflow.hpp"

#include <gurobi_c++.h>
#include <map>

GRBModel generate_arcflow_model(GRBEnv &env, std::map<Arc, GRBVar> &va, const Instance &inst, const Arcflow &afg);

//int solve_arcflow_model(GRBModel &model, float &time, int &obj_val, int &best_bound);
int solve_arcflow_model(GRBModel &model, int &obj_val, int &best_bound);

int get_arcflow_solution(const Instance &inst, Arcflow &afg, std::map<Arc, GRBVar> &va);

#endif // VPS_UTILS_HPP
