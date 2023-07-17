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

#include <vector>
#include <algorithm>
#include <cmath>
#include "vps_utils.hpp"
#include "arcflowsol.hpp"

using namespace std;

GRBModel generate_arcflow_model(GRBEnv &env, std::map<Arc, GRBVar> &va, const Instance &inst, const Arcflow &afg)
{
    char vtype = inst.vtype;
    
    GRBModel model = GRBModel(env);
    model.set(GRB_StringAttr_ModelName, "flow");

    vector<Arc> As(afg.A);
    sort(all(As));
    //map<Arc, GRBVar> va;
    int lastv = afg.Ts[0]-1;
    for (int i = 0; i < inst.nbtypes; i++) {
        lastv = min(lastv, afg.Ts[i]-1);
    }
    for (int i = 0; i < 3; i++) {
        for (const Arc &a : As) {
            if (i == 1 && a.u != afg.S) {
                continue;
            } else if (i == 2 && a.v <= lastv) {
                continue;
            } else if (i == 0 && (a.u == afg.S || a.v > lastv)) {
                continue;
            }

            if (a.label == afg.LOSS || inst.relax_domains) {
                va[a] = model.addVar(
                    0.0, inst.n, 0, vtype);
            } else {
                va[a] = model.addVar(
                    0.0, inst.items[a.label].demand, 0, vtype);
            }
        }
    }
    model.update();

    for (int i = 0; i < inst.nbtypes; i++) {
        GRBVar &feedback = va[Arc(afg.Ts[i], afg.S, afg.LOSS)];
        feedback.set(GRB_DoubleAttr_Obj, inst.Cs[i]);
        if (inst.Qs[i] >= 0) {
            feedback.set(GRB_DoubleAttr_UB, inst.Qs[i]);
        }
    }

    vector<vector<Arc>> Al(inst.nsizes);
    vector<vector<Arc>> in(afg.NV);
    vector<vector<Arc>> out(afg.NV);

    for (const Arc &a : As) {
        if (a.label != afg.LOSS) {
            Al[a.label].push_back(a);
        }
        out[a.u].push_back(a);
        in[a.v].push_back(a);
    }

    for (int i = 0; i < inst.m; i++) {
        GRBLinExpr lin = 0;
        for (int it = 0; it < inst.nsizes; it++) {
            if (inst.items[it].type == i) {
                for (const Arc &a : Al[it]) {
                    lin += va[a];
                }
            }
        }
        if (inst.ctypes[i] == '>' || inst.relax_domains) {
            model.addConstr(lin >= inst.demands[i]);
        } else {
            model.addConstr(lin == inst.demands[i]);
        }
    }

    for (int u = 0; u < afg.NV; u++) {
        GRBLinExpr lin = 0;
        for (const Arc &a : in[u]) {
            lin += va[a];
        }
        for (const Arc &a : out[u]) {
            lin -= va[a];
        }
        model.addConstr(lin == 0);
    }

    Al.clear();
    in.clear();
    out.clear();

    return model;
}


//int solve_arcflow_model(GRBModel &model, float &time, int &obj_val, int &best_bound)
int solve_arcflow_model(GRBModel &model, int &obj_val, int &best_bound)
{
    obj_val = -1;
    best_bound = -1;

    model.optimize();
    
    int status = model.get(GRB_IntAttr_Status);
    obj_val = model.get(GRB_DoubleAttr_ObjVal);
    best_bound = model.get(GRB_DoubleAttr_ObjBound);

    return status;
}

int get_arcflow_solution(const Instance &inst, Arcflow &afg, map<Arc, GRBVar> &va)
{
    map<Arc, int> flow;
    for (const auto &a : va) {
        double x = a.second.get(GRB_DoubleAttr_X);
        int rx = static_cast<int>(round(x));
        //assert(x - rx <= EPS);
        if (rx > 0) {
            int u = a.first.u;
            int v = a.first.v;
            int lbl = a.first.label;
            Arc a(u, v, lbl);
            flow[a] = rx;
        }
    }
    ArcflowSol solution(inst, flow, afg.S, afg.Ts, afg.LOSS);
    return solution.get_solution_value();
}
