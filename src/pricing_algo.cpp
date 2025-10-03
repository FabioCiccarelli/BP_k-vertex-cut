#include "pricing_algo.h"
#include "gurobi_c++.h"
#include <iostream>
#include <vector>
#include <string>
using namespace std;

GRBModel* init_LP_pricer(int nedges, int nnodes, int* tail, int* head) {
    try {
        GRBEnv* env = new GRBEnv(true);
        env->start();

        GRBModel* model = new GRBModel(*env);
        model->set(GRB_IntAttr_ModelSense, -1); // Maximization

        // Suppress output
        model->set(GRB_IntParam_OutputFlag, 0);

        // Set number of threads to 1 for reproducibility
        model->set(GRB_IntParam_Threads, 1);

        // Create variables x_v for nodes
        vector<GRBVar> x_vars(nnodes);
        for (int v = 0; v < nnodes; ++v) {
            x_vars[v] = model->addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "x_" + to_string(v));
        }

        cout << "Variables x_v initialized successfully" << endl;

        // Create variables y_e for edges
        vector<GRBVar> y_vars;
        for (int e = 0; e < nedges; ++e) {
            y_vars.push_back(model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "y_" + to_string(e)));
        }

        cout << "Variables y_e initialized successfully" << endl;

        // Add constraints: y_e <= x_u and y_e <= x_v for each edge e=(u,v)
        for (int e = 0; e < nedges; ++e) {
            int u = tail[e] - 1;
            int v = head[e] - 1;
            model->addConstr(y_vars[e] <= x_vars[u], "y_le_xu_" + to_string(e));
            model->addConstr(y_vars[e] <= x_vars[v], "y_le_xv_" + to_string(e));
        }

        cout << "Constraints added successfully" << endl;

        model->update();

        cout << "Gurobi model initialized successfully" << endl;

        return model;
    } catch (GRBException &e) {
        cerr << "Gurobi error: " << e.getMessage() << endl;
        return nullptr;
    }

}


void set_coeffs(GRBModel* model, double pi_val, double* edge_coeffs, double* node_coeffs, int nedges, int nnodes) {
    for (int e = 0; e < nedges; ++e) {
        GRBVar var = model->getVarByName("y_" + to_string(e));
        var.set(GRB_DoubleAttr_Obj, edge_coeffs[e]);
    }
    for (int v = 0; v < nnodes; ++v) {
        GRBVar var = model->getVarByName("x_" + to_string(v));
        var.set(GRB_DoubleAttr_Obj, node_coeffs[v] + pi_val);
    }
    model->update();
}


PoolResult solve_pricing(GRBModel* model, int nnodes, int nedges, double pi_val) {
    PoolResult result;
    result.obj_val = -1;
    result.num_solutions = 0;

    try {

        // write the model in a LP file for debugging
        model->optimize();

        double obj = model->get(GRB_DoubleAttr_ObjVal);

        if (obj > 1e-6) { // Positive objective
            result.num_solutions = 1;
            result.obj_val = obj;
            vector<double> sol;
            for (int v = 0; v < nnodes; ++v) {
                GRBVar var = model->getVarByName("x_" + to_string(v));
                sol.push_back(var.get(GRB_DoubleAttr_X));
            }
            result.solutions.push_back(sol);
            return result;
        }

        // Objective is zero, try fixing x variables
        for (int v = 0; v < nnodes; ++v) {
            GRBVar xvar = model->getVarByName("x_" + to_string(v));
            double old_obj = xvar.get(GRB_DoubleAttr_Obj);
            xvar.set(GRB_DoubleAttr_Obj, old_obj - pi_val);
            model->update();
            model->optimize();

            double new_obj = model->get(GRB_DoubleAttr_ObjVal);
            if (new_obj > 1e-9) {
                result.num_solutions = 1;
                result.obj_val = new_obj;
                vector<double> sol;
                for (int u = 0; u < nnodes; ++u) {
                    GRBVar var = model->getVarByName("x_" + to_string(u));
                    sol.push_back(var.get(GRB_DoubleAttr_X));
                }
                result.solutions.push_back(sol);
                xvar.set(GRB_DoubleAttr_Obj, old_obj);
                model->update();
                return result;
            }
            xvar.set(GRB_DoubleAttr_Obj, old_obj);
            model->update();
        }

        // No positive solution found
        result.obj_val = -1;
        result.num_solutions = 0;
        result.solutions.clear();
        return result;

    } catch (GRBException &e) {
        cerr << "Gurobi error: " << e.getMessage() << endl;
        result.obj_val = -1;
        result.num_solutions = 0;
        result.solutions.clear();
        return result;
    }
}


void free_LP_pricer(GRBModel* model) {
    if (model) {
        delete model;
    }
}