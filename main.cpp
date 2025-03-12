#include <cmath>
#include <iostream>
#include <fstream>
#include <string.h>
#include <filesystem>
#include <vector>
#include <stdio.h>
#include <bits/stdc++.h>

using namespace std;

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;


vector<string> split(const string& str, const char del) {
    vector<string> tokens;
    stringstream ss(str);
    string word;
    while (!ss.eof()) {
        getline(ss, word, del);
        tokens.push_back(word);
    }
    return tokens;
}


struct EquationSolver {
    Matrix global_H;
    Matrix global_C;
    Vector global_P;
    Vector t0;
    Matrix HCdT;
    Vector PCdTt0;

    EquationSolver(int nN, int InitialTemp) {
        global_H = Matrix(nN, Vector(nN, 0));
        global_C = Matrix(nN, Vector(nN, 0));
        global_P = Vector(nN, 0);
        HCdT = Matrix(nN, Vector(nN, 0));
        PCdTt0 = Vector(nN, 0);
        t0 = Vector(nN, InitialTemp);
    }

    void set_to(Vector t0) {
        this->t0 = t0;
    }

    void print_global_H() {
        cout << endl << "GLOBAL H: " << endl;
        for (int i = 0; i < global_H.size(); i++) {
            for (int j = 0; j < global_H[i].size(); j++) {
                cout << global_H[i][j] << " ";
            }
            cout << endl;
        }
    }

    void print_global_C() {
        cout << endl << "GLOBAL C: " << endl;
        for (int i = 0; i < global_C.size(); i++) {
            for (int j = 0; j < global_C[i].size(); j++) {
                cout << global_C[i][j] << " ";
            }
            cout << endl;
        }
    }

    void print_HCdT() {
        cout << endl << "GLOBAL HCdT: " << endl;
        for (int i = 0; i < HCdT.size(); i++) {
            for (int j = 0; j < HCdT[i].size(); j++) {
                cout << HCdT[i][j] << " ";
            }
            cout << endl;
        }
    }

    void print_PCdTt0() {
        cout << endl << "GLOBAL PCdTt0: " << endl;
        for (int i = 0; i < PCdTt0.size(); i++) {
            cout << PCdTt0[i] << " ";
        }
        cout << endl;
    }

    void print_global_P() {
        cout << endl << "GLOBAL P: " << endl;
        for (int i = 0; i < global_P.size(); i++) {
            cout << global_P[i] << " ";
        }
        cout << endl;
    }

    void print_t0() {
        cout << endl << "GLOBAL t0: " << endl;
        for (int i = 0; i < t0.size(); i++) {
            cout << t0[i] << " ";
        }
        cout << endl;
    }

    void clear_HCdT() {
        for (int i = 0; i < HCdT.size(); i++) {
            for (int j = 0; j < HCdT[i].size(); j++) {
                HCdT[i][j] = 0;
            }
            // cout << endl;
        }
    }

    void clear_PCdTt0() {
        for (int i = 0; i < PCdTt0.size(); i++) {
            PCdTt0[i] = 0;
        }
        // cout << endl;
    }

    Vector gaussSolve(const Matrix& A, const Vector& b) {
        int n = A.size();

        // Augmented matrix
        Matrix augmented(n, Vector(n + 1));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                augmented[i][j] = A[i][j];
            }
            augmented[i][n] = b[i];
        }

        // Forward elimination
        for (int i = 0; i < n; ++i) {
            // Pivoting
            int maxRow = i;
            for (int k = i + 1; k < n; ++k) {
                if (std::abs(augmented[k][i]) > std::abs(augmented[maxRow][i])) {
                    maxRow = k;
                }
            }
            std::swap(augmented[i], augmented[maxRow]);

            // Check for singular matrix
            if (std::abs(augmented[i][i]) < 1e-9) {
                throw std::runtime_error("Matrix is singular or nearly singular.");
            }

            // Eliminate column i
            for (int k = i + 1; k < n; ++k) {
                double factor = augmented[k][i] / augmented[i][i];
                for (int j = i; j <= n; ++j) {
                    augmented[k][j] -= factor * augmented[i][j];
                }
            }
        }

        // Back substitution
        Vector x(n);
        for (int i = n - 1; i >= 0; --i) {
            x[i] = augmented[i][n] / augmented[i][i];
            for (int k = i - 1; k >= 0; --k) {
                augmented[k][n] -= augmented[k][i] * x[i];
            }
        }

        return x;
    }
};

struct Surface {
    vector<vector<float>> N;

    Surface(vector<vector<float>> N) {
        this->N = N;
    };



    void print_N() {
        for (int i = 0; i < N.size(); i++) {
            for (int j = 0; j < N[i].size(); j++) {
                cout << N[i][j] << " ";
            }
            cout << endl;
        }
    }
};

struct ElemUniv {
    int npc{};
    int npc_s{};
    vector<vector<float>> dNdksi;
    vector<vector<float>> dNdeta;
    vector<vector<float>> N_C;
    vector<Surface> surface;
    ElemUniv(int npc, int npc_s) {
        vector<vector<double>> pc;
        vector<vector<double>> pc_s;

        //PUNKTY CAŁKOWANIA W ŚRODKU ELEMENTU
        cout << "npc: " << npc << endl;
        if (npc == 4) {
            pc = {{-1/sqrt(3), -1/sqrt(3)}, {1/sqrt(3), -1/sqrt(3)}, {1/sqrt(3), 1/sqrt(3)}, {-1/sqrt(3), 1/sqrt(3)},  };
        }
        if (npc == 9) {
            double x_pom = 3.0/5.0;
            pc = {
                {-sqrt(x_pom), -sqrt(x_pom)}, {0, -sqrt(x_pom)}, {sqrt(x_pom), -sqrt(x_pom)},
                {-sqrt(x_pom), 0}, {0, 0}, {sqrt(x_pom), 0},
                {-sqrt(x_pom), sqrt(x_pom)}, {0, sqrt(x_pom)}, {sqrt(x_pom), sqrt(x_pom)}
            };
        }
        if (npc == 16) {
            double a = sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0));
            double b = sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0));
            pc = {
                {-a, -a}, {-b, -a}, {b, -a}, {a, -a},
                {-a, -b}, {-b, -b}, {b, -b}, {a, -b},
                {-a, b}, {-b, b}, {b, b}, {a, b},
                {-a, a}, {-b, a}, {b, a}, {a, a}
            };
        }

        //PUNKTY CAŁKOWANIA NA KRAWĘDZIACH ELEMENTU
        cout << "npc_s: " << npc_s << endl;
        if (npc_s == 4) {
            pc_s = { {-1/sqrt(3), -1}, {1/sqrt(3), -1},
                        {1, -1/sqrt(3)}, {1, 1/sqrt(3)},
                        {1/sqrt(3), 1}, {-1/sqrt(3), 1},
                        {-1, 1/sqrt(3)}, {-1, -1/sqrt(3)}};
        }
        if (npc_s == 9) {
            double x_pom = 3.0/5.0;
            pc_s = {
                {-sqrt(x_pom), -1}, {0, -1}, {sqrt(x_pom), -1},
                {1, -sqrt(x_pom)}, {1, 0}, {1, sqrt(x_pom)},
                {sqrt(x_pom), 1}, {0, 1}, {-sqrt(x_pom), 1},
                {-1, sqrt(x_pom)}, {-1, 0}, {-1, -sqrt(x_pom)}
            };
        }
        if (npc_s == 16) {
            double a = sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0));
            double b = sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0));
            pc_s = {
                {-a, -1}, {-b, -1}, {b, -1}, {a, -1},
                {1, -a}, {1, -b}, {1, b}, {1, a},
                {a, 1}, {b, 1}, {-b, 1}, {-a, 1},
                {-1, a}, {-1, b}, {-1, -b}, {-1, -a}
            };
        }
        this->npc = npc;
        this->npc_s = npc_s;
        vector<vector<float>> dNdksi(npc, vector<float>(4, 0));
        this->dNdksi = dNdksi;
        vector<vector<float>> dNdeta(npc, vector<float>(4, 0));
        this->dNdeta = dNdeta;
        vector<vector<float>> N_C(npc, vector<float>(4, 0));
        this->N_C = N_C;


        // this->dNdksi[0][0] =

        // for (int i = 0; i < npc; i++) {
        //     cout << "pc [" << i << "]" << pc[i][0] << " " << pc[i][1] <<  endl;
        // }
        for (int j = 0; j < 4; j++) {
            vector<vector<float>> N(sqrt(npc_s), vector<float>(4, 0));
            // cout << "MEGA N" << endl;
            for (int i = 0; i < sqrt(npc_s); i++) {
                N[i][0] = 0.25 * (1 - pc_s[i+(j*sqrt(npc_s))][0]) * (1 - pc_s[i+(j*sqrt(npc_s))][1]);
                N[i][1] = 0.25 * (1 + pc_s[i+(j*sqrt(npc_s))][0]) * (1 - pc_s[i+(j*sqrt(npc_s))][1]);
                N[i][2] = 0.25 * (1 + pc_s[i+(j*sqrt(npc_s))][0]) * (1 + pc_s[i+(j*sqrt(npc_s))][1]);
                N[i][3] = 0.25 * (1 - pc_s[i+(j*sqrt(npc_s))][0]) * (1 + pc_s[i+(j*sqrt(npc_s))][1]);

                // cout << " " << N[i][0] << " " << N[i][1] << " " << N[i][2] << " " << N[i][3] << endl;

                // cout << "HERE" << endl;
                // cout << pc_s[i+(j*sqrt(npc_s))][0] << endl;
            }

            // cout << endl;

            Surface s = Surface(N);
            this->surface.push_back(s);
        }

        for (int i = 0; i < npc; i++) {
            // cout << "wazne: " << pc[i][1];
            this->dNdksi[i][0] = -0.25 * (1 - pc[i][1]);
            this->dNdksi[i][1] = 0.25 * (1 - pc[i][1]);
            this->dNdksi[i][2] = 0.25 * (1 + pc[i][1]);
            this->dNdksi[i][3] = -0.25 * (1 + pc[i][1]);

            // cout << this->dNdksi[i][0] << " " << this->dNdksi[i][1] << " " << this->dNdksi[i][2] << " " << this->dNdksi[i][3] << endl;
            // cout << "pc: " << -0.25 * (1 - pc[i][0]) << endl;
            this->dNdeta[i][0] = -0.25 * (1 - pc[i][0]);
            this->dNdeta[i][1] = -0.25 * (1 + pc[i][0]);
            this->dNdeta[i][2] = 0.25 * (1 + pc[i][0]);
            this->dNdeta[i][3] = 0.25 * (1 - pc[i][0]);

            this->N_C[i][0] = 0.25 * (1 - pc[i][0]) * (1 - pc[i][1]);
            this->N_C[i][1] = 0.25 * (1 + pc[i][0]) * (1 - pc[i][1]);
            this->N_C[i][2] = 0.25 * (1 + pc[i][0]) * (1 + pc[i][1]);
            this->N_C[i][3] = 0.25 * (1 - pc[i][0]) * (1 + pc[i][1]);

            // cout << this->dNdeta[i][0] << " " << this->dNdeta[i][1] << " " << this->dNdeta[i][2] << " " << this->dNdeta[i][3] << endl;
        }

        // for (int i = 0; i < npc; i++) {
        //     cout << this->dNdksi[i][0] << " " << this->dNdksi[i][1] << " " << this->dNdksi[i][2] << " " << this->dNdksi[i][3] << endl;
        // }



        // cout << endl;

        // for (int i = 0; i < npc; i++) {
        //     cout << this->dNdeta[i][0] << " " << this->dNdeta[i][1] << " " << this->dNdeta[i][2] << " " << this->dNdeta[i][3] << endl;
        // }
        // cout << "TUTAJ" << endl;
        // for (int i = 0; i < npc; i++) {
        //     cout << this->N_C[i][0] << " " << this->N_C[i][1] << " " << this->N_C[i][2] << " " << this->N_C[i][3] << endl;
        // }
        //
        // cout << endl;
    }

    ElemUniv();

    vector<vector<float>> getDNdksi() {
        return dNdksi;
    }

    vector<vector<float>> getDNdeta() {
        return dNdeta;
    }

    void print_surface() {
        for (int i = 0; i < 4; i++) {
            surface[i].print_N();
            cout << endl;
        }
    }

};

struct jakobian {
    vector<vector<float>> J;
    vector<vector<float>> J1;
    float detJ{};

    jakobian(vector<float> dNdksi, vector<float> dNdeta, vector<vector<float>> finiteElement) {
        float dksi{};
        float ddeta{};
        vector<vector<float>> J(2, vector<float>(2));
        this->J = J;

        vector<vector<float>> J1(2, vector<float>(2));
        this->J1 = J1;


        for (int i = 0; i < 2; i++) {
            dksi = 0;
            ddeta = 0;
            for (int j = 0; j < 4; j++) {
                // cout << dNdksi[j] << " * " << finiteElement[i][j] << " + " << endl;
                dksi += finiteElement[i][j] * dNdksi[j];
                ddeta += finiteElement[i][j] * dNdeta[j];
            }
            // cout << " = " << dksi << endl;
            this->J[0][i] = dksi;
            this->J[1][i] = ddeta;
        }

        this->J1[0][0] = this->J[1][1];
        this->J1[1][1] = this->J[0][0];
        this->J1[0][1] = -this->J[0][1];
        this->J1[1][0] = -this->J[1][0];

        detJ = this->J[0][0] * this->J[1][1] - this->J[0][1] * this->J[1][0];
    }

    void print() {
        cout << "Jakobian: " << endl;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                cout << J[i][j] << " ";
            }
            cout << endl;
        }
        cout << "Odwrocony Jakobian" << endl;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                cout << J1[i][j] << " ";
            }
            cout << endl;
        }
        cout << "detJ: " << detJ << endl;
    }
    vector<vector<float>> get_J() {
        return J;
    };
    vector<vector<float>> get_J1() {
        return J1;
    }

    float get_detJ() {
        return detJ;
    }


};

struct node {
    friend ostream& operator<<(ostream& str, const node& n) {
        str << "x: " << n.x << " y: " << n.y;
        return str;
    }
    float x{}, y{};
    bool BC;


    node(const float x, const float y, const bool BC = 0) : x(x), y(y) {
        this -> x = x;
        this -> y = y;
        this -> BC = BC;
    }

    void setBC(bool BC) {
        this->BC = BC;
    }

    node() = default;
};


struct element {
    friend ostream& operator<<(ostream& str, const element& e) {
        str << e.ID[0] << " " << e.ID[1] << " " << e.ID[2] << " " << e.ID[3];
        return str;
    }
    int npc{};
    int npc_s{};
    vector<jakobian> Jakobian{};
    vector<vector<float>> finiteElement;
    vector<vector<float>> H_matrix;

    int ID[4]{};
    element(const int x, const int y, const int z, const int w, const int npc, const int npc_s) {
        ID[0] = x;
        ID[1] = y;
        ID[2] = z;
        ID[3] = w;
        this->npc = npc;
        this->npc_s = npc_s;
        finiteElement = vector<vector<float>>(2);
        H_matrix = vector<vector<float>>(4, vector<float>(4));
        // this->Jakobian = Jakobian;
        // vector<jakobian> Jakobian(npc);
        // this->Jakobian = Jakobian;
    }

    void set_Jakobian(ElemUniv eu) {

        for (int i = 0; i < npc; i++) {
            Jakobian.push_back(jakobian(eu.getDNdksi()[i], eu.getDNdeta()[i], finiteElement));
        }
    }

    void print_Jakobian() {

        for (int i = 0; i < npc; i++) {
            this->Jakobian[i].print();
        }
    }

    void makeFiniteElement (vector<node> nodes) {
        for (int i = 0; i < 4; i++) {
            finiteElement[0].push_back(nodes[this->ID[i]-1].x);
            finiteElement[1].push_back(nodes[this->ID[i]-1].y);
        }
    }

    void print_finiteElement() {
        // cout << "LOOK HERE" << endl;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 4; j++) {
                cout << finiteElement[i][j] << " ";
            }
            cout << endl;
        }
    }

    void H(ElemUniv eu, int Conductivity, int Alpha, int Tot, int Density, int SpecificHeat, int SimulationStepTime,  vector<node> nodes,  EquationSolver* eq) {
        vector<vector<float>> dNdksi = eu.getDNdksi();
        vector<vector<float>> dNdeta = eu.getDNdeta();
        vector<vector<vector<float>>> Hpc1(npc, vector<vector<float>>(4, vector<float>(4)));
        vector<vector<vector<float>>> Cpc1(npc, vector<vector<float>>(4, vector<float>(4)));
        vector<vector<float>> Hbcp(4, vector<float>(4,0));
        // vector<vector<float>> Hpc2(4, vector<float>(4));
        // vector<vector<float>> Hpc3(4, vector<float>(4));
        // vector<vector<float>> Hpc4(4, vector<float>(4));
        Vector P(4, 0);
        vector<vector<float>> H(4, vector<float>(4));
        vector<vector<float>> C(4, vector<float>(4));
        vector<vector<float>> Hbc(4, vector<float>(4,0));
        vector<vector<float>> bc_pom1(4, vector<float>(4));
        vector<vector<float>> bc_pom2(4, vector<float>(4));
        vector<vector<float>> pom1(4, vector<float>(4));
        vector<vector<float>> pom2(4, vector<float>(4));
        vector<vector<float>> dNdx(npc, vector<float>(4));
        vector<vector<float>> dNdy(npc, vector<float>(4));
        // cout << "siema" << Jakobian.size()  << endl;
        // cout << "npc" << this->npc << endl;

        for (int i = 0; i < npc; i++) {
            for (int j = 0; j < 4; j++) {
                vector<vector<float>> reJak = {{Jakobian[i].get_J1()[0][0] * (1/Jakobian[i].get_detJ()), Jakobian[i].get_J1()[0][1] * (1/Jakobian[i].get_detJ())}, {Jakobian[i].get_J1()[1][0] * (1/Jakobian[i].get_detJ()), Jakobian[i].get_J1()[1][1] * (1/Jakobian[i].get_detJ())}};
                // cout << reJak[0][0] << " " << reJak[0][1] << endl;
                // cout << reJak[1][0] << " " << reJak[1][1] << endl;

                dNdx[i][j] = reJak[0][0] * dNdksi[i][j] + reJak[0][1] * dNdeta[i][j];
                dNdy[i][j] = reJak[1][0] * dNdksi[i][j] + reJak[1][1] * dNdeta[i][j];
                // cout << "DNDX: " << dNdx[i][j] << " = " << Jakobian[i].get_J1()[0][0] << " * " << (1/Jakobian[i].get_detJ()) << " * " << dNdksi[i][j] << " + " << Jakobian[i].get_J1()[0][1] << " * " << (1/Jakobian[i].get_detJ()) << " * " << dNdeta[i][j] << endl;

            }
            // cout << endl;
        }
        // cout << "MOMENT" << endl;
        // // cout << npc << endl;
        // // cout << endl;
        // for (int i = 0; i < npc; i++) {
        //     for (int j = 0; j < 4; j++) {
        //         cout << dNdx[i][j] << " ";
        //     }
        //     cout << endl;
        // }
        // cout << endl;
        // for (int i = 0; i < npc; i++) {
        //     for (int j = 0; j < 4; j++) {
        //         cout << dNdy[i][j] << " ";
        //     }
        //     // cout << "super" << i;
        //     cout << endl;
        // }
        // cout << "Hpc1: " << Hpc1.size() << " " << Hpc1[0].size() << " " << Hpc1[0][0].size() << endl;
        for (int i = 0; i < npc; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    pom1[j][k] = dNdx[i][j] * dNdx[i][k];
                    pom2[j][k] = dNdy[i][j] * dNdy[i][k];
                    // cout << "POMOC: " <<
                    // cout << "suma: 30*(" << pom1[j][k] << " + " << pom2[j][k] << ")*" << Jakobian[i].get_detJ() << endl;
                    H[j][k] =
                    Hpc1[i][j][k] = Conductivity*(pom1[j][k] + pom2[j][k])*Jakobian[i].get_detJ();
                    // cout <<"MEGA: "<< Hpc1[i][j][k] << " = " << Conductivity << " * (" << pom1[j][k] << " + " << pom2[j][k] << ") * " << Jakobian[i].get_detJ() << endl;
                    Cpc1[i][j][k] = Density * SpecificHeat * (eu.N_C[i][j] * eu.N_C[i][k]) * Jakobian[i].get_detJ();
                    // cout << "suma: " << SpecificHeat << " * " << Density << " * (" << eu.N_C[i][j] * eu.N_C[i][k] << ")*" << Jakobian[i].get_detJ() <<  endl;
                }
            }




            // cout << "MOMENT2" << endl;
            //
            // for (int i = 0; i < 4; i++) {
            //     for (int j = 0; j < 4; j++) {
            //         cout << pom1[i][j] << " ";
            //     }
            //     cout << endl;
            // }
            //
            // cout << endl;
            //
            // for (int i = 0; i < 4; i++) {
            //     for (int j = 0; j < 4; j++) {
            //         cout << pom2[i][j] << " ";
            //     }
            //     cout << endl;
            // }
        }

        // cout << "MOMENT3" << endl;

        // for (int i = 0; i < 4; i++) {
        //     for (int j = 0; j < 4; j++) {
        //         cout << pom1[i][j] << " ";
        //     }
        //     cout << endl;
        // }
        //
        // cout << endl;
        //
        // for (int i = 0; i < 4; i++) {
        //     for (int j = 0; j < 4; j++) {
        //         cout << pom2[i][j] << " ";
        //     }
        //     cout << endl;
        // }


        // for (int i = 0; i < 4; i++) {
        //     for (int j = 0; j < 4; j++) {
        //         pom1[i][j] = dNdx[1][i] * dNdx[1][j];
        //         pom2[i][j] = dNdy[1][i] * dNdy[1][j];
        //         Hpc2[i][j] = 30*(pom1[i][j] + pom2[i][j])*Jakobian[i].get_detJ();//nie jestem pewny czy ten Jakobian
        //     }
        // }
        //
        // for (int i = 0; i < npc; i++) {
        //     for (int j = 0; j < 4; j++) {
        //         pom1[i][j] = dNdx[2][i] * dNdx[2][j];
        //         pom2[i][j] = dNdy[2][i] * dNdy[2][j];
        //         Hpc3[i][j] = 30*(pom1[i][j] + pom2[i][j])*Jakobian[i].get_detJ();//nie jestem pewny czy ten Jakobian
        //     }
        // }
        //
        // for (int i = 0; i < npc; i++) {
        //     for (int j = 0; j < 4; j++) {
        //         pom1[i][j] = dNdx[3][i] * dNdx[3][j];
        //         pom2[i][j] = dNdy[3][i] * dNdy[3][j];
        //         Hpc4[i][j] = 30*(pom1[i][j] + pom2[i][j])*Jakobian[i].get_detJ();//nie jestem pewny czy ten Jakobian
        //     }
        // }




        for(int i = 0; i < 4; i++) {
            float detJ = (sqrt(pow((nodes[ID[i]-1].x-nodes[ID[(i+1)%4]-1].x),2)+pow((nodes[ID[i]-1].y-nodes[ID[(i+1)%4]-1].y),2)))/2;
            // cout << "DETJ" << detJ << endl;
            if((nodes[ID[i]-1].BC) && (nodes[ID[(i+1)%4]-1].BC)) {

                // cout << "detJ: " << detJ << endl;
                for(int j =0; j < 4; j++) {
                    float w1 = 0;
                    float w2 = 0;
                    float w3 = 0;
                    float w4 = 0;
                    if(npc_s == 4) {
                        w1 = w2 = 1;
                        P[j] += Alpha * detJ * ((w1 *(eu.surface[i].N[0][j]*Tot)) + (w2 * (eu.surface[i].N[1][j]*Tot)));
                        (*eq).global_P[ID[j]-1] += Alpha * detJ * ((w1 *(eu.surface[i].N[0][j]*Tot)) + (w2 * (eu.surface[i].N[1][j]*Tot)));
                        (*eq).PCdTt0[ID[j]-1] += Alpha * detJ * ((w1 *(eu.surface[i].N[0][j]*Tot)) + (w2 * (eu.surface[i].N[1][j]*Tot)));

                        // cout << Alpha <<  " * " << detJ << " * (" << w1 << " *(" << eu.surface[i].N[0][j] << " * " << Tot << ") + " << w2 << " * (" << eu.surface[i].N[1][j] << " * " << Tot << "))" << " = " << (Alpha * detJ * ((w1 *(eu.surface[i].N[0][j]*Tot)) + (w2 * (eu.surface[i].N[1][j]*Tot)))) << endl;
                    }else if(npc_s == 9) {
                        w1 = 5.0/9.0;
                        w2 = 8.0/9.0;
                        w3 = 5.0/9.0;
                        float formula = Alpha * detJ * ((w1 *(eu.surface[i].N[0][j]*Tot)) + (w2 * (eu.surface[i].N[1][j]*Tot)) + (w3 * (eu.surface[i].N[2][j]*Tot)));
                        P[j] += formula;
                        (*eq).global_P[ID[j]-1] += formula;
                        (*eq).PCdTt0[ID[j]-1] += formula;

                    }else if(npc_s == 16) {
                        w1 = (18 - sqrt(30)) / 36; // waga dla punktu 1
                        w2 = (18 + sqrt(30)) / 36; // waga dla punktu 2
                        w3 = (18 + sqrt(30)) / 36; // waga dla punktu 3
                        w4 = (18 - sqrt(30)) / 36; // waga dla punktu 4
                        float formula = Alpha * detJ * ((w1 *(eu.surface[i].N[0][j]*Tot)) + (w2 * (eu.surface[i].N[1][j]*Tot)) + (w3 * (eu.surface[i].N[2][j]*Tot)) + (w4 * (eu.surface[i].N[3][j]*Tot)));
                        P[j] += formula;
                        (*eq).global_P[ID[j]-1] += formula;
                        (*eq).PCdTt0[ID[j]-1] += formula;

                    }

                    for(int k =0; k < 4; k++) {
                        w1 = 0;
                        w2 = 0;
                        w3 = 0;
                        w4 = 0;
                        if(npc_s == 4) {
                            w1 = w2 = 1;
                            Hbc[j][k] += ((w1 * Alpha * eu.surface[i].N[0][k] * eu.surface[i].N[0][j]) +(w2 * Alpha * eu.surface[i].N[1][k] * eu.surface[i].N[1][j])) * detJ;
                            // cout << "(("
                            //      << w1 << " * " << Alpha << " * "
                            //      << eu.surface[i].N[0][k] << " * " << eu.surface[i].N[0][j] << ") + ("
                            //      << w2 << " * " << Alpha << " * "
                            //      << eu.surface[i].N[1][k] << " * " << eu.surface[i].N[1][j] << ")) * "
                            //      << detJ
                            //      << std::endl;
                            (*eq).global_H[ID[j]-1][ID[k]-1] += ((Alpha * eu.surface[i].N[0][k] * eu.surface[i].N[0][j]) +(Alpha * eu.surface[i].N[1][k] * eu.surface[i].N[1][j])) * detJ;
                            (*eq).HCdT[ID[j]-1][ID[k]-1] += ((Alpha * eu.surface[i].N[0][k] * eu.surface[i].N[0][j]) +(Alpha * eu.surface[i].N[1][k] * eu.surface[i].N[1][j])) * detJ;

                            // cout << H[j][k] << " += " << Hbc[j][k] << endl;
                            H[j][k] += ((Alpha * eu.surface[i].N[0][k] * eu.surface[i].N[0][j]) +(Alpha * eu.surface[i].N[1][k] * eu.surface[i].N[1][j])) * detJ;
                            // cout << "((" <<  Alpha  << " * " << eu.surface[i].N[0][k] << " * " << eu.surface[i].N[0][j] << ") + " << "(" <<  Alpha  << " * " << eu.surface[i].N[1][k] << " * " << eu.surface[i].N[1][j] << ")) * " << detJ <<  endl;
                        }else if(npc_s == 9) {
                            w1 = 5.0/9.0;
                            w2 = 8.0/9.0;
                            w3 = 5.0/9.0;
                            float formula = ((w1 * Alpha * eu.surface[i].N[0][k] * eu.surface[i].N[0][j]) +(w2 * Alpha * eu.surface[i].N[1][k] * eu.surface[i].N[1][j]) + (w3 * Alpha * eu.surface[i].N[2][k] * eu.surface[i].N[2][j])) * detJ;
                            Hbc[j][k] += formula;
                            (*eq).global_H[ID[j]-1][ID[k]-1] += formula;
                            (*eq).HCdT[ID[j]-1][ID[k]-1] += formula;

                            // cout << H[j][k] << " += " << Hbc[j][k] << endl;
                            H[j][k] += formula;
                            // cout << "((" <<  Alpha  << " * " << eu.surface[i].N[0][k] << " * " << eu.surface[i].N[0][j] << ") + " << "(" <<  Alpha  << " * " << eu.surface[i].N[1][k] << " * " << eu.surface[i].N[1][j] << ")) * " << detJ <<  endl;
                        }else if(npc_s == 16){
                            w1 = (18 - sqrt(30)) / 36; // waga dla punktu 1
                            w2 = (18 + sqrt(30)) / 36; // waga dla punktu 2
                            w3 = (18 + sqrt(30)) / 36; // waga dla punktu 3
                            w4 = (18 - sqrt(30)) / 36; // waga dla punktu 4
                            float formula = ((w1 * Alpha * eu.surface[i].N[0][k] * eu.surface[i].N[0][j]) +(w2 * Alpha * eu.surface[i].N[1][k] * eu.surface[i].N[1][j]) + (w3 * Alpha * eu.surface[i].N[2][k] * eu.surface[i].N[2][j]) + (w4 * Alpha * eu.surface[i].N[3][k] * eu.surface[i].N[3][j])) * detJ;
                            Hbc[j][k] += formula;
                            (*eq).global_H[ID[j]-1][ID[k]-1] += formula;
                            (*eq).HCdT[ID[j]-1][ID[k]-1] += formula;

                            // cout << H[j][k] << " += " << Hbc[j][k] << endl;
                            H[j][k] += formula;
                            // cout << "((" <<  Alpha  << " * " << eu.surface[i].N[0][k] << " * " << eu.surface[i].N[0][j] << ") + " << "(" <<  Alpha  << " * " << eu.surface[i].N[1][k] << " * " << eu.surface[i].N[1][j] << ")) * " << detJ <<  endl;
                        }

                    }
                }
                // cout << endl;
            }
        }

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                float w1 = 0;
                float w2 = 0;
                float w3 = 0;
                float w4 = 0;
                if (npc == 4) {
                    w1 = 1;
                    w2 = 1;

                    float Hpom = Hpc1[0][i][j] * w1 * w1 + Hpc1[1][i][j] * w2 * w1 + Hpc1[2][i][j] * w1 * w2 + Hpc1[3][i][j] * w2 * w2;
                    float Cpom = Cpc1[0][i][j] * w1 * w1 + Cpc1[1][i][j] * w2 * w1 + Cpc1[2][i][j] * w1 * w2 + Cpc1[3][i][j] * w2 * w2;

                    H[i][j] = Hpom;
                    C[i][j] = Cpom;

                    (*eq).global_H[ID[i]-1][ID[j]-1] += Hpom;
                    (*eq).global_C[ID[i]-1][ID[j]-1] += Cpom;
                    (*eq).HCdT[ID[i]-1][ID[j]-1] += Hpom+(Cpom/SimulationStepTime);
                    (*eq).PCdTt0[ID[i]-1] += (Cpom/SimulationStepTime) * (*eq).t0[ID[j]-1];
                }else if (npc == 9) {
                    w1 = 5.0/9.0;
                    w2 = 8.0/9.0;
                    w3 = 5.0/9.0;

                    float Hpom = Hpc1[0][i][j] * w1 * w1
                            + Hpc1[1][i][j] * w1 * w2
                            + Hpc1[2][i][j] * w1 * w3
                            + Hpc1[3][i][j] * w2 * w1
                            + Hpc1[4][i][j] * w2 * w2
                            + Hpc1[5][i][j] * w2 * w3
                            + Hpc1[6][i][j] * w3 * w1
                            + Hpc1[7][i][j] * w3 * w2
                            + Hpc1[8][i][j] * w3 * w3;

                    float Cpom = Cpc1[0][i][j] * w1 * w1
                            + Cpc1[1][i][j] * w1 * w2
                            + Cpc1[2][i][j] * w1 * w3
                            + Cpc1[3][i][j] * w2 * w1
                            + Cpc1[4][i][j] * w2 * w2
                            + Cpc1[5][i][j] * w2 * w3
                            + Cpc1[6][i][j] * w3 * w1
                            + Cpc1[7][i][j] * w3 * w2
                            + Cpc1[8][i][j] * w3 * w3;



                    (*eq).global_H[ID[i]-1][ID[j]-1] += Hpom;
                    (*eq).global_C[ID[i]-1][ID[j]-1] += Cpom;
                    (*eq).HCdT[ID[i]-1][ID[j]-1] += Hpom+(Cpom/SimulationStepTime);
                    (*eq).PCdTt0[ID[i]-1] += (Cpom/SimulationStepTime) * (*eq).t0[ID[j]-1];


                    H[i][j] = Hpom;
                    C[i][j] = Cpom;

                }else if (npc == 16) {
                    w1 = (18 - sqrt(30)) / 36; // waga dla punktu 1
                    w2 = (18 + sqrt(30)) / 36; // waga dla punktu 2
                    w3 = (18 + sqrt(30)) / 36; // waga dla punktu 3
                    w4 = (18 - sqrt(30)) / 36; // waga dla punktu 4

                    float Hpom = Hpc1[0][i][j] * w1 * w1
                            + Hpc1[1][i][j] * w1 * w2
                            + Hpc1[2][i][j] * w1 * w3
                            + Hpc1[3][i][j] * w1 * w4
                            + Hpc1[4][i][j] * w2 * w1
                            + Hpc1[5][i][j] * w2 * w2
                            + Hpc1[6][i][j] * w2 * w3
                            + Hpc1[7][i][j] * w2 * w4
                            + Hpc1[8][i][j] * w3 * w1
                            + Hpc1[9][i][j] * w3 * w2
                            + Hpc1[10][i][j] * w3 * w3
                            + Hpc1[11][i][j] * w3 * w4
                            + Hpc1[12][i][j] * w4 * w1
                            + Hpc1[13][i][j] * w4 * w2
                            + Hpc1[14][i][j] * w4 * w3
                            + Hpc1[15][i][j] * w4 * w4;

                    float Cpom = Cpc1[0][i][j] * w1 * w1
                            + Cpc1[1][i][j] * w1 * w2
                            + Cpc1[2][i][j] * w1 * w3
                            + Cpc1[3][i][j] * w1 * w4
                            + Cpc1[4][i][j] * w2 * w1
                            + Cpc1[5][i][j] * w2 * w2
                            + Cpc1[6][i][j] * w2 * w3
                            + Cpc1[7][i][j] * w2 * w4
                            + Cpc1[8][i][j] * w3 * w1
                            + Cpc1[9][i][j] * w3 * w2
                            + Cpc1[10][i][j] * w3 * w3
                            + Cpc1[11][i][j] * w3 * w4
                            + Cpc1[12][i][j] * w4 * w1
                            + Cpc1[13][i][j] * w4 * w2
                            + Cpc1[14][i][j] * w4 * w3
                            + Cpc1[15][i][j] * w4 * w4;

                    (*eq).global_H[ID[i]-1][ID[j]-1] += Hpom;
                    (*eq).global_C[ID[i]-1][ID[j]-1] += Cpom;
                    (*eq).HCdT[ID[i]-1][ID[j]-1] += Hpom+(Cpom/SimulationStepTime);
                    (*eq).PCdTt0[ID[i]-1] += (Cpom/SimulationStepTime) * (*eq).t0[ID[j]-1];


                    H[i][j] = Hpom;
                    C[i][j] = Cpom;
                }
            }
        }

        // for (int i = 1; i < 16; i++) {
        //     for (int j = 1; j < 16; j++) {
        //         (*eq).PCdTt0[i] += ((*eq).global_C[i][j]/SimulationStepTime) * (*eq).t0[j];
        //     }
        // }

        // cout << "MEGA" << endl;
        // for (int i = 0; i < 4; i++) {
        //     for (int j = 0; j < 4; j++) {
        //         cout << Hbc[i][j] << " ";
        //     }
        //     cout << endl;
        // }

        // cout << "MEGA" << endl;
        // for (int i = 0; i < 4; i++) {
        //     for (int j = 0; j < 4; j++) {
        //         for (int k = 0; k < 4; k++) {
        //             cout << Hbcp[i][j][k] << " ";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        // cout << "MOMENT4" << endl;


        // for (int i = 0; i < npc; i++) {
        //     for (int j = 0; j < 4; j++) {
        //         cout << pom1[i][j] << " ";
        //     }
        //     cout << endl;
        // }
        //
        // cout << endl;
        //
        // for (int i = 0; i < npc; i++) {
        //     for (int j = 0; j < 4; j++) {
        //         cout << pom2[i][j] << " ";
        //     }
        //     cout << endl;
        // }
        //
        // cout << endl;
        // for (int i = 0; i < npc; i++) {
        //     for (int j = 0; j < 4; j++) {
        //         for (int k = 0; k < 4; k++) {
        //             cout << Hpc1[i][j][k] << " ";
        //         }
        //         cout << endl;
        //     }
        //     cout << endl;
        // }
        //
        //
        // cout << endl;

        // for (int i = 0; i < npc; i++) {
        //     for (int j = 0; j < 4; j++) {
        //         cout << Hpc2[i][j] << " ";
        //     }
        //     cout << endl;
        // }
        //
        // cout << endl;
        //
        // for (int i = 0; i < npc; i++) {
        //     for (int j = 0; j < 4; j++) {
        //         cout << Hpc3[i][j] << " ";
        //     }
        //     cout << endl;
        // }
        //
        // cout << endl;
        //
        // for (int i = 0; i < npc; i++) {
        //     for (int j = 0; j < 4; j++) {
        //         cout << Hpc4[i][j] << " ";
        //     }
        //     cout << endl;
        // }

        // cout << "MACIERZ C:" << endl;
        //
        // for (int i = 0; i < 4; i++) {
        //     for (int j = 0; j < 4; j++) {
        //         // this->H_matrix[i][j] = H[i][j];
        //         cout << C[i][j] << " ";
        //     }
        //     cout << endl;
        // }
        //
        // cout << "MACIERZ H:" << endl;
        //
        // for (int i = 0; i < 4; i++) {
        //     for (int j = 0; j < 4; j++) {
        //         // this->H_matrix[i][j] = H[i][j];
        //         cout << H[i][j] << " ";
        //     }
        //     cout << endl;
        // }
        // cout << "BC: ";
        // for (int i = 0; i < 4; i++) {
        //     cout << P[i] << " ";
        // }
        // cout << endl;
    }
    element() = default;
};

struct grid {
    int nN{}; //Liczba węzłów
    int nE{}; //Liczba elementów
    int npc{};
    int npc_s{};
    int Conductivity;
    int Alpha;
    int Tot;
    int Density;
    int SpecificHeat;
    int SimulationTime;
    int SimulationStepTime;
    vector<node> nodes;
    vector<element> elements;
    vector<int> BC;



    grid(const string& fileName_, const int nN, const int nE, int Conductivity, int Alpha, int Tot, int Density, int SpecificHeat, int SimulationTime, int SimulationStepTime, int npc, int npc_s) {
        this->npc = npc;
        this->nN = nN;
        this->nE = nE;
        this->Conductivity = Conductivity;
        this->Alpha = Alpha;
        this->npc_s = npc_s;
        this->Tot = Tot;
        this->Density = Density;
        this-> SpecificHeat = SpecificHeat;
        this->SimulationTime = SimulationTime;
        this->SimulationStepTime = SimulationStepTime;
        ifstream MyReadFile("..\\" + fileName_);
        string MyText;
        int i = -1;
        while (getline(MyReadFile, MyText)) {
            i++;
            if (i < 10) {
                continue;
            }
            if (i > 10 && i < 11 + nN) {
                vector<string> tokens = split(MyText, ',');
                nodes.push_back(node(stof(tokens[1]), stof(tokens[2]), 0));
            }
            if (i > 11 + nN && i < 12 + nN + nE) {
                vector<string> tokens = split(MyText, ',');
                elements.push_back(element(stoi(tokens[1]), stoi(tokens[2]), stoi(tokens[3]), stoi(tokens[4]), this->npc, this->npc_s));
            }

            if (i == 12 + nN + nE + 1 ) {
                // cout << "TUTAJ" << MyText << endl;
                vector<string> tokens = split(MyText, ',');
                for (int i = 0; i < tokens.size(); i++) {
                    BC.push_back(stoi(tokens[i]));
                }
                // elements.push_back(element(stoi(tokens[1]), stoi(tokens[2]), stoi(tokens[3]), stoi(tokens[4]), this->npc));
            }


        }


    }

    void setLoop(ElemUniv eu, EquationSolver* eq) {

        for (int i = 0; i < this->SimulationTime; i+=this->SimulationStepTime) {
            (*eq).clear_HCdT();
            (*eq).clear_PCdTt0();
            for (int j = 0; j < nE; j++) {
                elements[j].makeFiniteElement(this->nodes);
                elements[j].set_Jakobian(eu);
                // elements[j].print_Jakobian();
                elements[j].H(eu, this->Conductivity, this->Alpha, this->Tot, this->Density, this->SpecificHeat, this->SimulationStepTime, this->nodes, eq);
            }
            cout << endl << endl << "THIS IS " << i+this->SimulationStepTime << " ITERATION" ;
            // (*eq).print_global_H();
            // (*eq).print_global_P();
            // (*eq).print_global_C();
            // (*eq).print_HCdT();
            // (*eq).print_PCdTt0();
            Vector solve = (*eq).gaussSolve((*eq).HCdT, (*eq).PCdTt0);
            (*eq).set_to(solve);

            cout << endl << "SOLVED: " << endl;
            for (int i = 0; i < solve.size(); i++) {
                cout << solve[i] << " ";
            }
        }

    }

    void setBC() {
        for (int i = 0; i < nodes.size(); i++) {
            // cout << (find(BC.begin(), BC.end(), nodes[i].BC) != BC.end()) << endl;
            // cout << nodes[i].BC << endl;
            if (find(BC.begin(), BC.end(), i+1) != BC.end()) {
                // cout << "cokolwiek" << endl;
                nodes[i].BC = 1;
            }
        }
    }

    void printNodes() {
        for (int i = 0; i < nN; i++) {
            cout << nodes[i] << " BC_flag: " << nodes[i].BC << endl;
        }
    }

    void printElements() {
        for (int i = 0; i < nE; i++) {
            cout << elements[i] << endl;
        }
    }
};

struct GlobalData {
    int SimulationTime;
    int SimulationStepTime{};
    int Conductivity{};
    int Alfa{};
    int Tot{};
    int InitialTemp{};
    int Density{};
    int SpecificHeat{};
    int nN{}; //Liczba węzłów
    int nE{}; //Liczba elementów
    explicit GlobalData(const string& fileName_) {
        int variables[10];
        ifstream MyReadFile("..\\" + fileName_);
        string MyText;
        int i = 0;
        while (getline(MyReadFile, MyText) && i < 10) {
            vector<string> tokens = split(MyText, ' ');
            variables[i] = stoi(tokens.at(tokens.size()-1));
            i++;
        }
        this->SimulationTime = variables[0];
        this->SimulationStepTime = variables[1];
        this->Conductivity = variables[2];
        this->Alfa = variables[3];
        this->Tot = variables[4];
        this->InitialTemp = variables[5];
        this->Density = variables[6];
        this->SpecificHeat = variables[7];
        this->nN = variables[8];
        this->nE = variables[9];


    }

    void PrintData() const {
        cout << "Simulation Time: " << SimulationTime << endl;
        cout << "Simulation Step Time: " << SimulationStepTime << endl;
        cout << "Conductivity: " << Conductivity << endl;
        cout << "Alfa: " << Alfa << endl;
        cout << "Tot: " << Tot << endl;
        cout << "SpecificHeat: " << SpecificHeat << endl;
        cout << "NN: " << nN << endl;
        cout << "NE: " << nE << endl;
    }

    int getnN() const {
        return nN;
    }
    int getnE() const {
        return nE;
    }
};

int main() {

    int npc = 16;
    int npc_s = 9;
    ElemUniv eu(npc, npc_s);
    // eu.print_surface();

    GlobalData data_1("Test2_4_4_MixGrid.txt");
    // GlobalData data_1("Test1_4_4.txt");
    // GlobalData data_1("Test3_31_31_kwadrat.txt");
    EquationSolver eq(data_1.nN, data_1.InitialTemp);
    cout << &eq;
    data_1.PrintData();
    grid myGrid("Test2_4_4_MixGrid.txt", data_1.getnN(), data_1.getnE(), data_1.Conductivity, data_1.Alfa, data_1.Tot, data_1.Density, data_1.SpecificHeat, data_1.SimulationTime, data_1.SimulationStepTime, npc, npc_s);
    // grid myGrid("Test1_4_4.txt", data_1.getnN(), data_1.getnE(), data_1.Conductivity, data_1.Alfa, data_1.Tot, data_1.Density, data_1.SpecificHeat, data_1.SimulationTime, data_1.SimulationStepTime, npc, npc_s);
    // grid myGrid("Test3_31_31_kwadrat.txt", data_1.getnN(), data_1.getnE(), data_1.Conductivity, data_1.Alfa, data_1.Tot, data_1.Density, data_1.SpecificHeat, data_1.SimulationTime, data_1.SimulationStepTime, npc, npc_s);

    myGrid.setBC();
    myGrid.printNodes();
    myGrid.printElements();
    myGrid.setLoop(eu, &eq);



    // eq.print_t0();




    // vector<vector<float>> dNdksi = eu.getDNdksi();
    // vector<vector<float>> dNdeta = eu.getDNdeta();
    // vector<vector<float>> finiteElement = {{0, 0.025, 0.025, 0}, {0, 0, 0.025, 0.025}};
    // vector<jakobian> Jakobian;
    //
    // for (int i = 0; i < npc; i++) {
    //     Jakobian.push_back(jakobian(dNdksi[i], dNdeta[i], finiteElement));
    // }
    //
    // for (int i = 0; i < npc; i++) {
    //     Jakobian[i].print();
    // }

    // jakobian jk(dNdksi[0], dNdeta[0], finiteElement);
    // jakobian jk2(dNdksi[1], dNdeta[1], finiteElement);
    // jakobian jk3(dNdksi[2], dNdeta[2], finiteElement);
    // jakobian jk4(dNdksi[3], dNdeta[3], finiteElement);
    // Jakobian = {jk, jk2, jk3, jk4};
    // jk.print();
    // jk2.print();
    // jk3.print();
    // jk4.print();

    // element e(0, 0.025, 0.025, 0, npc);
    // e.set_Jakobian(Jakobian);
    //
    // e.H(eu, k);

    return 0;
}
