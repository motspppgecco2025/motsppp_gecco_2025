#include "gurobi_c++.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include <string>
#include <tuple>
#include <set>
#include <cmath>
#include <chrono>

using namespace std;

int s_initial = 0;

// Classe de callback para limitar o tempo de execução
class TimeLimitCallback : public GRBCallback {
public:
    std::chrono::steady_clock::time_point start_time;
    int time_limit_seconds;
    bool stop_optimization = false;

    TimeLimitCallback(int time_limit) : time_limit_seconds(time_limit) {
        start_time = std::chrono::steady_clock::now();
    }

protected:
    void callback() override {
        if (where == GRB_CB_MIP || where == GRB_CB_MIPSOL || where == GRB_CB_MIPNODE) {
            auto current_time = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed_seconds = current_time - start_time;

            if (elapsed_seconds.count() >= time_limit_seconds) {
                stop_optimization = true;
                abort();  // Interrompe a otimização
            }
        }
    }
};

struct Passenger {
    int wk;  // Contribuição do passageiro
    int ok;  // Cidade de origem do passageiro
    int dk;  // Cidade de destino do passageiro
    int tk;  // Tempo máximo da rota do passageiro
};

struct InputData {
    int n;                  // Número de cidades
    int l;                  // Número de passageiros
    int vehicle_capacity;   // Capacidade do veículo
    vector<vector<int>> c;  // Matriz de custosnObjectives
    vector<vector<int>> y;  // Matriz de tempo
    vector<int> b;          // Bônus de cada cidade
    vector<int> g;          // Tempo de coleta de bônus de cada cidade
    vector<Passenger> passengers; // Dados dos passageiros
};

// Estrutura para armazenar os pontos não-dominados
struct Solution {
    double obj_cost;
    double obj_time;
    double obj_bonus;
    vector<int> tour;
    vector<int> collectedBonuses;
    vector<int> passengers;
};

// leitura dos dados do 'input.txt'
InputData readInput(const string &filename) {
    ifstream file(filename);
    InputData data;

    if (file.is_open()) {
        int quota;
        file >> data.n >> data.l >> data.vehicle_capacity;

        data.c.resize(data.n, vector<int>(data.n));
        data.y.resize(data.n, vector<int>(data.n));
        data.b.resize(data.n);
        data.g.resize(data.n);
        data.passengers.resize(data.l);

        for (int i = 0; i < data.n; i++) {
            for (int j = 0; j < data.n; j++) {
                file >> data.c[i][j];
            }
        }

        for (int i = 0; i < data.n; i++) {
            for (int j = 0; j < data.n; j++) {
                file >> data.y[i][j];
            }
        }

        for (int i = 0; i < data.l; i++) {
            int wk, ok, dk, tk;
            file >> wk >> ok >> dk >> tk;
            data.passengers[i] = {wk, ok, dk, tk};
        }

        file >> quota; // não utilizada;

        for (int i = 0; i < data.n; i++) {
            int cityIndex, bonus, collectTime;
            file >> cityIndex >> bonus >> collectTime;
            data.b[cityIndex] = bonus;
            data.g[cityIndex] = collectTime;
        }

        file.close();
    } else {
        cerr << "Não foi possível abrir o arquivo " << filename << endl;
    }

    return data;
}

// Função para extrair o tour da solução com base na variável x
vector<int> extractTour(GRBModel &model, vector<vector<GRBVar>> &x, int n) {
    vector<int> tour;

    try {
        vector<bool> visited(n, false);
        int current = 0;
        tour.push_back(current);
        visited[current] = true;

        while (tour.size() < n) {
            bool found = false;
            for (int next = 0; next < n; next++) {
                //if (!visited[next] && x[current][next].get(GRB_DoubleAttr_X) > .5) {
                if (!visited[next] && fabs(x[current][next].get(GRB_DoubleAttr_X)-1.0) < 0.000001){
                    tour.push_back(next);
                    visited[next] = true;
                    current = next;
                    found = true;
                    break;
                }
            }
            if (!found) {
                break; // Não há mais cidades para visitar
            }
        }

        // Verifica se o tour retornou ao ponto inicial
        if (current != 0) {
            tour.push_back(0); // Retorna à cidade inicial
        }

    } catch (const GRBException &e) {
        cerr << "Erro: " << e.getMessage() << endl;
        cerr << "Código do erro: " << e.getErrorCode() << endl;
    } catch (const exception &e) {
        cerr << "Exceção: " << e.what() << endl;
    } catch (...) {
        cerr << "Exceção desconhecida." << endl;
    }

    return tour;
};

// Função para coletar bônus ao longo do tour
vector<int> collectBonuses(const vector<int>& tour, const InputData &data, GRBModel &model,  vector<GRBVar> &p) {
    vector<int> collectedBonuses;

    try {
        for (int city : tour) {
            if (fabs(p[city].get(GRB_DoubleAttr_X)-1.0) < 0.000001) {
                collectedBonuses.push_back(city);
            }
        }
    } catch (const GRBException &e) {
        cerr << "Erro: " << e.getMessage() << endl;
        cerr << "Código do erro: " << e.getErrorCode() << endl;
    } catch (const exception &e) {
        cerr << "Exceção: " << e.what() << endl;
    } catch (...) {
        cerr << "Exceção desconhecida." << endl;
    }

    return collectedBonuses;
}

// Função para coletar os passageiros viáveis
vector<int> collectPassengers(GRBModel &model, vector<vector<vector<GRBVar>>> &h, int n, int m) {
    vector<int> passengers;

    try {
        for (int k = 0; k < m; k++) {
            bool collected = false;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (h[k][i][j].get(GRB_DoubleAttr_X) > 0.5) {
                        passengers.push_back(k);
                        collected = true;
                        break;
                    }
                }
                if (collected) break;
            }
        }
    } catch (const GRBException &e) {
        cerr << "Erro: " << e.getMessage() << endl;
        cerr << "Código do erro: " << e.getErrorCode() << endl;
    } catch (const exception &e) {
        cerr << "Exceção: " << e.what() << endl;
    } catch (...) {
        cerr << "Exceção desconhecida." << endl;
    }
    return passengers;
}

// Função para salvar os resultados intermediários
void saveProgress(const vector<Solution> &solutions, const string &path_output, const string &file_name_output) {
    ofstream outFile(path_output + file_name_output);
    if (outFile.is_open()) {
        outFile << "Fo(custo)\tFo(tempo)\tFo(bônus)\n";
        outFile << endl;

        for (const auto &solution : solutions) {
            outFile << fabs(solution.obj_cost) << "\t" << fabs(solution.obj_time) << "\t" << solution.obj_bonus << "\n";
            outFile << "tour: ";
            for (int city : solution.tour) {
                outFile << city << " ";
            }
            outFile << "\n";
            outFile << "bonus: ";
            for (int bonus : solution.collectedBonuses) {
                outFile << bonus << " ";
            }
            outFile << "\n";
            outFile << "passengers: ";
            for (int passenger : solution.passengers) {
                outFile << passenger << " ";
            }
            outFile << "\n\n";
        }
        outFile.close();
        cout << "Resultados intermediários gravados!" << endl;
    } else {
        cerr << "Não foi possível abrir o arquivo de resultados intermediários." << endl;
    }
}

// Função para configurar e resolver o modelo Gurobi
void solveModel(GRBEnv &env, const InputData &data, vector<Solution> &solutions, string fragment_file) {
    try {
        bool stop = false;
        int count_interaction = 0;
        int count_file = 1;

        int time_limit_seconds = 7200;  // 2 horas em segundos
        TimeLimitCallback cb(time_limit_seconds);
        
        auto start_time = std::chrono::steady_clock::now();
        auto last_save_time = start_time;  // Para controle de 30 minutos

        while(!stop) {
            GRBModel model = GRBModel(env);
            model.setCallback(&cb);
            model.set(GRB_DoubleParam_MIPGap, 1e-9);
            model.set(GRB_DoubleParam_MIPGapAbs, 1e-9);


            int n = data.n;
            int m = data.l;
            double M = 1000000;
            double epsilon = 0.001;
            cout << M << endl;
            cout << epsilon << endl;
            int nObjectives = 3;

            // +----------------------+
            // + Variáveis de decisão +
            // +----------------------+
            
            vector<vector<GRBVar>> x(n, vector<GRBVar>(n));                                      // Variável binária que indica se a aresta de cidade i para cidade j é percorrida;
            vector<GRBVar> u(n);                                                                 // Variável inteira que indica a ordem da cidade i na rota do caixeiro viajante, usada nas restrições MTZ para eliminação de subtours
            vector<GRBVar> p(n);                                                                 // Variável binária que indica se o bônus na cidade i é coletado
            vector<vector<vector<GRBVar>>> h(m, vector<vector<GRBVar>>(n, vector<GRBVar>(n)));   // Variável binária que indica se o k-ésimo passageiro está no veículo na aresta (i, j)
            vector<vector<vector<GRBVar>>> psi(n, vector<vector<GRBVar>>(n, vector<GRBVar>(m)));   // Variável continua para linearização das restrições envolvendo alpha // by Felipe 28/08/2024
            vector<vector<GRBVar>> alpha(n, vector<GRBVar>(n));                                  // Variável contínua que indica o inverso do número de ocupantes do carro na aresta (i, j);
            vector<vector<GRBVar>> y(solutions.size(), vector<GRBVar>(nObjectives));             // Ajuste o tamanho conforme necessário

            for (int i = 0; i < n; i++) {
                p[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "p_" + to_string(i));
                
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + to_string(i) + "_" + to_string(j));
                        alpha[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "alpha_" + to_string(i) + "_" + to_string(j));
                    }
                }
            }
            
            for (int i = 0; i < n; i++) {
                u[i] = model.addVar(1.0, n - 1, 0.0, GRB_INTEGER, "u_" + to_string(i));
            }

            for (int k = 0; k < m; k++){
                for (int i = 0; i < n; i++){
                    for (int j = 0; j < n; j++){
                        h[k][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "h_" + to_string(k) + "_" + to_string(i) + "_" + to_string(j));
                        psi[i][j][k] = model.addVar(0, 1, 0.0, GRB_CONTINUOUS, "psi" + std::to_string(i) + std::to_string(j) + std::to_string(k)); // by Felipe (28/08/2024)

                    }
                }
            }

            for (int t = 0; t < solutions.size(); t++) {
                for (int j = 0; j < nObjectives; j++) {
                    y[t][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_" + to_string(t) + "_" + to_string(j));
                }
            }


            // +-----------------+
            // + Funcao objetivo +
            // +-----------------+

            GRBQuadExpr obj_cost = 0;
            GRBLinExpr obj_time  = 0;
            GRBLinExpr obj_bonus = 0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        obj_cost -= x[i][j] * data.c[i][j] * alpha[i][j];
                        obj_time -= data.y[i][j] * x[i][j];
                    }
                }
                obj_time -= p[i] * data.g[i];
                obj_bonus += p[i] * data.b[i];
            }

            model.setObjective(obj_bonus + (epsilon * obj_cost) + (epsilon * obj_time), GRB_MAXIMIZE);

            

            // +------------------------------+
            // + Restrições Lokman e Koksalan +
            // +------------------------------+

            for (int t = 0; t < solutions.size(); t++){
                model.addQConstr(obj_cost >= ((solutions[t].obj_cost + 1) * y[t][0]) - (M * (1 - y[t][0])), "comp1_" + to_string(t) + "_" + to_string(0));
                model.addQConstr(obj_time >= ((solutions[t].obj_time + 1) * y[t][1]) - (M * (1 - y[t][1])), "comp1_" + to_string(t) + "_" + to_string(1));
            }

            for (int t = 0; t < solutions.size(); t++) {
                GRBLinExpr sumY = 0;
                for (int j = 0; j < nObjectives; j++){
                    if (j != 2) {
                        sumY += y[t][j];
                    }
                }
                model.addConstr(sumY == 1, "compSum_" + to_string(t));
            }

            // +------------------------+
            // + Restrições do problema +
            // +------------------------+

            // Restrição 4: Cada cidade, exceto a origem, deve ser visitada exatamente uma vez
            // Restrição 5: Apenas uma saída de cada cidade
            GRBLinExpr expr4 = 0;
            GRBLinExpr expr5 = 0;
            for (int i = 0; i < n; i++){
                if (i != s_initial){
                    expr4 += x[s_initial][i];
                    expr5 += x[i][s_initial];
                }
            }
            model.addConstr(expr4 == 1, "visit_" + to_string(s_initial));
            model.addConstr(expr5 == 1, "departure_" + to_string(s_initial));

            // Restrição 6: Cada cidade deve ser visitada no máximo uma vez
            // Restrição 7: Somatório de entradas na cidade j deve ser no máximo 1
            for (int j = 0; j < n; j++) {
                if (j != s_initial) {
                    GRBLinExpr expr6 = 0;
                    GRBLinExpr expr7 = 0;
                    for (int i = 0; i < n; i++) {
                        if (i != j) {
                            expr6 += x[i][j];
                            expr7 += x[j][i];
                        }
                    }
                    model.addConstr(expr6 <= 1, "visit_once_" + to_string(j));
                    model.addConstr(expr7 <= 1, "entry_once_" + to_string(j));
                }
            }

            // Restrição 8: Fluxo de entrada igual ao fluxo de saída
            for (int j = 0; j < n; j++){
                if (j != s_initial) {
                    GRBLinExpr expr81 = 0;
                    GRBLinExpr expr82 = 0;
                    for (int i = 0; i < n; i++) {
                        if (i != j) {
                            expr81 += x[i][j];
                            expr82 += x[j][i];
                        }
                    }
                    model.addConstr(expr81 - expr82 == 0, "flow_balance_" + to_string(j));
                }
            }

            // Restrição 9: Restrições MTZ para eliminação de subtours
            for (int i = 1; i < n; i++) {
                for (int j = 1; j < n; j++) {
                    if (i != j) {
                        model.addConstr(u[i] - u[j] + (n - 1) * x[i][j] <= n - 2, "mtz_" + to_string(i) + "_" + to_string(j));
                    }
                }
            }

            // Restrição 10: Somatório de x <= p para todas as cidades j exceto a cidade de origem
            for (int i = 0; i < n; i++){
                GRBLinExpr expr10 = 0;
                for (int j = 0; j < n; j++){
                    if (j != i){
                        expr10 += x[i][j];
                    }
                }
                model.addConstr(expr10 >= p[i], "max_visits_" + to_string(i));   
            }

            // Restrição 11: Duração máxima do tour para cada passageiro
            for (int k = 0; k < m; k++) {
                GRBQuadExpr expr11 = 0;

                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                            expr11 += h[k][i][j] * data.y[i][j];
                            expr11 += h[k][i][j] * p[i] * data.g[i];
                            
                    }
                }
                for (int i = 0; i < n; i++) {
                    if (i != data.passengers[k].dk)
                        expr11 += h[k][i][data.passengers[k].dk] * p[data.passengers[k].dk] * data.g[data.passengers[k].dk];
                }

                model.addQConstr(expr11 <= data.passengers[k].tk, "max_tour_duration_passenger_" + to_string(k));
            }

            // Restrição 12: Capacidade do veículo
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                    if (i != j){
                        GRBLinExpr expr12 = 0;
                        for (int k = 0; k < m; k++){
                            expr12 += h[k][i][j];
                        }
                        model.addConstr(expr12 <= data.vehicle_capacity * x[i][j], "capacity_" + to_string(i) + "_" + to_string(j));
                    }
                }
            }

            //RESTRICAO 19 - LINEARIZADA - CALCULO DA TAXA DE ALPHA
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        GRBLinExpr sum_hk = alpha[i][j];
                        for (int k = 0; k < m; k++) {
                            sum_hk += psi[i][j][k];
                        }
                        model.addConstr(sum_hk == 1.0, "alpha_" + to_string(i) + "_" + to_string(j));
                    }
                }
            }

            //RESTRICAO 13  LINEARIZADA: sum(i ∈ V) sum(j ∈ V - {i}) ((hkij * c) / (1 + sum(e ∈ L (heij))) <= wk para todo k ∈ L
            for (int k = 0; k < m; k++) {
                GRBLinExpr  expr13 = 0;
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        if (i != j) {
                            expr13 += data.c[i][j] * psi[i][j][k];
                            model.addConstr(psi[i][j][k] <= h[k][i][j], "droga1_"+ to_string(k)+ to_string(i)+to_string(j));
                            model.addConstr(psi[i][j][k] <= alpha[i][j], "droga2_" + to_string(k) + to_string(i) + to_string(j));
                            model.addConstr(psi[i][j][k] >= alpha[i][j] - 1 + h[k][i][j], "droga3_" + to_string(k) + to_string(i) + to_string(j));

                        }
                    }
                }
                model.addConstr(expr13 <= data.passengers[k].wk, "max_fare_" + to_string(k));

            }

            // Restrição 14: Cada passageiro deve ser deixado na cidade de destino e levado de sua origem:
            for (int k = 0; k < m; k++) {
                GRBLinExpr expr14 = 0;
                for (int i = 0; i < n; i++) {
                    if (i != data.passengers[k].ok){
                        expr14 += h[k][i][data.passengers[k].ok];
                    }
                    if (i != data.passengers[k].dk){
                        expr14 += h[k][data.passengers[k].dk][i];
                    }
                }

                model.addConstr(expr14 == 0, "pckup_dropoff_balance_" + to_string(k));
            }
        
            // Restrição 15: Se um passageiro é transportado de i para j, o veículo deve ir de i para j
            for (int k = 0; k < m; k++){
                if (data.passengers[k].ok != s_initial) {
                    GRBLinExpr expr15 = 0;
                    for (int i = 0; i < n; i++){
                        if (i != s_initial){
                            expr15 += h[k][s_initial][i];
                        }
                    }
                    model.addConstr(expr15 == 0, "hkis_" + to_string(k));
                }
            }

            // Restrição 16: Fluxo de entrada igual ao fluxo de saída para cada passageiro k em cada cidade i
            for (int k = 0; k < m; k++){
                for (int i = 0; i < n; i++){
                    if(i != data.passengers[k].ok && i != data.passengers[k].dk) {
                        GRBLinExpr expr161 = 0;
                        GRBLinExpr expr162 = 0;
                        for (int j = 0; j < n; j++){
                            if (i != j){
                                expr161 += h[k][j][i];  // fluxo de entrada
                                expr162 += h[k][i][j]; // fluxo de saída
                            }
                        }
                        model.addConstr(expr161 == expr162, "flow_balance_passenger_" + to_string(k) + "_" + to_string(i));
                    }
                }
            }


            model.optimize();
            
           // Verificar viabilidade e ignorar se o tempo foi excedido
            if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL && !cb.stop_optimization) {
                cout << "obj_cost: " << obj_cost.getValue() << " obj_time: " << obj_time.getValue() << " obj_bonus: " << obj_bonus.getValue() << endl;
                cout << endl << endl;
                Solution newSolution = {
                    obj_cost.getValue(),
                    obj_time.getValue(),
                    obj_bonus.getValue(),
                    extractTour(model, x, n),
                    collectBonuses(newSolution.tour, data, model, p),
                    collectPassengers(model, h, n, m)
                };

                solutions.push_back(newSolution);

                 // Verificar se passaram 20 minutos para salvar o progresso
                auto current_time = std::chrono::steady_clock::now();
                std::chrono::duration<double> elapsed_minutes = current_time - last_save_time;
                if (elapsed_minutes.count() >= 1200) {  // 20 minutos
                    saveProgress(solutions, "./", to_string(count_file) + "_" + fragment_file + "_file_20_min.txt");
                    last_save_time = current_time;  // Atualiza o último tempo de gravação
                    count_file++;
                }

            } else if (cb.stop_optimization) {
                cout << "Otimização interrompida devido ao limite de tempo. Solução ignorada." << endl;
                stop = cb.stop_optimization;
            } else if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
                cout << "Otimização finalizada!" << endl;
                stop = true;
            }

       }

    } catch (GRBException &e) {
        cerr << "Error code = " << e.getErrorCode() << endl;
        cerr << e.getMessage() << endl;
    } catch (...) {
        cerr << "Exception during optimization" << endl;
    }
}

// Função principal para o algoritmo
int main() {
    vector<string> typeInstancies = { "asymmetric", "symmetric" };

    vector<string> instances = {
        "10.1",
        "10.2",
        "10.3",
        "10.4",
        "10.5",
        "10.6",
        "20.1",
        "20.2",
        "20.3",
        "20.4",
        "20.5",
        "20.6",
    };

    for(string type : typeInstancies) {
        string outputFileType = "asym";
        if (type == "symmetric") {
            outputFileType = "sym";
        };

    for(string instance : instances)
    {
        auto start_time = chrono::high_resolution_clock::now();

        string path_input = "../instances/A3/"+type+"/"+instance+".in";
        string path_output = "../resultados/exato/a3/"+outputFileType+"/";
        string path_log = "../logs/exato/a3/"+outputFileType+"/";
        string file_name_output = instance+".txt";
        string file_log_name_output = instance+".log";
        string file_log_time_output = "time."+file_log_name_output;


        InputData data = readInput(path_input);

        GRBEnv env = GRBEnv(true);
        string log_file = path_log + file_log_name_output;
        env.set("LogFile", log_file);
        env.start();

        vector<Solution> solutions;

        solveModel(env, data, solutions, instance + "_" + type);

        ofstream outTimeLog(path_log + file_log_time_output);
        auto end_time = chrono::high_resolution_clock::now();
        chrono::duration<double> execution_time = end_time - start_time;

        outTimeLog << "Tempo de execução: " << execution_time.count() << " segundos\n";


        if (execution_time.count() >= 7200.0) {
            file_name_output = instance + ".2h.txt";
        }

        // Escrever resultados em solutions.txt
        ofstream outFile(path_output + file_name_output);
        if (outFile.is_open()) {
            outFile << "Fo(custo)\tFo(tempo)\tFo(bônus)\n";
            outFile << endl;

            for (const auto &solution : solutions) {
                outFile << fabs(solution.obj_cost) << "\t" << fabs(solution.obj_time) << "\t" << solution.obj_bonus << "\n";
                outFile << "tour: ";
                for (int city : solution.tour) {
                    outFile << city << " ";
                }
                outFile << "\n";
                outFile << "bonus: ";
                for (int bonus : solution.collectedBonuses) {
                    outFile << bonus << " ";
                }
                outFile << "\n";
                outFile << "passengers: ";
                for (int passenger : solution.passengers) {
                    outFile << passenger << " ";
                }
                outFile << "\n\n";
            }
            outFile.close();
        } else {
            cerr << "Unable to open file solutions.txt";
        }
        
        remove("mip1.log");
    }

    }
    return 0;
}