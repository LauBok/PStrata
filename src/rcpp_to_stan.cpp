#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp17)]]
std::vector<std::string> split (const std::string &s, std::string delimiter) {
  std::size_t pos_start = 0, pos_end, delim_len = delimiter.length();
  std::string token;
  std::vector<std::string> res;
  
  while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
    token = s.substr (pos_start, pos_end - pos_start);
    pos_start = pos_end + delim_len;
    res.push_back (token);
  }
  
  res.push_back (s.substr (pos_start));
  return res;
}

class Data {
private:
  std::map<std::pair<unsigned int, unsigned int>, 
           std::set<std::pair<unsigned int, unsigned int>>> SZDG_table;
  unsigned int S_max = 0, Z_max = 0, D_max = 0, G_max = 0;
  unsigned int S_re = 0, Y_re = 0;
  std::string Y_type, Y_family, Y_link;
  std::string func_family, func_link;
  std::vector<std::pair<std::string, std::string>> parameters;
  std::vector<std::tuple<std::string, std::string, std::vector<double>>> prior;
public:
  Data(const std::string& file_name) {
    std::fstream input(file_name);
    std::string ctrl;
    while (input >> ctrl) {
      if (ctrl == "SZDG") {
        unsigned int s, z, d, g;
        input >> s >> z >> d >> g;
        SZDG_table[std::make_pair(z, d)].insert(std::make_pair(s, g));
        S_max = std::max(S_max, s);
        Z_max = std::max(Z_max, z);
        D_max = std::max(D_max, d);
        G_max = std::max(G_max, g);
      }
      else if (ctrl == "Y") {
        input >> this->Y_family >> this->Y_link;
      }
      else if (ctrl == "random") {
        std::string var;
        input >> var;
        if (var == "S") input >> this->S_re;
        else if (var == "Y") input >> this->Y_re;
      }
      else if (ctrl == "prior") {
        std::string type, func;
        input >> type >> func;
        std::vector<double> tmp_param;
        int n;
        input >> n;
        for (int i = 0; i < n; ++i){
          double param;
          input >> param;
          tmp_param.push_back(param);
        }
        prior.emplace_back(type, func, tmp_param);
      }
    }
    input.close();
    input = std::fstream("auto_generated_files/family_info.txt");
    std::string type, family, family_func, link, link_func;
    int param_count;
    while (input >> family >> type >> family_func >> param_count) {
      for (int k = 0; k < param_count; ++k) {
        std::string param_name, param_type;
        input >> param_name >> param_type;
        if (family == Y_family)
          parameters.emplace_back(param_name, param_type);
      }
      if (family == Y_family) {
        func_family = family_func;
        Y_type = type;
        break;
      }
    }
    input.close();
    
    input = std::fstream("auto_generated_files/link_info.txt");
    while (input >> family >> link >> link_func) {
      if (family == Y_family && link == Y_link) {
        func_link = link_func;
        break;
      }
    }
    input.close();
  }
  
  std::string to_stan_functions() const {
    std::fstream input = std::fstream("auto_generated_files/function_implement.txt");
    std::string str;
    bool aim = false;
    for (std::string line; std::getline(input, line); ){
      if (line == "<<<>>>") aim = false;
      if (!aim) {
        auto splitted = split(line, " ");
        if (splitted[0] == "<<<") {
          std::string func_name = splitted[1];
          if (func_name == func_family || func_name == func_link) {
            aim = true;
          }
        }
        continue;
      }
      else {
        str += "    " + line + "\n";
      }
    }
    input.close();
    if (!str.empty())
      str = "functions {\n" + str + "}\n";
    return str;
  }
  
  std::string to_stan_data() const {
    std::string str;
    str += "data {\n";
    str += "    int<lower=0> N;\n";
    str += "    int<lower=0> PS;\n";
    str += "    int<lower=0> PG;\n";
    str += "    int<lower=0, upper=" + std::to_string(Z_max) + "> Z[N];\n";
    str += "    int<lower=0, upper=" + std::to_string(D_max) + "> D[N];\n";
    std::string Y_str;
    if (Y_type == "real") Y_str = "real Y[N];";
    else if (Y_type == "positive") Y_str = "real<lower=0> Y[N];";
    else if (Y_type == "binary") Y_str = "int<lower=0, upper=1> Y[N];";
    else if (Y_type == "count") Y_str = "int<lower=0> Y[N];";
    str += "    " + Y_str + "\n";
    str += "    matrix[N, PS] XS;\n";
    str += "    matrix[N, PG] XG;\n";
    for (int i = 0; i < S_re; ++i) {
      std::string tmp_P = "PS_RE_" + std::to_string(i + 1);
      std::string tmp_N = "NS_RE_" + std::to_string(i + 1);
      std::string tmp_X = "XS_RE_" + std::to_string(i + 1);
      str += "    int<lower=0> " + tmp_P + ";\n";
      str += "    int<lower=0> " + tmp_N + ";\n";
      str += "    matrix[N, " + tmp_P + "*" + tmp_N + "] " + tmp_X + ";\n";
    }
    for (int i = 0; i < Y_re; ++i) {
      std::string tmp_P = "PG_RE_" + std::to_string(i + 1);
      std::string tmp_N = "NG_RE_" + std::to_string(i + 1);
      std::string tmp_X = "XG_RE_" + std::to_string(i + 1);
      str += "    int<lower=0> " + tmp_P + ";\n";
      str += "    int<lower=0> " + tmp_N + ";\n";
      str += "    matrix[N, " + tmp_P + "*" + tmp_N + "] " + tmp_X + ";\n";
    }
    str += "}\n";
    return (str);
  };
  
  std::string to_stan_transformed_data() const {
    std::string str;
    str += "transformed data {\n";
    str += "    int S[" + std::to_string(G_max + 1) + "];\n";
    std::set<std::pair<unsigned int, unsigned int>> tmpset;
    for (auto zdsg: SZDG_table) {
      for (auto sg : zdsg.second) {
        tmpset.insert(sg);
      }
    }
    for (auto sg: tmpset){
      unsigned int s = sg.first, g = sg.second;
      str += "    S[" + std::to_string(g + 1) + "] = " + std::to_string(s + 1) + ";\n";
    }
    str += "}\n";
    return (str);
  }
  std::string to_stan_parameters() const {
    std::string S_str = std::to_string(S_max); // always zero for S == 0, so omit that.
    std::string G_str = std::to_string(1 + G_max);
    std::string str;
    str += "parameters {\n";
    str += "    matrix[" + S_str + ", PS] beta_S;\n";
    str += "    matrix[" + G_str + ", PG] beta_G;\n";
    for (int i = 0; i < S_re; ++i) {
      std::string tmp_P = "PS_RE_" + std::to_string(i + 1);
      std::string tmp_N = "NS_RE_" + std::to_string(i + 1);
      std::string tmp_X = "XS_RE_" + std::to_string(i + 1);
      std::string tmp_beta = "beta_S_RE_" + std::to_string(i + 1);
      str += "    matrix[" + S_str + ", " + tmp_P + "*" + tmp_N + "] " + tmp_beta + ";\n";
      str += "    real<lower = 0> tau_S_RE_" + std::to_string(i + 1) + "[" + S_str + ", " + tmp_P + "];\n";
    }
    for (int i = 0; i < Y_re; ++i) {
      std::string tmp_P = "PG_RE_" + std::to_string(i + 1);
      std::string tmp_N = "NG_RE_" + std::to_string(i + 1);
      std::string tmp_X = "XG_RE_" + std::to_string(i + 1);
      std::string tmp_beta = "beta_G_RE_" + std::to_string(i + 1);
      str += "    matrix[" + G_str + ", " + tmp_P + "*" + tmp_N + "] " + tmp_beta + ";\n";
      str += "    real<lower = 0> tau_G_RE_" + std::to_string(i + 1) + "[" + G_str + ", " + tmp_P + "];\n";
    }
    for (auto &p : parameters) {
      std::string type_str;
      if (p.second == "real") type_str = "real";
      else if (p.second == "positive") type_str = "real<lower=0>";
      str += "    " + type_str + " " + p.first + "[" + G_str + "];\n";
    }
    str += "}\n";
    return str;
  }
  
  std::string to_stan_transformed_parameters() const {
    std::string S_str = std::to_string(S_max); // always zero for S == 0, so omit that.
    std::string G_str = std::to_string(1 + G_max);
    std::string str;
    str += "transformed parameters {\n";
    for (int i = 0; i < S_re; ++i) {
      std::string tmp_P = "PS_RE_" + std::to_string(i + 1);
      std::string tmp_N = "NS_RE_" + std::to_string(i + 1);
      std::string tmp_beta = "beta_S_RE_" + std::to_string(i + 1);
      str += "    matrix[" + tmp_N + ", " + tmp_P + "] M_" + tmp_beta + "[" + S_str + "];\n";
    }
    for (int i = 0; i < Y_re; ++i) {
      std::string tmp_P = "PG_RE_" + std::to_string(i + 1);
      std::string tmp_N = "NG_RE_" + std::to_string(i + 1);
      std::string tmp_beta = "beta_G_RE_" + std::to_string(i + 1);
      str += "    matrix[" + tmp_N + ", " + tmp_P + "] M_" + tmp_beta + "[" + G_str + "];\n";
    }
    for (int i = 0; i < S_re; ++i) {
      str += "    for (i in 1:" + S_str + "){\n";
      for (int j = 0; j < S_max; ++j) {
        std::string tmp_P = "PS_RE_" + std::to_string(i + 1);
        std::string tmp_N = "NS_RE_" + std::to_string(i + 1);
        std::string tmp_beta = "beta_S_RE_" + std::to_string(i + 1);
        str += "        M_" + tmp_beta + "[" + std::to_string(j + 1) + "] = to_matrix(" + tmp_beta + "[" + std::to_string(j + 1) + "]" + ", " + tmp_N + ", " + tmp_P + ", 0);\n";
      }
      str += "    }\n";
    }
    for (int i = 0; i < Y_re; ++i) {
      str += "    for (i in 1:" + S_str + "){\n";
      for (int j = 0; j <= G_max; ++j) {
        std::string tmp_P = "PG_RE_" + std::to_string(i + 1);
        std::string tmp_N = "NG_RE_" + std::to_string(i + 1);
        std::string tmp_beta = "beta_G_RE_" + std::to_string(i + 1);
        str += "        M_" + tmp_beta + "[" + std::to_string(j + 1) + "] = to_matrix(" + tmp_beta + "[" + std::to_string(j + 1) + "]" + ", " + tmp_N + ", " + tmp_P + ", 0);\n";
      }
      str += "    }\n";
    }
    str += "}\n\n";
    return str;
  }
  
  std::string to_stan_model() const {
    std::string str_sp1 = std::to_string(S_max + 1);
    std::string str;
    str += "model {\n";
    
    // random effect
    
    for (int i = 0; i < S_re; ++i) {
      std::string tmp_P = "PS_RE_" + std::to_string(i + 1);
      std::string tmp_N = "NS_RE_" + std::to_string(i + 1);
      std::string tmp_X = "XS_RE_" + std::to_string(i + 1);
      std::string tmp_tau = "tau_S_RE_" + std::to_string(i + 1);
      std::string tmp_beta = "beta_S_RE_" + std::to_string(i + 1);
      str += "    for (i in 1:" + std::to_string(S_max) + ") {\n";
      str += "        " + tmp_tau + "[i] ~ inv_gamma(1, 1);\n";
      str += "        for (j in 1:" + tmp_P + ") {\n";
      str += "            M_" + tmp_beta + "[i][:, j] ~ normal(0, " + tmp_tau + "[i, j]); \n";
      str += "        }\n";
      str += "    }\n";
    }
    
    for (int i = 0; i < Y_re; ++i) {
      std::string tmp_P = "PG_RE_" + std::to_string(i + 1);
      std::string tmp_N = "NG_RE_" + std::to_string(i + 1);
      std::string tmp_X = "XG_RE_" + std::to_string(i + 1);
      std::string tmp_tau = "tau_G_RE_" + std::to_string(i + 1);
      std::string tmp_beta = "beta_G_RE_" + std::to_string(i + 1);
      str += "    for (i in 1:" + std::to_string(G_max + 1) + ") {\n";
      str += "        " + tmp_tau + "[i] ~ inv_gamma(1, 1);\n";
      str += "        for (j in 1:" + tmp_P + ") {\n";
      str += "            M_" + tmp_beta + "[i][:, j] ~ normal(0, " + tmp_tau + "[i, j]); \n";
      str += "        }\n";
      str += "    }\n";
    }
    
    // prior
    for (auto& prior_info : prior) {
      std::string type = std::get<0>(prior_info);
      std::string func = std::get<1>(prior_info);
      std::vector<double> params = std::get<2>(prior_info);
      std::string params_str;
      for (auto param: params) {
        if (!params_str.empty())
          params_str += ", ";
        params_str += std::to_string(param);
      }
      if (func == "flat")
        continue;
      if (type == "intercept") {
        str += "    beta_S[:, 1] ~ " + func + "(" + params_str + ");\n";
        str += "    beta_G[:, 1] ~ " + func + "(" + params_str + ");\n";
      }
      else if (type == "coefficient") {
        str += "    to_vector(beta_S[:, 2:PS]) ~ " + func + "(" + params_str + ");\n";
        str += "    to_vector(beta_G[:, 2:PG]) ~ " + func + "(" + params_str + ");\n";
      }
      else {
        bool found = false;
        for (auto& param : parameters)
          if (type == param.first) {
            found = true;
            break;
          }
          if (found) {
            str += "    " + type + " ~ " + func + "(" + params_str + ");\n";
          }
      }
    }
    
    str += "    for (n in 1:N) {\n";
    str += "        int length;\n";
    str += "        real log_prob[" + str_sp1 + "];\n";
    str += "        log_prob[1] = 0;\n";
    str += "        for (s in 2:" + str_sp1 + ") {\n";
    str += "            log_prob[s] = XS[n] * beta_S[s - 1]'";
    for (int i = 0; i < S_re; ++i) {
      str += " + XS_RE_" + std::to_string(i + 1) + "[n] * beta_S_RE_" + std::to_string(i + 1) + "[s - 1]'";
    }
    str += ";\n";
    str += "        }\n";
    bool b_else = false;
    for (auto &p : SZDG_table) {
      auto z = p.first.first, d = p.first.second;
      auto sg_set = p.second;
      str += "        " + (b_else?std::string("else "):std::string("")) + "if (" +
        "Z[n] == " + std::to_string(z) + " && D[n] == " + std::to_string(d) + ")\n";
      str += "            length = " + std::to_string(sg_set.size()) + ";\n";
      b_else = true;
    }
    str += "        {\n";
    str += "            real log_l[length];\n";
    b_else = false;
    for (auto &p : SZDG_table) {
      auto z = p.first.first, d = p.first.second;
      auto sg_set = p.second;
      str += "            " + (b_else?std::string("else "):std::string("")) + "if (" +
        "Z[n] == " + std::to_string(z) + " && D[n] == " + std::to_string(d) + ") {\n";
      str += "                // strata:";
      for (auto &sg : sg_set) str += " " + std::to_string(sg.first);
      str += "\n";
      int i = 0;
      for (auto &sg : sg_set) {
        unsigned int s = sg.first, g = sg.second;
        std::string str_param;
        for (auto &p : parameters) {
          str_param += ", " + p.first + "[" + std::to_string(g + 1) + "]";
        }
        str += "                log_l[" + std::to_string(++i) + "] = log_prob[" +
          std::to_string(s + 1) + "] + " +
          func_family + "(Y[n] | " + func_link + "(" +
          "XG[n] * beta_G[" + std::to_string(g + 1) + "]'";
        for (int i = 0; i < Y_re; ++i) {
          str += " + XG_RE_" + std::to_string(i + 1) + "[n] * beta_G_RE_" + std::to_string(i + 1) + "[" + std::to_string(g + 1) + "]'";
        }
        str += ")" + str_param + ");\n";
      }
      str += "            }\n";
      b_else = true;
    }
    str += "            target += log_sum_exp(log_l) - log_sum_exp(log_prob);\n";
    str += "        }\n";
    str += "    }\n";
    str += "}\n";
    return str;
  }
  std::string to_stan_generated_quantities() const {
    std::string S_str = std::to_string(S_max + 1);
    std::string G_str = std::to_string(G_max + 1);
    std::string str = "generated quantities {\n";
    str += "    vector[" + S_str + "] strata_prob;\n";
    str += "    vector[" + G_str + "] mean_effect;\n";
    str += "    {\n";
    str += "        matrix[N, " + S_str + "] log_prob;\n";
    // str += "        vector[" + S_str + "] denom;\n";
    str += "        vector[" + G_str + "] numer;\n";
    str += "        matrix[N, " + G_str + "] expected_mean;\n";
    str += "        for (i in 1:N)\n";
    str += "            for (j in 1:" + G_str + ")\n";
    str += "                expected_mean[i, j] = " + func_link + 
      "(XG[i] * beta_G[j]'";
    for (int i = 0; i < Y_re; ++i) {
      str += " + XG_RE_" + std::to_string(i + 1) + "[i] * beta_G_RE_" + std::to_string(i + 1) + "[j]'";
    }
    str += ");\n";
    str += "        log_prob[:, 1] = rep_vector(0, N);\n";
    str += "        log_prob[:, 2:" + S_str + "] = XS * beta_S'";
    for (int i = 0; i < S_re; ++i) {
      str += " + XS_RE_" + std::to_string(i + 1) + " * beta_S_RE_" + std::to_string(i + 1) + "'";
    }
    str += ";\n";
    str += "        for (n in 1:N) {\n";
    str += "            log_prob[n] -= log_sum_exp(log_prob[n]);\n";
    str += "        }\n";
    str += "        for (s in 1:" + S_str + ") strata_prob[s] = mean(exp(log_prob[:, s]));\n";
    str += "        for (g in 1:" + G_str + ") {\n";
    str += "            numer[g] = mean(expected_mean[:, g] .* exp(log_prob[:, S[g]]));\n";
    str += "            mean_effect[g] = numer[g] / strata_prob[S[g]];\n";
    str += "        }\n";
    str += "    }\n";
    str += "}\n";
    return (str);
  }
  void print_info() const {
    std::cout << S_max << ' ' << Z_max << ' ' << D_max << ' ' << G_max << std::endl;
    std::cout << Y_type << ' ' << Y_family << ' ' << Y_link << std::endl;
  }
};

// [[Rcpp::export]]
int to_stan(const std::string& name) {
  Data data(name + ".pso");
  std::fstream output(name + ".stan", std::fstream::out | std::fstream::trunc);
  output << data.to_stan_functions() << std::endl;
  output << data.to_stan_data() << std::endl;
  output << data.to_stan_transformed_data() << std::endl;
  output << data.to_stan_parameters() << std::endl;
  output << data.to_stan_transformed_parameters() << std::endl;
  output << data.to_stan_model() << std::endl;
  output << data.to_stan_generated_quantities() << std::endl;
  output.close();
  return 0;
}