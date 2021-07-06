#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp17)]]

template <class Iterable, class func>
std::string concat(Iterable container, std::string delim, func f){
  std::string res;
  for (auto p = container.begin(); p != container.end(); ++p){
    if (p != container.begin())
      res += delim;
    res += f(*p);
  }
  return res;
}

std::string replace(const std::string& text, bool single = false){
  single = false;
  std::string res;
  std::string tmp;
  char cur_state = 0;
  bool replaceable = true;
  for (char c: text){
    if (c == '|')
      replaceable = false;
    if (c == '.' && replaceable){
      res += "Y[n]";
      continue;
    }
    if (cur_state != 0){
      if (c >= '0' && c <= '9') {
        tmp += c;
        continue;
      }
      else{
        if (cur_state == '@')
          res += "par_" + tmp;
        else if (single)
          res += "X[n]";
        else
          res += "X[n, " + tmp + "]";
        tmp = "";
      }
    }
    if (c == '@'){
      cur_state = '@';
    }
    else if (c == '$'){
      cur_state = '$';
    }
    else {
      cur_state = 0;
      res += c;
    }
  }
  if (cur_state == '@')
    res += "par_" + tmp;
  else if (cur_state == '$' && single)
    res += "X[n]";
  else if (cur_state == '$')
    res += "X[n, " + tmp + "]";
  return res;
}

struct Object{
  std::string Y_type;
  std::vector<std::string> X_name;
  std::vector<std::string> P_name;
  std::vector<std::string> P_type;
  std::vector<std::string> P_model;
  std::vector<int> S;
  std::map<std::pair<int, int>, std::vector<int>> _ZDS;
  std::map<int, std::string> S_model;
  std::map<std::pair<int, int>, std::string> Y_model;
  std::string to_string_data() const;
  std::string to_string_parameters() const;
  std::string _to_string_prior() const;
  std::string _to_string_S() const;
  std::string _to_string_len() const;
  std::string _to_string_latent() const;
  std::string to_string() const;
  Object ();
  void write_to_file(std::string filename) const;
};

std::string Object::to_string_data() const {
  std::vector<std::string> tmp;
  tmp.emplace_back("int<lower=1> N;");
  tmp.emplace_back("int<lower=0, upper=1> Z[N];");
  tmp.emplace_back("int<lower=0, upper=1> D[N];");
  if (X_name.size() >= 1)
    tmp.push_back(
      std::string("real X[N, ") +
        std::to_string(this->X_name.size()) +
        "];");
  else if (X_name.size() == 1)
    tmp.push_back(
      std::string("real X[N];"));
  std::string t;
  if (this->Y_type == "continuous")
    t = "real";
  else if (this->Y_type == "positive")
    t = "real<lower=0>";
  else if (this->Y_type == "binary")
    t = "int<lower=0, upper=1>";
  else if (this->Y_type == "count")
    t = "int<lower=0>";
  tmp.push_back(t + " Y[N];");
  std::string res("data {\n");
  for (auto & s : tmp){
    res += std::string(4, ' ') + s + "\n";
  }
  res += "}\n";
  return res;
}

Object::Object() {
  for (int z : {0, 1})
    for (int d: {0, 1})
      _ZDS[std::make_pair(z, d)] = std::vector<int>();
}

std::string Object::to_string_parameters() const {
  std::vector<std::string> tmp;
  for (int i = 0; i < this->P_name.size(); ++i){
    std::string t;
    if (this->P_type[i] == "continuous")
      t = "real";
    else if (this->P_type[i] == "positive")
      t = "real<lower=0>";
    tmp.push_back(t + " par_" + std::to_string(i + 1) + ';');
  }
  std::string res("parameters {\n");
  for (auto & s : tmp){
    res += std::string(4, ' ') + s + "\n";
  }
  res += "}\n";
  return res;
}

std::string Object::_to_string_prior() const {
  std::string res;
  for (int i = 0; i < this->P_name.size(); ++i){
    res += std::string(4, ' ') +
      "par_" + std::to_string(i + 1) + " ~ " +
      this->P_model[i] + ";\n";
  }
  return res;
}

std::string Object::_to_string_S() const {
  std::string res;
  for (auto& s: this->S_model){
    res += std::string(4, ' ') +
      "real log_prob_" + std::to_string(s.first) + " = " +
      s.second + ";\n";
  }
  res += std::string(4, ' ') +
    "real log_prob_all = log_sum_exp({" +
    concat(this->S_model, ", ",
           [](std::pair<int, std::string> p){
             return std::string("log_prob_") + std::to_string(p.first);
           }) +
             "});\n";
  return res;
}

std::string Object::_to_string_len() const {
  std::string blank(8, ' ');
  std::string res;
  res += blank + "int length;\n";
  bool _else = false;
  for (int z: {0, 1}){
    for (int d: {0, 1}){
      std::vector<int> tmp;
      if (this->_ZDS.find(std::make_pair(z, d)) != this->_ZDS.end())
        tmp = this->_ZDS.at(std::make_pair(z, d));
      if (!tmp.empty()){
        res += blank + (_else ? "else if ": "if ") + "(" +
          "Z[n] == " + std::to_string(z) + " && " +
          "D[n] == " + std::to_string(d) + ")\n";
        res += blank + std::string(4, ' ') +
          "length = " + std::to_string(tmp.size()) + ";\n";
      }
    }
  }
  return res;
}

std::string Object::_to_string_latent() const {
  std::string blank(12, ' ');
  std::string res;
  res += std::string(8, ' ') + "{\n";
  res += blank + "real log_l[length];\n";
  bool _else = false;
  for (int z: {0, 1}){
    for (int d: {0, 1}){
      std::vector<int> tmp;
      if (this->_ZDS.find(std::make_pair(z,d)) != this->_ZDS.end())
        tmp = this->_ZDS.at(std::make_pair(z, d));
      if (!tmp.empty()){
        res += blank + (_else ? "else if ": "if ") + "(" +
          "Z[n] == " + std::to_string(z) + " && " +
          "D[n] == " + std::to_string(d) + ") {\n";
        res += blank + std::string(4, ' ') + "// strata: " +
          concat(tmp, ", ", [](int c) {return std::to_string(c);}) +
          "\n";
        for (int i = 0; i < tmp.size(); ++i){
          res += blank + std::string(4, ' ') +
            "log_l[" + std::to_string(i + 1) + "] = " +
            "log_prob_" + std::to_string(tmp[i]) + " + " +
            this->Y_model.at(std::make_pair(tmp[i], z)) + ";\n";
        }
        res += blank + "}\n";
      }
    }
  }
  res += blank + "target += log_sum_exp(log_l) - log_prob_all;\n";
  res += std::string(8, ' ') + "}\n";
  return res;
}

std::string Object::to_string() const {
  std::string res;
  res += this->to_string_data();
  res += this->to_string_parameters();
  res += "model {\n";
  res += this->_to_string_prior();
  res += std::string(4, ' ') + "for (n in 1:N) {\n";
  res += this->_to_string_S();
  res += this->_to_string_len();
  res += this->_to_string_latent();
  res += std::string(4, ' ') + "}\n";
  res += "}\n";
  return res;
}

void Object::write_to_file(std::string filename) const {
  std::ofstream file(filename);
  file << this->to_string();
  file.close();
}

std::vector<std::string> readfile(std::string filename){
  std::vector<std::string> buffer;
  std::ifstream file(filename);
  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      buffer.push_back(line);
    }
    file.close();
  }
  return buffer;
}

std::vector<std::string> split(const std::string& str, char delim){
  std::vector<std::string> tokens;
  auto l = str.begin();
  auto r = str.begin();
  while (true) {
    while (*l == ' ' || *l == '\t' || *l == '\n')
      ++l, ++r;
    while (r != str.end() && *r != delim) {
      ++r;
    }
    tokens.emplace_back(l, r);
    if (*r == delim)
      l = ++r;
    else
      return tokens;
  }
}

template <typename T>
T parse_covariate(T it, Object& obj){
  while (*(++it) != "}"){
    obj.X_name.push_back(split(*it, ':')[1]);
  }
  return it;
}

template <typename T>
T parse_parameter(T it, Object& obj){
  while (*(++it) != "}"){
    auto tokens = split(split(*it, ':')[1], ' ');
    auto type = tokens[0].substr(1, tokens[0].size() - 2);
    obj.P_name.push_back(tokens[1]);
    obj.P_type.push_back(type);
  }
  return it;
}

template <typename T>
T parse_prior(T it, Object& obj) {
  while (*(++it) != "}") {
    obj.P_model.push_back(split(*it, '~')[1]);
  }
  return it;
}

template <typename T>
T parse_strata(T it, Object& obj) {
  while (*(++it) != "}") {
    auto tokens = split(*it, ':');
    obj.S_model[std::stoi(tokens[0])] = replace(tokens[1], obj.X_name.size() == 1);
  }
  return it;
}

template <typename T>
T parse_outcome(T it, Object& obj) {
  while (*(++it) != "}") {
    auto tokens = split(*it, ':');
    auto token_ids = split(tokens[0], ',');
    obj.Y_model[
    std::make_pair(std::stoi(token_ids[0]), std::stoi(token_ids[1]))
    ] = replace(tokens[1], obj.X_name.size() == 1);
  }
  return it;
}

Object parse(const std::vector<std::string>& buffer){
  Object obj;
  auto it = buffer.begin();
  while (it != buffer.end()) {
    if (it->rfind("Y:", 0) == 0) {
      // Y type
      obj.Y_type = split(*it, ':')[1];
    } else if (it->rfind("S:", 0) == 0) {
      // strata
      auto tokens = split(split(*it, ':')[1], ' ');
      for (auto token : tokens) {
        int num = std::stoi(token);
        obj.S.push_back(num);
        for (int z : {0, 1}){
          obj._ZDS[std::make_pair(z, z == 0? num >> 1: num & 1)].push_back(num);
        }
      }
    } else if (it->rfind("covariate", 0) == 0) {
      it = parse_covariate(it, obj);
    } else if (it->rfind("parameter", 0) == 0) {
      it = parse_parameter(it, obj);
    } else if (it->rfind("prior", 0) == 0){
      it = parse_prior(it, obj);
    } else if (it->rfind("strata", 0) == 0){
      it = parse_strata(it, obj);
    } else if (it->rfind("outcome", 0) == 0){
      it = parse_outcome(it, obj);
    }
    ++it;
  }
  return obj;
}

// [[Rcpp::export]]
int main() {
  std::string filename("model.pso");
  auto buffer = readfile(filename);
  auto obj = parse(buffer);
  obj.write_to_file("model.stan");
  return 0;
}