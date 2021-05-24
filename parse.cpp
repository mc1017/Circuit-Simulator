#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <complex>

class ImpedanceDevice{
public:
  virtual std::complex<double> get_impedance(double omega) const = 0;

  virtual std::string show_nodeinfo() const = 0;

  virtual ~ImpedanceDevice() { }
};

class Resistor : public ImpedanceDevice{

public:

  Resistor(int n1, int n2, double r) : node1(n1), node2(n2), resistance(r) { }

  std::string show_nodeinfo() const {
    return "Nodal Coordinates: (" + std::to_string(node1) + ", " + std::to_string(node2) + ")";
  }

  std::complex<double> get_impedance(double omega) const {
    std::complex<double> impedance(resistance);

    return impedance;
  }

private:
  int node1;
  int node2;
  double resistance;
};

class Capacitor : public ImpedanceDevice{
public:

  Capacitor(int n1, int n2, double c) : node1(n1), node2(n2), capacitance(c) { }

  std::string show_nodeinfo() const {
    return "Nodal Coordinates: (" + std::to_string(node1) + ", " + std::to_string(node2) + ")";
  }
  std::complex<double>get_impedance(double omega) const {
    std::complex<double> impedance(0, - 1/(omega * capacitance));

    return impedance;
  }
 
private :
  int node1;
  int node2;
  double capacitance;
 
};

class Inductor : public ImpedanceDevice{
public:
  Inductor(int n1, int n2, double l) : node1(n1), node2(n2), inductance(l) { }

  std::string show_nodeinfo() const {
    return "Nodal Coordinates: (" + std::to_string(node1) + ", " + std::to_string(node2) + ")";
  }

  std::complex<double>get_impedance(double omega) const {
    std::complex<double> impedance(0, (omega * inductance));

    return impedance;
  }

private:
  int node1;
  int node2;
  double inductance;
};

class Source{
public:
  virtual double get_value() const = 0;

  virtual std::string show_nodeinfo() const = 0;

  virtual std::string get_type() const = 0;

  virtual ~Source() { }
};

class DCVSource : public Source{
public:
  DCVSource(int n_p, int n_m, double v) : node_plus(n_p), node_minus(n_m), voltage(v) { }

  std::string show_nodeinfo() const {
    return "Nodal Coordinates: (" + std::to_string(node_plus) + ", " + std::to_string(node_minus) + ")";
  }

  double get_value() const {
    return voltage;
  }

  std::string get_type() const {
    return "DC V";
  }

private:
  int node_plus;
  int node_minus;
  double voltage;
};

class DCISource : public Source{
public:
  DCISource(int n_in, int n_out, double i) : node_in(n_in), node_out(n_out), current(i){ }

  std::string show_nodeinfo() const {
    return "Nodal Coordinates: (" + std::to_string(node_in) + ", " + std::to_string(node_out) + ")";
  }

  double get_value() const {
    return current;
  }

  std::string get_type() const {
    return "DC I";
  }

private:
  int node_in;
  int node_out;
  double current;
};

int node_to_number(std::string node){
  std::string node_label;
  int node_number;

  if(node != "0"){

    for(int i = 1; i < node.size(); i++){
      node_label.push_back(node[i]);
    }

    node_number = std::stoi(node_label);

    return node_number;
  }

  return 0;
}

double prefix_convertor(std::string value){
  std::string quantity;
  std::string prefix;
  double num_quantity;

  for(int i = 0; i < value.size(); i++){
    if(std::isdigit(value[i])|| value[i] == '.'){
      quantity.push_back(value[i]);
    }
  }

  for(int i = 0; i < value.size(); i++){
    if(std::isdigit(value[i]) == false){
      prefix.push_back(value[i]);
    }
  }

  num_quantity = std::stod(quantity);

  if(prefix == "p"){
    num_quantity = num_quantity * pow(10, -12);
  }
  else if(prefix == "n"){
    num_quantity = num_quantity * pow(10, -9);
  }
  else if(prefix == "\u00B5"){
    num_quantity = num_quantity * pow(10, -6);
  }
  else if(prefix == "m"){
    num_quantity = num_quantity * pow(10, -3);
  }
  else if(prefix == "k"){
    num_quantity = num_quantity * pow(10, 3);
  }
  else if(prefix == "Meg"){
    num_quantity = num_quantity * pow(10, 6);
  }
  else if(prefix == "G"){
    num_quantity = num_quantity * pow(10, 9);
  }
  else{
    return num_quantity;
  }

  return num_quantity;
}

int main(){
  std::ifstream infile; 
  infile.open("testlist.txt");
 
  if(!infile.is_open()){
    std::cout << "error opening file" << std::endl;
    return EXIT_FAILURE;
  }
 
  std::string component;
  std::vector<std::string> substrs;
  std::vector<ImpedanceDevice*> impedance_devices; 
  ImpedanceDevice* tmp_id;
  std::vector<Source*> sources; 
  Source* tmp_s;
 
  while(std::getline(infile, component)){
    std::stringstream line(component);

    while(line.good()){
      std::string substr;
      std::getline(line, substr, ' ');

      substrs.push_back(substr);
    }

    if(substrs[0][0] == 'R'){
      tmp_id = new Resistor(node_to_number(substrs[1]), node_to_number(substrs[2]), prefix_convertor(substrs[3]));
            
      impedance_devices.push_back(tmp_id);
    }
    else if(substrs[0][0] == 'C'){
      tmp_id = new Capacitor(node_to_number(substrs[1]), node_to_number(substrs[2]), prefix_convertor(substrs[3]));

      impedance_devices.push_back(tmp_id);
    }
    else if(substrs[0][0] == 'L'){
      tmp_id = new Inductor(node_to_number(substrs[1]), node_to_number(substrs[2]), prefix_convertor(substrs[3]));

      impedance_devices.push_back(tmp_id);
    }
    else if(substrs[0][0] == 'V'){
      tmp_s = new DCVSource(node_to_number(substrs[1]), node_to_number(substrs[2]), prefix_convertor(substrs[3]));

      sources.push_back(tmp_s);
    }
    else if(substrs[0][0] == 'I'){
      tmp_s = new DCISource(node_to_number(substrs[1]), node_to_number(substrs[2]), prefix_convertor(substrs[3]));

      sources.push_back(tmp_s);
    }
        

    substrs.clear();
  }

  infile.close();
    
  for(int i = 0; i < impedance_devices.size(); i++){
    std::cout << impedance_devices[i]->show_nodeinfo() << std::endl;
    std::cout << "Impedance: " << impedance_devices[i]->get_impedance(1) << std::endl;
  }

  for(int i = 0; i < sources.size(); i++){
    std::cout << sources[i]->show_nodeinfo() << std::endl;
    std::cout << "Source Value: " << sources[i]->get_value() << std::endl;
    std::cout << "Source Type: " << sources[i]->get_type() << std::endl;
  }

}

