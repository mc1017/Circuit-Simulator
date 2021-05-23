#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>


class Resistor{

public:

    Resistor(int n1, int n2, double val) : node1(n1), node2(n2), resistance(val){

    }

    std::string show_nodeinfo() const {
        return "Nodal Coordinates: (" + std::to_string(node1) + ", " + std::to_string(node2) + ")";
    }

    double get_resistance() const {
        return resistance;
    }

private:
    int node1;
    int node2;
    double resistance;
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
    std::vector<Resistor> Two_T; 
 
    while(std::getline(infile, component)){
        std::stringstream line(component);

        while(line.good()){
          std::string substr;

          std::getline(line, substr, ' ');

          substrs.push_back(substr);
        }

        if(substrs[0][0] == 'R'){
            Resistor R(node_to_number(substrs[1]), node_to_number(substrs[2]), prefix_convertor(substrs[3]));
            
            Two_T.push_back(R);
        }
        else if(substrs[0][0] == 'C'){
            Resistor R(node_to_number(substrs[1]), node_to_number(substrs[2]), prefix_convertor(substrs[3]));

            Two_T.push_back(R);
        }

        substrs.clear();
    }

    std::cout << Two_T[0].show_nodeinfo() << std::endl;
    std::cout << "Reistance: " << Two_T[0].get_resistance() << std::endl;
    std::cout << Two_T[1].show_nodeinfo() << std::endl;
    std::cout << "Reistance: " << Two_T[1].get_resistance() << std::endl;

    infile.close();
}

