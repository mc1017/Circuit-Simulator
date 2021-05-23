#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>


struct Resistor{
  int node1;
  int node2;
  double resistance;
};

int node_to_number(std::string node){
  std::string number;
  int hehe;

  if(node != "0"){

    for(int i = 1; i < node.size(); i++){
      number.push_back(node[i]);
    }

      hehe = std::stoi(number);

    return hehe;
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
    std::vector<std::string> haha;
    std::vector<Resistor> Two_T; 
 
    while(std::getline(infile, component)){
        std::stringstream hehe(component);

        while(hehe.good()){
          std::string substr;

          std::getline(hehe, substr, ' ');

          haha.push_back(substr);
        }

        if(haha[0][0] == 'R'){
            Resistor R;
            R.node1 = node_to_number(haha[1]);
            R.node2 = node_to_number(haha[2]);
            R.resistance = prefix_convertor(haha[3]);
            
            Two_T.push_back(R);
        }
        else if(haha[0][0] == 'C'){
            Resistor R;
            R.node1 = node_to_number(haha[1]);
            R.node2 = node_to_number(haha[2]);
            R.resistance = prefix_convertor(haha[3]);

            Two_T.push_back(R);
        }

        haha.clear();
    }

    std::cout << Two_T[0].node1 << std::endl;
    std::cout << Two_T[0].node2 << std::endl;
    std::cout << Two_T[0].resistance << std::endl;
    std::cout << Two_T[1].node1 << std::endl;
    std::cout << Two_T[1].node2 << std::endl;
    std::cout << Two_T[1].resistance << std::endl;

    infile.close();
}

