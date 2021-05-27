double extract_double(std::string label){
    std::string double_string;

    for(int i = 0; i < label.size(); i++){
        if(std::isdigit(label[i]) || label[i] == '.'){
         double_string.push_back(label[i]);
        }
    }

    return std::stod(double_string);
}

double prefix_convertor(std::string value){
    std::string prefix;
    double num_quantity;

    for(int i = 0; i < value.size(); i++){
        if(std::isdigit(value[i]) == false && value[i] != '.'){
        prefix.push_back(value[i]);
        }
    }

    num_quantity = extract_double(value);

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