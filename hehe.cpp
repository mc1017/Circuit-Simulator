#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <complex>
#include <float.h>
#include "library/Eigen/Dense"

struct NodePoint{
    int x;
    int y;
};

struct NodeTri{
    int x;
    int y;
    int z;
};

class ImpedanceDevice{
public:
    virtual std::string show_nodeinfo() const = 0;

    virtual std::complex<double> get_impedance(double omega) const = 0;

    virtual std::complex<double> get_conductance(double omega) const = 0;

    virtual NodePoint give_nodeinfo() const = 0;

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

    std::complex<double> get_conductance(double omega) const {
        std::complex<double> conductance(1/resistance);

        return conductance;
    }

    NodePoint give_nodeinfo() const {
        NodePoint N;

        N.x = node1;
        N.y = node2;

        return N;
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

    std::complex<double> get_conductance(double omega) const {
        std::complex<double> conductance(0, omega * capacitance);

        return conductance;
    }

    NodePoint give_nodeinfo() const {
        NodePoint N;

        N.x = node1;
        N.y = node2;

        return N;
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

    std::complex<double> get_conductance(double omega) const {
        std::complex<double> conductance(0, -1/(omega * inductance));

        return conductance;
    }

    NodePoint give_nodeinfo() const {
        NodePoint N;

        N.x = node1;
        N.y = node2;

        return N;
    }

private:
    int node1;
    int node2;
    double inductance;
};

class Source{
public:
    virtual std::string show_nodeinfo() const = 0;

    virtual double get_magnitude() const = 0;

    virtual double get_phase() const = 0;

    virtual double get_gm() const = 0;

    virtual std::string get_type() const = 0;

    virtual NodePoint give_nodeinfo() const = 0;

    virtual NodePoint give_controlinfo() const = 0;

    virtual ~Source() { }
};

class DCVSource : public Source{
public:

    DCVSource(int n_p, int n_m, double v) : node_plus(n_p), node_minus(n_m), voltage(v) { }

    std::string show_nodeinfo() const {
        return "Nodal Coordinates: (" + std::to_string(node_plus) + ", " + std::to_string(node_minus) + ")";
    }

    double get_magnitude() const {
        return voltage;
    }

    double get_phase() const { 
        return 0;
    }

    double get_gm() const {
        return 0;
    }

    std::string get_type() const {
        return "DC V";
    }

    NodePoint give_nodeinfo() const {
        NodePoint N;

        N.x = node_plus;
        N.y = node_minus;

        return N;
    }

    NodePoint give_controlinfo() const {
        NodePoint C;

        C.x = 0;
        C.y = 0;

        return C;
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

    double get_magnitude() const {
        return current;
    }

    double get_phase() const { 
        return 0;
    }

    double get_gm() const {
        return 0;
    }

    std::string get_type() const {
        return "DC I";
    }

    NodePoint give_nodeinfo() const {
        NodePoint N;

        N.x = node_in;
        N.y = node_out;

        return N;
    }

    NodePoint give_controlinfo() const {
        NodePoint C;

        C.x = 0;
        C.y = 0;

        return C;
    }

private:
    int node_in;
    int node_out;
    double current;
};

class ACVSource : public Source{
public:

    ACVSource(int n_p, int n_m, double v_m, double v_p) : node_plus(n_p), node_minus(n_m), amplitude(v_m), phase(v_p){ }

    std::string show_nodeinfo() const {
        return "Nodal Coordinates: (" + std::to_string(node_plus) + ", " + std::to_string(node_minus) + ")";
    }

    double get_magnitude() const {
        return amplitude;
    }

    double get_phase() const {
        return phase;
    }

    double get_gm() const {
        return 0;
    }

    std::string get_type() const {
        return "AC V";
    }

    NodePoint give_nodeinfo() const {
        NodePoint N;

        N.x = node_plus;
        N.y = node_minus;

        return N;
    }

    NodePoint give_controlinfo() const {
        NodePoint C;

        C.x = 0;
        C.y = 0;

        return C;
    }

private:
    int node_plus;
    int node_minus;
    double amplitude;
    double phase;
};

class ACISource : public Source{
public:

    ACISource(int n_in, int n_out, double i_m, double i_p) : node_in(n_in), node_out(n_out), amplitude(i_m), phase(i_p){ }

    std::string show_nodeinfo() const {
        return "Nodal Coordinates: (" + std::to_string(node_in) + ", " + std::to_string(node_out) + ")";
    }

    double get_magnitude() const {
        return amplitude;
    }

    double get_phase() const {
        return phase;
    }

    double get_gm() const {
        return 0;
    }

    std::string get_type() const {
        return "AC I";
    }

    NodePoint give_nodeinfo() const {
        NodePoint N;

        N.x = node_in;
        N.y = node_out;

        return N;
    }

    NodePoint give_controlinfo() const {
        NodePoint C;

        C.x = 0;
        C.y = 0;

        return C;
    }

private:
    int node_in;
    int node_out;
    double amplitude;
    double phase;
};

class VCCSource : public Source{
public:

    VCCSource(int n_p, int n_m, int c_p, int c_m, double g_m) : node_plus(n_p), node_minus(n_m), control_plus(c_p), control_minus(c_m), transconductance(g_m){ }

    std::string show_nodeinfo() const {
        return "Nodal Coordinates: (" + std::to_string(node_plus) + ", " + std::to_string(node_minus) + ")";
    }

    double get_magnitude() const {
        return 0;
    }

    double get_phase() const {
        return 0;
    }

    double get_gm() const {
        return transconductance;
    }

    std::string get_type() const {
        return "VCCS";
    }

    NodePoint give_nodeinfo() const {
        NodePoint N;

        N.x = node_plus;
        N.y = node_minus;

        return N;
    }

    NodePoint give_controlinfo() const {
        NodePoint C;

        C.x = control_plus;
        C.y = control_minus;

        return C;
    }

private:
    int node_plus;
    int node_minus;
    int control_plus;
    int control_minus;
    double transconductance;
};

class NonLinearDevice{
public:
    virtual std::string show_nodeinfo() const = 0;

    virtual NodePoint give_binodeinfo() const = 0;

    virtual NodeTri give_trinodeinfo() const = 0;

    virtual std::string get_model() const = 0;

    virtual ~NonLinearDevice() { }
};

class Diode: public NonLinearDevice{
public:

    Diode(int n_a, int n_c, std::string mod) : node_an(n_a), node_cat(n_c), model(mod){ }
    
    std::string show_nodeinfo() const {
        return "Nodal Coordinates: (" + std::to_string(node_an) + ", " + std::to_string(node_cat) +")";
    }

    NodePoint give_binodeinfo() const {
        NodePoint N;

        N.x = node_an;
        N.y = node_cat;

        return N;
    }

    NodeTri give_trinodeinfo() const {
        NodeTri N;

        N.x = 0;
        N.y = 0;
        N.z = 0;

        return N;
    }

    std::string get_model() const {
        return model;
    }

private:
    int node_an;
    int node_cat;
    std::string model;
};

class BJT: public NonLinearDevice{
public: 

    BJT(int n_c, int n_b, int n_e, std::string mod) : node_c(n_c), node_b(n_b), node_e(n_e), model(mod) { }

    std::string show_nodeinfo() const {
        return "Nodal Coordinates: (" + std::to_string(node_c) + ", " + std::to_string(node_b) + ", " + std::to_string(node_e) +")";
    }

    NodePoint give_binodeinfo() const {
        NodePoint N;

        N.x = 0;
        N.y = 0;

        return N;
    }

    NodeTri give_trinodeinfo() const {
        NodeTri N;

        N.x = node_c;
        N.y = node_b;
        N.z = node_e;

        return N;
    }

    std::string get_model() const {
        return model;
    }

private:
    int node_c;
    int node_b;
    int node_e;
    std::string model;
};

class MOSFET: public NonLinearDevice{
public: 

    MOSFET(int n_d, int n_g, int n_s, std::string mod) : node_d(n_d), node_g(n_g), node_s(n_s), model(mod) { }

    std::string show_nodeinfo() const {
        return "Nodal Coordinates: (" + std::to_string(node_d) + ", " + std::to_string(node_g) + ", "  + std::to_string(node_s) + ")";
    }

    NodePoint give_binodeinfo() const {
        NodePoint N;

        N.x = 0;
        N.y = 0;

        return N;
    }

    NodeTri give_trinodeinfo() const {
        NodeTri N;

        N.x = node_d;
        N.y = node_g;
        N.z = node_s;

        return N;
    }

    std::string get_model() const {
        return model;
    }

private:
    int node_d;
    int node_g;
    int node_s;
    std::string model;
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

double get_AC_magnitude(std::string ac_param){
    std::string magnitude;

    magnitude = std::to_string(extract_double(ac_param)) + ac_param[ac_param.size()-1];

    return prefix_convertor(magnitude);
}

int max_node_number(int node_max, int node1, int node2, int node3, int node4){
    
    if((node1 > node_max) && (node1 >= node2) && (node1 >= node3) && (node1 >= node4)){
        return node1;
    }
    else if((node2 > node_max) && (node2 >= node1) && (node2 >= node3) && (node2 >= node4)){
        return node2;
    }
    else if((node3 > node_max) && (node3 >= node1) && (node3 >= node2) && (node3 >= node4)){
        return node3;
    }
    else if((node4 > node_max) && (node4 >= node1) && (node4 >= node2) && (node4 >= node3)){
        return node4;
    }

    return node_max;
}

double return_tf_magnitude(std::complex<double> source, std::complex<double> output_node){
    double source_mag, output_mag, gain;

    source_mag = std::sqrt(pow(std::real(source), 2) + pow(std::imag(source),2));
    output_mag = std::sqrt(pow(std::real(output_node), 2) + pow(std::imag(output_node),2));
    gain = output_mag / source_mag;

    return gain;
}

double return_tf_phase(std::complex<double> source, std::complex<double> output_node){
    double phase_change;

    phase_change = (std::arg(output_node) - std::arg(source)) * 180 / M_PI;

    return phase_change;
}

using namespace Eigen;

int main(){
    std::ifstream infile; 
    infile.open("lpf.txt");
 
    if(!infile.is_open()){
        std::cout << "error opening file" << std::endl;
        return EXIT_FAILURE;
    }
 
    std::string component;
    std::vector<std::string> substrs;

    std::vector<ImpedanceDevice*> impedance_devices, ss_impedance_devices; 
    ImpedanceDevice* tmp_id;
    ImpedanceDevice* tmp_id2;

    std::vector<Source*> sources, ss_sources; 
    Source* tmp_s;

    std::vector<NonLinearDevice*> non_linear_devices;
    NonLinearDevice* tmp_nld;
    // frequency step parameters, we assume ac analysis always done in decades
    double f_start, f_stop, n_ppd, f, omega;
    int n_max = 0;
 
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
            ss_impedance_devices.push_back(tmp_id);
            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), 0, 0);
        }
        else if(substrs[0][0] == 'C'){
            tmp_id = new Capacitor(node_to_number(substrs[1]), node_to_number(substrs[2]), prefix_convertor(substrs[3]));

            impedance_devices.push_back(tmp_id);
            ss_impedance_devices.push_back(tmp_id);
            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), 0, 0);
        }
        else if(substrs[0][0] == 'L'){
            tmp_id = new Inductor(node_to_number(substrs[1]), node_to_number(substrs[2]), prefix_convertor(substrs[3]));

            impedance_devices.push_back(tmp_id);
            ss_impedance_devices.push_back(tmp_id);
            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), 0, 0);
        }
        else if(substrs[0][0] == 'V' && substrs[3][0] != 'A'){
            tmp_s = new DCVSource(node_to_number(substrs[1]), node_to_number(substrs[2]), prefix_convertor(substrs[3]));
            tmp_id = new Resistor(node_to_number(substrs[1]), node_to_number(substrs[2]), DBL_MIN);
            //ss equivalent of DC Voltage Source is short circuit, represented by resistor of least possible value in C++
            sources.push_back(tmp_s);
            ss_impedance_devices.push_back(tmp_id);
            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), 0, 0);
        }
        else if(substrs[0][0] == 'V' && substrs[3][0] == 'A'){
            tmp_s = new ACVSource(node_to_number(substrs[1]), node_to_number(substrs[2]), get_AC_magnitude(substrs[3]), extract_double(substrs[4]));

            sources.push_back(tmp_s);
            ss_sources.push_back(tmp_s);
            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), 0, 0);
        }
        else if(substrs[0][0] == 'I' && substrs[3][0] != 'A'){
            tmp_s = new DCISource(node_to_number(substrs[1]), node_to_number(substrs[2]), prefix_convertor(substrs[3]));
            
            sources.push_back(tmp_s);
            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), 0, 0);
        }
        else if(substrs[0][0] == 'I' && substrs[3][0] == 'A'){
            tmp_s = new ACISource(node_to_number(substrs[1]), node_to_number(substrs[2]), get_AC_magnitude(substrs[3]), extract_double(substrs[4]));

            sources.push_back(tmp_s);
            ss_sources.push_back(tmp_s);
            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), 0, 0);
        }
        else if(substrs[0][0] == 'G'){
            tmp_s = new VCCSource(node_to_number(substrs[1]), node_to_number(substrs[2]), node_to_number(substrs[3]), node_to_number(substrs[4]), prefix_convertor(substrs[5]));

            sources.push_back(tmp_s);
            //what is ss equavalent of VCCS?
            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), node_to_number(substrs[3]), node_to_number(substrs[4]));
        }
        else if(substrs[0][0] == 'D'){
            tmp_nld = new Diode(node_to_number(substrs[1]), node_to_number(substrs[2]), substrs[3]);

            non_linear_devices.push_back(tmp_nld);
            //include ss equivalent for diode
            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), 0, 0);
        }
        else if(substrs[0][0] == 'Q'){
            tmp_nld = new BJT(node_to_number(substrs[1]), node_to_number(substrs[2]), node_to_number(substrs[3]), substrs[4]);

            non_linear_devices.push_back(tmp_nld);
            //include ss equivalent for BJT
            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), node_to_number(substrs[3]), 0);
        }
        else if(substrs[0][0] == 'M'){
            tmp_nld = new MOSFET(node_to_number(substrs[1]), node_to_number(substrs[2]), node_to_number(substrs[3]), substrs[4]);

            non_linear_devices.push_back(tmp_nld);
            //include ss equivalent for MOSFET
            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), node_to_number(substrs[3]), 0);
        }
        else if(substrs[0] == ".ac"){
            n_ppd = prefix_convertor(substrs[2]);
            f_start = prefix_convertor(substrs[3]);
            f_stop = prefix_convertor(substrs[4]);
        }
        else if(substrs[0] == ".end"){
            infile.close();
        }      

        substrs.clear();
    }

    infile.close();
    
    //for(int i = 0; i < impedance_devices.size(); i++){
        //std::cout << impedance_devices[i]->show_nodeinfo() << std::endl;
        //std::cout << impedance_devices[i]->give_nodeinfo().x << std::endl;
        //std::cout << impedance_devices[i]->give_nodeinfo().y << std::endl;
        //std::cout << "Impedance: " << impedance_devices[i]->get_impedance(1) << std::endl;
        //std::cout << "Conductance: " << impedance_devices[i]->get_conductance(1) << std::endl;
        //std::cout << std::endl;
    //}

    //for(int i = 0; i < sources.size(); i++){
        //std::cout << sources[i]->show_nodeinfo() << std::endl;
        //std::cout << "Source Magnitude: " << sources[i]->get_magnitude() << std::endl;
        //std::cout << "Source Phase: " << sources[i]->get_phase() << std::endl;
        //std::cout << "Source Type: " << sources[i]->get_type() << std::endl;
        //if(sources[i]->get_type() == "VCCS"){
            //std::cout << "Control +: " << sources[i]->give_controlinfo().x << std::endl;
            //std::cout << "Control -: " << sources[i]->give_controlinfo().y << std::endl;
            //std::cout << "Transconductance: " << sources[i]->get_gm() << std::endl;
        //}
        //std::cout << std::endl;
    //}

    //for(int i = 0; i < non_linear_devices.size(); i++){
        //std::cout << non_linear_devices[i]->show_nodeinfo() << std::endl;
        //std::cout << "Model: " << non_linear_devices[i]->get_model() << std::endl;
        //std::cout << std::endl;
    //}

    //std::cout << "Number of points per decade: " << n_ppd << std::endl;
    //std::cout << "Start frequency: " << f_start << " Hz" << std::endl;
    //std::cout << "Stop frequency: " << f_stop << " Hz" << std::endl;
    //std::cout << "Total number of nodes: " << n_max << std::endl;

    MatrixXcd matrixA(n_max,n_max), matrixB(n_max, 1);
    matrixB.setZero();

    std::complex<double> negative(-1, 0), zero(0,0), one(1,0);
    std::complex<double> ACSource(ss_sources[0]->get_magnitude() * cos(ss_sources[0]->get_phase() * M_PI / 180), ss_sources[0]->get_magnitude() * sin(ss_sources[0]->get_phase() * M_PI / 180));
    matrixB(0,0) = ACSource;

    std::vector<double> frequencies;
    std::vector<double> magnitude;
    std::vector<double> phase;

    int n_output, hehe;
    std::cout << "Which node is the output node?" << std::endl;
    std::cin >> n_output;

    for(int n = 0; f < f_stop; n++){
        f = f_start * pow(10, n/n_ppd);
        frequencies.push_back(f);
        omega = 2 * M_PI * f;

        matrixA.setZero();

        for(int i = 0; i < ss_impedance_devices.size(); i++){

            if(ss_impedance_devices[i]->give_nodeinfo().x != 0 && ss_impedance_devices[i]->give_nodeinfo().y != 0){
                matrixA(ss_impedance_devices[i]->give_nodeinfo().x - 1, ss_impedance_devices[i]->give_nodeinfo().y - 1) = matrixA(ss_impedance_devices[i]->give_nodeinfo().x - 1, ss_impedance_devices[i]->give_nodeinfo().y - 1) - ss_impedance_devices[i]->get_conductance(omega);
                matrixA(ss_impedance_devices[i]->give_nodeinfo().y - 1, ss_impedance_devices[i]->give_nodeinfo().x - 1) = matrixA(ss_impedance_devices[i]->give_nodeinfo().y - 1, ss_impedance_devices[i]->give_nodeinfo().x - 1) - ss_impedance_devices[i]->get_conductance(omega);
                matrixA(ss_impedance_devices[i]->give_nodeinfo().x - 1, ss_impedance_devices[i]->give_nodeinfo().x - 1) = matrixA(ss_impedance_devices[i]->give_nodeinfo().x - 1, ss_impedance_devices[i]->give_nodeinfo().x - 1) + ss_impedance_devices[i]->get_conductance(omega);
                matrixA(ss_impedance_devices[i]->give_nodeinfo().y - 1, ss_impedance_devices[i]->give_nodeinfo().y - 1) = matrixA(ss_impedance_devices[i]->give_nodeinfo().y - 1, ss_impedance_devices[i]->give_nodeinfo().y - 1) + ss_impedance_devices[i]->get_conductance(omega);
            }
            else if(ss_impedance_devices[i]->give_nodeinfo().x == 0 && ss_impedance_devices[i]->give_nodeinfo().y != 0){
                matrixA(ss_impedance_devices[i]->give_nodeinfo().y - 1, ss_impedance_devices[i]->give_nodeinfo().y - 1) = matrixA(ss_impedance_devices[i]->give_nodeinfo().y - 1, ss_impedance_devices[i]->give_nodeinfo().y - 1) + ss_impedance_devices[i]->get_conductance(omega);
            }
            else if(ss_impedance_devices[i]->give_nodeinfo().x != 0 && ss_impedance_devices[i]->give_nodeinfo().y == 0){
                matrixA(ss_impedance_devices[i]->give_nodeinfo().x - 1, ss_impedance_devices[i]->give_nodeinfo().x - 1) = matrixA(ss_impedance_devices[i]->give_nodeinfo().x - 1, ss_impedance_devices[i]->give_nodeinfo().x - 1) + ss_impedance_devices[i]->get_conductance(omega);
            }
        
        }

        for(int i = 0; i < n_max; i++){
            matrixA(ss_sources[0]->give_nodeinfo().x - 1,i) = zero;
        }

        matrixA(ss_sources[0]->give_nodeinfo().x - 1,ss_sources[0]->give_nodeinfo().x - 1) = one;

        MatrixXcd matrixX = matrixA.fullPivLu().solve(matrixB);

        magnitude.push_back(return_tf_magnitude(ACSource, matrixX(n_output - 1, 0)));
        phase.push_back(return_tf_phase(ACSource, matrixX(n_output - 1, 0)));

    }

    for(int i = 0; i < magnitude.size(); i++){
        std::cout << magnitude[i] << std::endl;
    }



    //for(int i = 0; i < frequencies.size(); i++){
        //std::cout << frequencies[i] << std::endl;
    //}

    //for(int i = 0; i < phase.size(); i++){
        //std::cout << phase[i] << std::endl;
    //}

}

