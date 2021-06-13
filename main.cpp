#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <complex>
#include "library/Eigen/Dense"
#include <sys/time.h>

using namespace Eigen;

typedef unsigned long long timemarker;

    static timemarker get_timestamp (){
      struct timeval now;
      gettimeofday (&now, NULL);
      return  now.tv_usec + (timemarker)now.tv_sec * 1000000;
    }

//define parsing structure for two terminal components
struct NodePoint{
    int x;
    int y;
};


//define parsing structure for three terminal componenets
struct NodeTri{
    int x;
    int y;
    int z;
};


//class for impedance decives (Resistor, Capacitor, Inductor)
class ImpedanceDevice{
public:
    virtual std::string show_nodeinfo() const = 0;
    
    virtual std::complex<double> get_impedance(double omega) const = 0;

    virtual std::complex<double> get_conductance(double omega) const = 0;

    virtual NodePoint give_nodeinfo() const = 0;
    
    virtual ~ImpedanceDevice() { }
};


//inheritance class for resistor
class Resistor : public ImpedanceDevice{
public:

    Resistor(int n1, int n2, double r) : node1(n1), node2(n2), resistance(r) { }

    std::string show_nodeinfo() const {
        return "Nodal Coordinates: (" + std::to_string(node1) + ", " + std::to_string(node2) + ")";
    }

    std::complex<double> get_impedance(double omega) const {
        std::complex<double> impedance(resistance,0);

        return impedance;
    }

    std::complex<double> get_conductance(double omega) const {
        std::complex<double> conductance(1/resistance,0);

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


//inheritance class for capacittor
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


//inheritance class for inductor
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


//class for sources (DC Voltage, DC Current, AC Voltage, AC Current, Voltage Controlled Current)
class Source{
public:
    virtual std::string show_nodeinfo() const = 0;

    virtual double get_magnitude() const { 
        return 0;
    }

    virtual double get_phase() const { 
        return 0;
    }

    virtual double get_gm() const {
        return 0;
    }

    virtual std::string get_type() const = 0;

    virtual std::string get_source_label() const = 0;

    virtual NodePoint give_nodeinfo() const = 0;

    virtual NodePoint give_controlinfo() const {
        NodePoint C;

        C.x = 0;
        C.y = 0;

        return C;
    }

    virtual ~Source() { }
};


//inheritance class for DC Voltage Source
class DCVSource : public Source{
public:

    DCVSource(int n_p, int n_m, double v, std::string label) : node_plus(n_p), node_minus(n_m), voltage(v), source_label(label){ }

    std::string show_nodeinfo() const {
        return "Nodal Coordinates: (" + std::to_string(node_plus) + ", " + std::to_string(node_minus) + ")";
    }

    double get_magnitude() const {
        return voltage;
    }

    std::string get_type() const {
        return "DC V";
    }

    std::string get_source_label() const {
        return source_label;
    }

    NodePoint give_nodeinfo() const {
        NodePoint N;

        N.x = node_plus;
        N.y = node_minus;

        return N;
    }

private:
    int node_plus;
    int node_minus;
    double voltage;
    std::string source_label;
};


//inheritance class for DC current Source
class DCISource : public Source{
public:

    DCISource(int n_in, int n_out, double i, std::string label) : node_in(n_in), node_out(n_out), current(i), source_label(label){ }

    std::string show_nodeinfo() const {
        return "Nodal Coordinates: (" + std::to_string(node_in) + ", " + std::to_string(node_out) + ")";
    }

    double get_magnitude() const {
        return current;
    }

    std::string get_type() const {
        return "DC I";
    }

    std::string get_source_label() const {
        return source_label;
    }

    NodePoint give_nodeinfo() const {
        NodePoint N;

        N.x = node_in;
        N.y = node_out;

        return N;
    }

private:
    int node_in;
    int node_out;
    double current;
    std::string source_label;
};


//inheritance class for AC Voltage Source
class ACVSource : public Source{
public:

    ACVSource(int n_p, int n_m, double v_m, double v_p, std::string label) : node_plus(n_p), node_minus(n_m), amplitude(v_m), phase(v_p), source_label(label){ }

    std::string show_nodeinfo() const {
        return "Nodal Coordinates: (" + std::to_string(node_plus) + ", " + std::to_string(node_minus) + ")";
    }

    double get_magnitude() const {
        return amplitude;
    }

    double get_phase() const {
        return phase;
    }

    std::string get_type() const {
        return "AC V";
    }

    std::string get_source_label() const {
        return source_label;
    }

    NodePoint give_nodeinfo() const {
        NodePoint N;

        N.x = node_plus;
        N.y = node_minus;

        return N;
    }

private:
    int node_plus;
    int node_minus;
    double amplitude;
    double phase;
    std::string source_label;
};


//inheritance class for AC Current Source
class ACISource : public Source{
public:

    ACISource(int n_in, int n_out, double i_m, double i_p, std::string label) : node_in(n_in), node_out(n_out), amplitude(i_m), phase(i_p), source_label(label){ }

    std::string show_nodeinfo() const {
        return "Nodal Coordinates: (" + std::to_string(node_in) + ", " + std::to_string(node_out) + ")";
    }

    double get_magnitude() const {
        return amplitude;
    }

    double get_phase() const {
        return phase;
    }

    std::string get_type() const {
        return "AC I";
    }

    std::string get_source_label() const {
        return source_label;
    }

    NodePoint give_nodeinfo() const {
        NodePoint N;

        N.x = node_in;
        N.y = node_out;

        return N;
    }

private:
    int node_in;
    int node_out;
    double amplitude;
    double phase;
    std::string source_label;
};


//inheritance class for Voltage Controlled Current Source
class VCCSource : public Source{
public:

    VCCSource(int n_p, int n_m, int c_p, int c_m, double g_m, std::string label) : node_plus(n_p), node_minus(n_m), control_plus(c_p), control_minus(c_m), transconductance(g_m), source_label(label){ }

    std::string show_nodeinfo() const {
        return "Nodal Coordinates: (" + std::to_string(node_plus) + ", " + std::to_string(node_minus) + ")";
    }

    double get_gm() const {
        return transconductance;
    }

    std::string get_type() const {
        return "VCCS";
    }

    std::string get_source_label() const {
        return source_label;
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
    std::string source_label;
};


//class for non-linear devices (diode, bjt, mosfet)
class NonLinearDevice{
public:
    virtual std::string show_nodeinfo() const = 0;

    virtual NodePoint give_binodeinfo() const {
        NodePoint N;

        N.x = 0;
        N.y = 0;

        return N;
    }

    virtual NodeTri give_trinodeinfo() const {
        NodeTri N;

        N.x = 0;
        N.y = 0;
        N.z = 0;

        return N;
    }

    virtual std::string get_model() const = 0;

    virtual ~NonLinearDevice() { }
};


//inheritance class for diode
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

    std::string get_model() const {
        return model;
    }

private:
    int node_an;
    int node_cat;
    std::string model;
};


//inheritance class for bjt
class BJT: public NonLinearDevice{
public: 

    BJT(int n_c, int n_b, int n_e, std::string mod) : node_c(n_c), node_b(n_b), node_e(n_e), model(mod) { }

    std::string show_nodeinfo() const {
        return "Nodal Coordinates: (" + std::to_string(node_c) + ", " + std::to_string(node_b) + ", " + std::to_string(node_e) +")";
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


//inheritance class for mosfet
class MOSFET: public NonLinearDevice{
public: 

    MOSFET(int n_d, int n_g, int n_s, std::string mod) : node_d(n_d), node_g(n_g), node_s(n_s), model(mod) { }

    std::string show_nodeinfo() const {
        return "Nodal Coordinates: (" + std::to_string(node_d) + ", " + std::to_string(node_g) + ", "  + std::to_string(node_s) + ")";
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


//ignores non-numeric characters in value part and extract the magnitude
double extract_double(const std::string& label){
    std::string double_string;

    for(int i = 0; i < label.size(); i++){
        if(std::isdigit(label[i]) || label[i] == '.' || label[i] == '-'){
         double_string.push_back(label[i]);
        }
    }

    return std::stod(double_string);
}


//change input string node into number
int node_to_number(const std::string& node){
    int node_number;
    double node_double;

    node_double = extract_double(node);
    node_number = (int) node_double;

    return node_number;
}


//convert prefix into value 
double prefix_convertor(const std::string& value){
    std::string prefix;
    double num_quantity;
    //loop needed as not all prefixes are 1 character long
    for(int i = 0; i < value.size(); i++){
        if(std::isdigit(value[i]) == false && value[i] != '.' && value[i] != '-' && value[i] != '(' && value[i] != 'A' && value[i] != 'C'){
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


//find the number of nodes in the circuit
// can limit to three node inputs because VCCS control nodes are presumably connected to other components present in the circuit that will be considered
int max_node_number(int node_max, const int& node1, const int& node2, const int& node3){
    
    if((node1 > node_max) && (node1 >= node2) && (node1 >= node3)){
        return node1;
    }
    else if((node2 > node_max) && (node2 >= node1) && (node2 >= node3)){
        return node2;
    }
    else if((node3 > node_max) && (node3 >= node1) && (node3 >= node2)){
        return node3;
    }

    return node_max;
}


bool detect_input_source(const std::vector<Source*>& sources, const std::string& inputsource_label){

    for(int i = 0; i < sources.size(); i++){
        if(sources[i]->get_source_label() == inputsource_label){
            return true;
        }
    }

    return false;
}


//calculate the magnitude of the transfer function
double return_tf_magnitude(const std::complex<double>& source, const std::complex<double>& output_node){
    double gain;

    gain = std::abs(output_node/source);

    return gain;
}


//calculate the phase of the transfer function
double return_tf_phase(const std::complex<double>& source, const std::complex<double>& output_node){
    double phase_change;

    phase_change = (std::arg(output_node) - std::arg(source)) * 180 / M_PI;

    return phase_change;
}


//construct conductanc matrix only with conductances (ignoring source rows)
MatrixXcd cons_conductance_matrix(MatrixXcd A, const std::vector<ImpedanceDevice*>& impedances, const double& omega){

    for(int i = 0; i < impedances.size(); i++){

        if(impedances[i]->give_nodeinfo().x != 0 && impedances[i]->give_nodeinfo().y != 0){
            A(impedances[i]->give_nodeinfo().x - 1, impedances[i]->give_nodeinfo().y - 1) = A(impedances[i]->give_nodeinfo().x - 1, impedances[i]->give_nodeinfo().y - 1) - impedances[i]->get_conductance(omega);
            A(impedances[i]->give_nodeinfo().y - 1, impedances[i]->give_nodeinfo().x - 1) = A(impedances[i]->give_nodeinfo().y - 1, impedances[i]->give_nodeinfo().x - 1) - impedances[i]->get_conductance(omega);
            A(impedances[i]->give_nodeinfo().x - 1, impedances[i]->give_nodeinfo().x - 1) = A(impedances[i]->give_nodeinfo().x - 1, impedances[i]->give_nodeinfo().x - 1) + impedances[i]->get_conductance(omega);
            A(impedances[i]->give_nodeinfo().y - 1, impedances[i]->give_nodeinfo().y - 1) = A(impedances[i]->give_nodeinfo().y - 1, impedances[i]->give_nodeinfo().y - 1) + impedances[i]->get_conductance(omega);
        }
        else if(impedances[i]->give_nodeinfo().x == 0 && impedances[i]->give_nodeinfo().y != 0){
            A(impedances[i]->give_nodeinfo().y - 1, impedances[i]->give_nodeinfo().y - 1) = A(impedances[i]->give_nodeinfo().y - 1, impedances[i]->give_nodeinfo().y - 1) + impedances[i]->get_conductance(omega);
        }
        else if(impedances[i]->give_nodeinfo().x != 0 && impedances[i]->give_nodeinfo().y == 0){
            A(impedances[i]->give_nodeinfo().x - 1, impedances[i]->give_nodeinfo().x - 1) = A(impedances[i]->give_nodeinfo().x - 1, impedances[i]->give_nodeinfo().x - 1) + impedances[i]->get_conductance(omega);
        }
        
    }

    return A;
}


//considers both grounded and floating current sources
void isource_analysis(MatrixXcd& Bref, Source* source_i, const std::complex<double>& isource){
    std::complex<double> negative(-1,0);

    if(source_i->give_nodeinfo().x != 0){
        Bref(source_i->give_nodeinfo().x - 1,0) = negative * isource;
        //negative due to orientation of current source
    }

    if(source_i->give_nodeinfo().y != 0){
        Bref(source_i->give_nodeinfo().y - 1,0) = isource;
    }
}


//considers both grounded and floating VCCS
void vccsource_analysis(MatrixXcd& A, Source* source_vcc, const std::complex<double>& gm){

    if(source_vcc->give_nodeinfo().x != 0){
        A(source_vcc->give_nodeinfo().x - 1, source_vcc->give_controlinfo().x - 1) = A(source_vcc->give_nodeinfo().x - 1, source_vcc->give_controlinfo().x - 1) + gm;
        A(source_vcc->give_nodeinfo().x - 1, source_vcc->give_controlinfo().y - 1) = A(source_vcc->give_nodeinfo().x - 1, source_vcc->give_controlinfo().y - 1) - gm;
    }

    if(source_vcc->give_nodeinfo().y != 0){
        A(source_vcc->give_nodeinfo().y - 1, source_vcc->give_controlinfo().x - 1) = A(source_vcc->give_nodeinfo().y - 1, source_vcc->give_controlinfo().x - 1) - gm;
        A(source_vcc->give_nodeinfo().y - 1, source_vcc->give_controlinfo().y - 1) = A(source_vcc->give_nodeinfo().y - 1, source_vcc->give_controlinfo().y - 1) + gm;
    }

}


//forms supernode row by adding the rows of the 2 nodes that form the supernode
void fvsource_analysis1(MatrixXcd& A, MatrixXcd& B, const MatrixXcd& G, const MatrixXcd& Bref, Source* source_fv, const int& n_max){

    for(int k = 0; k < n_max; k++){
        A(source_fv->give_nodeinfo().y - 1, k) = G(source_fv->give_nodeinfo().y - 1, k) + G(source_fv->give_nodeinfo().x - 1, k);
    }

    B(source_fv->give_nodeinfo().y - 1, 0) = Bref(source_fv->give_nodeinfo().y - 1, 0) + Bref(source_fv->give_nodeinfo().x - 1, 0);
}


//inserts 0,1,-1 in A and Vsrc in B
void fvsource_analysis2(MatrixXcd& A, MatrixXcd& B, Source* source_fv, const std::complex<double>& fvsource, const int& n_max){
    std::complex<double> zero(0,0), one(1,0), negative(-1,0);

    //sets row representing floating source to all zero first
    for(int k = 0; k < n_max; k++){
        A(source_fv->give_nodeinfo().x - 1,k) = zero;
    }
    //inserts 1 and -1 into row representing voltage source
    A(source_fv->give_nodeinfo().x - 1,source_fv->give_nodeinfo().x - 1) = one;
    A(source_fv->give_nodeinfo().x - 1,source_fv->give_nodeinfo().y - 1) = negative;
    //sets correct entry of B matrix to represent the source
    B(source_fv->give_nodeinfo().x - 1,0) = fvsource;
}


//inserts 0,1 in A and Vsrc or -Vsrc in B
void gvsource_analysis(MatrixXcd& A, MatrixXcd& B, Source* source_gv, const std::complex<double>& gvsource, const int& n_max){
    std::complex<double> zero(0,0), one(1,0), negative(-1,0);

    if(source_gv->give_nodeinfo().x != 0){
        B(source_gv->give_nodeinfo().x - 1,0) = gvsource;

        for(int k = 0; k < n_max; k++){
            A(source_gv->give_nodeinfo().x - 1,k) = zero;
        }

        A(source_gv->give_nodeinfo().x - 1,source_gv->give_nodeinfo().x - 1) = one;
    }
    else{
        B(source_gv->give_nodeinfo().y - 1,0) = negative * gvsource;
        //account for polarity of voltage source
        for(int k = 0; k < n_max; k++){
            A(source_gv->give_nodeinfo().y - 1,k) = zero;
        }

        A(source_gv->give_nodeinfo().y - 1,source_gv->give_nodeinfo().y - 1) = one;
    }
}


//resets matrices after every iteration of loop for DC and AC analysis
void reset_matrices(MatrixXcd& A, MatrixXcd& B, MatrixXcd& Bref, MatrixXcd& G, MatrixXcd& X, std::vector<ImpedanceDevice*> impedance_devices, double omega){
    A.setZero();
    A = cons_conductance_matrix(A, impedance_devices, omega);
    B.setZero();
    Bref.setZero();
    G.setZero();
    G = cons_conductance_matrix(G, impedance_devices, omega);
    X.setZero();
}


int main(){
    timemarker t0 = get_timestamp();
    std::ifstream infile; 
    std::string input_file_name;

    std::cout << "What is the name of the netlist input file?" << std::endl;
    timemarker t1 = get_timestamp();
    std::cin >> input_file_name;
    timemarker t2 = get_timestamp();
    infile.open(input_file_name);
 
    if(!infile.is_open()){
        std::cout << "error opening file" << std::endl;
        return EXIT_FAILURE;
    }
 
    std::string component;
    std::vector<std::string> substrs;
    //can we delete impedance_devices and sources? not needed for analysis portion only needed to test parse
    std::vector<ImpedanceDevice*> ss_impedance_devices, dc_impedance_devices; 
    ImpedanceDevice* tmp_id;
    ImpedanceDevice* tmp_id2;

    std::vector<Source*> ss_sources, dc_sources; 
    Source* tmp_s;
    Source* tmp_s2;

    //nonlinear values
    std::vector<NonLinearDevice*> non_linear_devices;
    NonLinearDevice* tmp_nld;
    
    // frequency step parameters, we assume ac analysis always done in decades
    double f_start, f_stop, n_ppd, f=0, omega = 0;
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
            
            ss_impedance_devices.push_back(tmp_id);
            dc_impedance_devices.push_back(tmp_id);//new

            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), 0);
        }
        else if(substrs[0][0] == 'C'){
            tmp_id = new Capacitor(node_to_number(substrs[1]), node_to_number(substrs[2]), prefix_convertor(substrs[3]));
            
            ss_impedance_devices.push_back(tmp_id);

            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), 0);
        }
        else if(substrs[0][0] == 'L'){
            tmp_id = new Inductor(node_to_number(substrs[1]), node_to_number(substrs[2]), prefix_convertor(substrs[3]));

            tmp_id2 = new Resistor(node_to_number(substrs[1]), node_to_number(substrs[2]), 0.001);

            ss_impedance_devices.push_back(tmp_id);
            dc_impedance_devices.push_back(tmp_id2);

            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), 0);
        }
        else if(substrs[0][0] == 'V' && substrs[3][0] != 'A'){
            tmp_s = new DCVSource(node_to_number(substrs[1]), node_to_number(substrs[2]), prefix_convertor(substrs[3]), substrs[0]);
            tmp_id = new Resistor(node_to_number(substrs[1]), node_to_number(substrs[2]), 0.001);
          
            //ss equivalent of DC Voltage Source is short circuit, represented by resistor of 1m
            dc_sources.push_back(tmp_s);//new
            ss_impedance_devices.push_back(tmp_id);

            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), 0);
        }
        else if(substrs[0][0] == 'V' && substrs[3][0] == 'A'){
            tmp_s = new ACVSource(node_to_number(substrs[1]), node_to_number(substrs[2]), prefix_convertor(substrs[3]), extract_double(substrs[4]), substrs[0]);
            tmp_id = new Resistor(node_to_number(substrs[1]), node_to_number(substrs[2]), 0.001);

            ss_sources.push_back(tmp_s);
            dc_impedance_devices.push_back(tmp_id);//new

            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), 0);
        }
        else if(substrs[0][0] == 'I' && substrs[3][0] != 'A'){
            tmp_s = new DCISource(node_to_number(substrs[1]), node_to_number(substrs[2]), prefix_convertor(substrs[3]), substrs[0]);
            
            dc_sources.push_back(tmp_s);//new

            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), 0);
        }
        else if(substrs[0][0] == 'I' && substrs[3][0] == 'A'){
            tmp_s = new ACISource(node_to_number(substrs[1]), node_to_number(substrs[2]), prefix_convertor(substrs[3]), extract_double(substrs[4]), substrs[0]);

            ss_sources.push_back(tmp_s);

            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), 0);
        }
        else if(substrs[0][0] == 'G'){
            tmp_s = new VCCSource(node_to_number(substrs[1]), node_to_number(substrs[2]), node_to_number(substrs[3]), node_to_number(substrs[4]), prefix_convertor(substrs[5]), substrs[0]);

            ss_sources.push_back(tmp_s);

            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), 0);
        }
        else if(substrs[0][0] == 'D'){
            tmp_nld = new Diode(node_to_number(substrs[1]), node_to_number(substrs[2]), substrs[3]);
            
            non_linear_devices.push_back(tmp_nld);

            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), 0);
        }
        else if(substrs[0][0] == 'Q'){
            tmp_nld = new BJT(node_to_number(substrs[1]), node_to_number(substrs[2]), node_to_number(substrs[3]), substrs[4]);

            non_linear_devices.push_back(tmp_nld);
            
            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), node_to_number(substrs[3]));
        }
        else if(substrs[0][0] == 'M'){
            tmp_nld = new MOSFET(node_to_number(substrs[1]), node_to_number(substrs[2]), node_to_number(substrs[3]), substrs[4]);

            non_linear_devices.push_back(tmp_nld);
            
            n_max = max_node_number(n_max, node_to_number(substrs[1]), node_to_number(substrs[2]), node_to_number(substrs[3]));
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

    MatrixXcd matrixA(n_max,n_max), matrixG(n_max, n_max), matrixB(n_max, 1), matrixBref(n_max, 1), matrixX(n_max, 1);

    std::complex<double> InputSource(0,0);

    int n_output;
    std::string s_input;
    std::cout << "Which node is the output node?" << std::endl;
    timemarker t3 = get_timestamp();
    std::cin >> n_output;

    if((n_output > n_max) || (n_output <= 0)){
        std::cout << "error, invalid output node" << std::endl;

        return EXIT_FAILURE;
    }

    timemarker t4 = get_timestamp();
    std::cout << "Which source is the input source?" << std::endl;
    timemarker t5 = get_timestamp();
    std::cin >> s_input;

    if(detect_input_source(ss_sources, s_input) == false){
        std::cout << "error, nominated input source does not exist" << std::endl;

        return EXIT_FAILURE;
    }

    timemarker t6 = get_timestamp();

    for(int i=0; i<non_linear_devices.size(); i++){
        double Geq, Ieq, Vd = 0.7, Id, Is_diode = 1 * pow(10, -14), Is_bjt = 1 * pow(10,-16), Vt = 25.865 * pow(10, -3), V1 = 0, V2 = 0, Vdlast = 1, beta = 100, Kp = 2 * pow(10,-5);
        
        while(std::abs(Vdlast - Vd)>=0.00000001){

            Vdlast = Vd;
            
            if(non_linear_devices[i]->get_model() == "D"){ 
                Id = Is_diode * (exp(Vd/Vt)-1);
                Geq = Is_diode/Vt * exp(Vd/Vt);
                Ieq = Id - Geq* Vd;
                tmp_s= new DCISource(non_linear_devices[i]->give_binodeinfo().x, non_linear_devices[0]->give_binodeinfo().y, Ieq, "NA" );
                tmp_id = new Resistor(non_linear_devices[i]->give_binodeinfo().x, non_linear_devices[0]->give_binodeinfo().y, 1/Geq);

                dc_sources.push_back(tmp_s);
                dc_impedance_devices.push_back(tmp_id);
            }
            else if(non_linear_devices[i]->get_model() == "NPN"){
                Id = Is_bjt/beta * (exp(Vd/Vt)-1);
                Geq = Is_bjt/Vt * exp(Vd/Vt);
                Ieq = Id - Geq* Vd;
                tmp_s = new DCISource(non_linear_devices[i]->give_trinodeinfo().y, non_linear_devices[0]->give_trinodeinfo().z, Ieq, "NA" );
                tmp_s2 = new DCISource(non_linear_devices[i]->give_trinodeinfo().x, non_linear_devices[0]->give_trinodeinfo().z, Id*beta, "NA" );
                tmp_id = new Resistor(non_linear_devices[i]->give_trinodeinfo().y, non_linear_devices[0]->give_trinodeinfo().z, 1/Geq);
                
                dc_sources.push_back(tmp_s);
                dc_sources.push_back(tmp_s2);
                dc_impedance_devices.push_back(tmp_id);
            }
            else if(non_linear_devices[i]->get_model() == "PNP"){
                Id = Is_bjt/beta * (exp(Vd/Vt)-1);
                Geq = Is_bjt/Vt * exp(Vd/Vt);
                Ieq = Id - Geq* Vd;
                tmp_s = new DCISource(non_linear_devices[i]->give_trinodeinfo().z, non_linear_devices[i]->give_trinodeinfo().y, Ieq, "NA");
                tmp_s2 = new DCISource(non_linear_devices[i]->give_trinodeinfo().z, non_linear_devices[i]->give_trinodeinfo().x, Id*beta, "NA");
                tmp_id = new Resistor(non_linear_devices[i]->give_trinodeinfo().z, non_linear_devices[i]->give_trinodeinfo().y, 1/Geq);
                
                dc_sources.push_back(tmp_s);
                dc_sources.push_back(tmp_s2);
                dc_impedance_devices.push_back(tmp_id);
            }
            else if(non_linear_devices[i]->get_model() == "NMOS"){
                Id = 0.5 * Kp * pow(Vd, 2);
                tmp_s = new DCISource(non_linear_devices[i]->give_trinodeinfo().x, non_linear_devices[i]->give_trinodeinfo().z, Id, "NA");

                dc_sources.push_back(tmp_s);
            }
            else if(non_linear_devices[i]->get_model() == "PMOS"){
                Id = 0.5 * Kp * pow(Vd, 2);
                tmp_s = new DCISource(non_linear_devices[i]->give_trinodeinfo().z, non_linear_devices[i]->give_trinodeinfo().x, Id, "NA");

                dc_sources.push_back(tmp_s);
            }

            reset_matrices(matrixA, matrixB, matrixBref, matrixG, matrixX, dc_impedance_devices, omega);

            for(int j = 0; j < dc_sources.size(); j++){
                std::complex<double> DCSource(dc_sources[j]->get_magnitude(), 0);

                matrixBref.setZero();

                if(dc_sources[j]->get_type() == "DC I"){
                    isource_analysis(matrixBref, dc_sources[j], DCSource);
                }

                matrixB = matrixB + matrixBref;                
            }

            matrixBref = matrixB;

            for(int j = 0; j < dc_sources.size(); j++){

                if(dc_sources[j]->get_type() == "VCCS"){
                    std::complex<double> Gm(dc_sources[j]->get_gm(), 0);

                    vccsource_analysis(matrixA, dc_sources[j], Gm);
                }
            }

            for(int j = 0; j < dc_sources.size(); j++){

                if(dc_sources[j]->get_type() == "DC V" && dc_sources[j]->give_nodeinfo().x != 0 && dc_sources[j]->give_nodeinfo().y != 0){
                    fvsource_analysis1(matrixA, matrixB, matrixG, matrixBref, dc_sources[j], n_max);
                }
            }

            for(int j = 0; j < dc_sources.size(); j++){
                std::complex<double> DCSource(dc_sources[j]->get_magnitude(), 0);

                if(dc_sources[j]->get_type() == "DC V" && dc_sources[j]->give_nodeinfo().x != 0 && dc_sources[j]->give_nodeinfo().y != 0){
                    fvsource_analysis2(matrixA, matrixB, dc_sources[j], DCSource, n_max);
                }
            }

            for(int j = 0; j < dc_sources.size(); j++){
                std::complex<double> DCSource(dc_sources[j]->get_magnitude(), 0);

                if(dc_sources[j]->get_type() == "DC V" && (dc_sources[j]->give_nodeinfo().x == 0 || dc_sources[j]->give_nodeinfo().y == 0)){
                    gvsource_analysis(matrixA, matrixB, dc_sources[j], DCSource, n_max);
                }
            }

            matrixX = matrixA.fullPivLu().solve(matrixB);
            
            if(non_linear_devices[i]->get_model()=="D"){
                if (non_linear_devices[i]->give_binodeinfo().x == 0){
                    V1 = 0;
                    V2 = std::real(matrixX(non_linear_devices[i]->give_binodeinfo().y -1,0));
                }
                else if (non_linear_devices[i]->give_binodeinfo().y == 0){
                    V1 = std::real(matrixX(non_linear_devices[i]->give_binodeinfo().x -1,0));
                    V2 =0;
                }
                else{
                    V1 = std::real(matrixX(non_linear_devices[i]->give_binodeinfo().x -1,0));
                    V2 = std::real(matrixX(non_linear_devices[i]->give_binodeinfo().y -1,0));
                }
                //conditional for grounding needed, if not the matrix coordinate will be negative
                dc_sources.pop_back();
                dc_impedance_devices.pop_back();
            }
            else if(non_linear_devices[i]->get_model() == "NPN"){
                if (non_linear_devices[i]->give_trinodeinfo().y == 0){
                    V1 = 0;
                    V2 = std::real(matrixX(non_linear_devices[i]->give_trinodeinfo().z -1,0));
                }
                else if (non_linear_devices[i]->give_trinodeinfo().z == 0){
                    V1 = std::real(matrixX(non_linear_devices[i]->give_trinodeinfo().y -1,0));
                    V2 = 0;
                }
                else{
                    V1 = std::real(matrixX(non_linear_devices[i]->give_trinodeinfo().y -1,0));
                    V2 = std::real(matrixX(non_linear_devices[i]->give_trinodeinfo().z -1,0));
                }

                dc_sources.pop_back();
                dc_sources.pop_back();
                dc_impedance_devices.pop_back(); 
            }
            else if(non_linear_devices[i]->get_model() == "PNP"){
                if(non_linear_devices[i]->give_trinodeinfo().y == 0){
                    V1 = std::real(matrixX(non_linear_devices[i]->give_trinodeinfo().z -1,0));
                    V2 = 0;
                }
                else if(non_linear_devices[i]->give_trinodeinfo().z == 0){
                    V1 = 0;
                    V2 = std::real(matrixX(non_linear_devices[i]->give_trinodeinfo().y -1,0));
                }
                else{
                    V1 = std::real(matrixX(non_linear_devices[i]->give_trinodeinfo().z -1,0));
                    V2 = std::real(matrixX(non_linear_devices[i]->give_trinodeinfo().y -1,0));
                }

                dc_sources.pop_back();
                dc_sources.pop_back();
                dc_impedance_devices.pop_back(); 
            }
            else if(non_linear_devices[i]->get_model() == "NMOS"){
                if (non_linear_devices[i]->give_trinodeinfo().y == 0){
                    V1 = 0;
                    V2 = std::real(matrixX(non_linear_devices[i]->give_trinodeinfo().z -1,0));
                }
                else if (non_linear_devices[i]->give_trinodeinfo().z == 0){
                    V1 = std::real(matrixX(non_linear_devices[i]->give_trinodeinfo().y -1,0));
                    V2 = 0;
                }
                else{
                    V1 = std::real(matrixX(non_linear_devices[i]->give_trinodeinfo().y -1,0));
                    V2 = std::real(matrixX(non_linear_devices[i]->give_trinodeinfo().z -1,0));
                }

                dc_sources.pop_back();
            }
            else if(non_linear_devices[i]->get_model() == "PMOS"){
                if (non_linear_devices[i]->give_trinodeinfo().y == 0){
                    V1 = std::real(matrixX(non_linear_devices[i]->give_trinodeinfo().z -1,0));
                    V2 = 0;
                }
                else if (non_linear_devices[i]->give_trinodeinfo().z == 0){
                    V1 = 0;
                    V2 = std::real(matrixX(non_linear_devices[i]->give_trinodeinfo().y -1,0));
                }
                else{
                    V1 = std::real(matrixX(non_linear_devices[i]->give_trinodeinfo().z -1,0));
                    V2 = std::real(matrixX(non_linear_devices[i]->give_trinodeinfo().y -1,0));
                }

                dc_sources.pop_back();
            }
            
            Vd = V1-V2;
        } 

        if(non_linear_devices[i]->get_model() == "D"){
            tmp_id = new Resistor(non_linear_devices[i]->give_binodeinfo().x, non_linear_devices[i]->give_binodeinfo().y,Vt/Id );
            
            ss_impedance_devices.push_back(tmp_id);
        }
        else if(non_linear_devices[i]->get_model() == "NPN"){
            tmp_id = new Resistor(non_linear_devices[i]->give_trinodeinfo().y, non_linear_devices[i]->give_trinodeinfo().z, Vt/Id);
            //No ro beacause early voltage is assumed infinite, which ro = VA/ IC
            tmp_s = new VCCSource(non_linear_devices[i]->give_trinodeinfo().x, non_linear_devices[i]->give_trinodeinfo().z, non_linear_devices[i]->give_trinodeinfo().y, non_linear_devices[i]->give_trinodeinfo().z, Id*beta/Vt, "Gnpn"); 
            //VCCS with Gm value of Id*beta/Vt
            ss_impedance_devices.push_back(tmp_id);
            ss_sources.push_back(tmp_s);
        }
        else if(non_linear_devices[i]->get_model() == "PNP"){
            tmp_id = new Resistor(non_linear_devices[i]->give_trinodeinfo().z, non_linear_devices[i]->give_trinodeinfo().y, Vt/Id);
            //No ro beacause early voltage is assumed infinite, which ro = VA/ IC
            tmp_s = new VCCSource(non_linear_devices[i]->give_trinodeinfo().z, non_linear_devices[i]->give_trinodeinfo().x, non_linear_devices[i]->give_trinodeinfo().z, non_linear_devices[i]->give_trinodeinfo().y, Id*beta/Vt, "Gpnp"); 
            //VCCS with Gm value of Id*beta/Vt
            ss_impedance_devices.push_back(tmp_id);
            ss_sources.push_back(tmp_s);
        }
        else if(non_linear_devices[i]->get_model() == "NMOS"){
            tmp_s = new VCCSource(non_linear_devices[i]->give_trinodeinfo().x, non_linear_devices[i]->give_trinodeinfo().z, non_linear_devices[i]->give_trinodeinfo().y, non_linear_devices[i]->give_trinodeinfo().z, 2 * std::sqrt(0.5 * Kp * Id), "Gnmos");

            ss_sources.push_back(tmp_s);
        }
        else if(non_linear_devices[i]->get_model() == "PMOS"){
            tmp_s = new VCCSource(non_linear_devices[i]->give_trinodeinfo().z, non_linear_devices[i]->give_trinodeinfo().x, non_linear_devices[i]->give_trinodeinfo().z, non_linear_devices[i]->give_trinodeinfo().y, 2 * std::sqrt(0.5 * Kp * Id), "Gpmos");

            ss_sources.push_back(tmp_s);
        }
    }

    for(int i = 0; i < ss_sources.size(); i++){
        if(ss_sources[i]->get_source_label() == s_input){
            std::complex<double> input_s(ss_sources[i]->get_magnitude() * cos(ss_sources[i]->get_phase() * M_PI / 180), ss_sources[i]->get_magnitude() * sin(ss_sources[i]->get_phase() * M_PI / 180));

            InputSource = input_s;
        }
    }

    std::ofstream outfile;

    outfile.open("output.txt");

    if(!outfile.is_open()){
        std::cout << "error opening file" << std::endl;
        return EXIT_FAILURE;
    }

    outfile << "frequency,magnitude,phase" << std::endl;

    double current_phase, last_phase;
    for(int n = 0; f < f_stop; n++){
        f = f_start * pow(10, n/n_ppd);
        omega = 2 * M_PI * f;

        reset_matrices(matrixA, matrixB, matrixBref, matrixG, matrixX, ss_impedance_devices, omega);

        for(int j = 0; j < ss_sources.size(); j++){
            std::complex<double> ACSource(ss_sources[j]->get_magnitude() * cos(ss_sources[j]->get_phase() * M_PI / 180), ss_sources[j]->get_magnitude() * sin(ss_sources[j]->get_phase() * M_PI / 180));

            matrixBref.setZero();

            if(ss_sources[j]->get_type() == "AC I"){
                isource_analysis(matrixBref, ss_sources[j], ACSource);
            }

            matrixB = matrixB + matrixBref;                
        }

        matrixBref = matrixB;

        for(int j = 0; j < ss_sources.size(); j++){
            std::complex<double> ACSource(ss_sources[j]->get_magnitude() * cos(ss_sources[j]->get_phase() * M_PI / 180), ss_sources[j]->get_magnitude() * sin(ss_sources[j]->get_phase() * M_PI / 180));

            if(ss_sources[j]->get_type() == "VCCS"){
                std::complex<double> Gm(ss_sources[j]->get_gm(), 0);

                vccsource_analysis(matrixA, ss_sources[j], Gm);
            }
        }

        for(int j = 0; j < ss_sources.size(); j++){
            std::complex<double> ACSource(ss_sources[j]->get_magnitude() * cos(ss_sources[j]->get_phase() * M_PI / 180), ss_sources[j]->get_magnitude() * sin(ss_sources[j]->get_phase() * M_PI / 180));

            if(ss_sources[j]->get_type() == "AC V" && ss_sources[j]->give_nodeinfo().x != 0 && ss_sources[j]->give_nodeinfo().y != 0){
                fvsource_analysis1(matrixA, matrixB, matrixG, matrixBref, ss_sources[j], n_max);
            }
        }

        for(int j = 0; j < ss_sources.size(); j++){
            std::complex<double> ACSource(ss_sources[j]->get_magnitude() * cos(ss_sources[j]->get_phase() * M_PI / 180), ss_sources[j]->get_magnitude() * sin(ss_sources[j]->get_phase() * M_PI / 180));

            if(ss_sources[j]->get_type() == "AC V" && ss_sources[j]->give_nodeinfo().x != 0 && ss_sources[j]->give_nodeinfo().y != 0){
                fvsource_analysis2(matrixA, matrixB, ss_sources[j], ACSource, n_max);
            }
        }

        for(int j = 0; j < ss_sources.size(); j++){
            std::complex<double> ACSource(ss_sources[j]->get_magnitude() * cos(ss_sources[j]->get_phase() * M_PI / 180), ss_sources[j]->get_magnitude() * sin(ss_sources[j]->get_phase() * M_PI / 180));

            if(ss_sources[j]->get_type() == "AC V" && (ss_sources[j]->give_nodeinfo().x == 0 || ss_sources[j]->give_nodeinfo().y == 0)){
                gvsource_analysis(matrixA, matrixB, ss_sources[j], ACSource, n_max);
            }
        }

        matrixX = matrixA.fullPivLu().solve(matrixB);

        current_phase = return_tf_phase(InputSource, matrixX(n_output - 1, 0));

        if(current_phase - last_phase > 180){
            if(current_phase > 0){
                current_phase = current_phase - 360;
            }
            else{
                current_phase = current_phase + 360;
            }
        }

        outfile << f << "," << 20 * log10(return_tf_magnitude(InputSource, matrixX(n_output - 1, 0))) << "," << current_phase << std::endl;

        last_phase = current_phase;
    }

    timemarker t7 = get_timestamp();
    std::cout << "Time taken by function: "<< (t1-t0) + (t3-t2) + (t5-t4) + (t7-t6) << " microseconds" << std::endl;
}

