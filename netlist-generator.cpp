#include <iostream>
#include <fstream>
#include <cstdlib>

int main(){
    std::ofstream testfile;
    
    double value = 10000;
    
    for (int s=0; s<65; s++){

        std::string name="netlist_timetest"+ std::to_string(s+1)+".txt";

        testfile.open(name);
        testfile<<"V1" << " " << "N001"<< " " << "0" << " " << "AC(5 0)"<<std::endl;
        testfile<<"L1" << " "<<"N002"<<" "<<"N001"<<" "<<"100m"<<std::endl;
        testfile<<"D1"<<" "<<"N002"<<' '<<"N003"<<" "<<"D"<<std::endl;
        testfile<<"C1"<<" "<<"N003"<<" " <<"0"<<" "<<"10"<<"\u00B5"<<std::endl;
        if (s==0){
            testfile<<"R1"<<" "<<"N002"<<" "<<"0"<<" "<<"10000"<<std::endl;
        }
        
        for (int q=0; q<4*s; q++){
            
            
            if ((q+1) ==4*s){
                testfile<<"R"<<std::to_string(1+q)<<" "<<"N00"<<std::to_string(3+q)<<" 0 "<<value/(4*s)<<std::endl;
            }
            else if (q==0){
                testfile<<"R1 N002 N004 "<<value/(4*s)<<std::endl;
            }
            else if (q==1){
                testfile<<"R2 N004 N005 "<<value/(4*s)<<std::endl;
            }
            else{
                testfile<<"R"<<std::to_string(1+q)<<" "<<"N00"<<std::to_string(3+q)<<" "<<"N00"<<std::to_string(4+q)<<" "<<value/(4*s)<<std::endl; 
            }
             
        }
        testfile<<".ac "<<"dec "<<"10 "<<"0.01 "<<"10k"<<std::endl;
        testfile<<".end"<<std::endl;       
    }
    testfile.close();
}