#include <iostream>
#include <fstream>
#include <cstdlib>

int main(){
    std::ofstream outfile;
    
    double value = 10000;
    
    for (int i=0; i<65; i++){

        std::string name="timetest"+ std::to_string(i+1)+".txt";

        outfile.open(name);
        outfile<<"V1" << " " << "N001"<< " " << "0" << " " << "AC(5 0)"<<std::endl;
        outfile<<"L1" << " "<<"N002"<<" "<<"N001"<<" "<<"100m"<<std::endl;
        outfile<<"D1"<<" "<<"N002"<<' '<<"N003"<<" "<<"D"<<std::endl;
        outfile<<"C1"<<" "<<"N003"<<" " <<"0"<<" "<<"10"<<"\u00B5"<<std::endl;
        if (i==0){
            outfile<<"R1"<<" "<<"N002"<<" "<<"0"<<" "<<"10000"<<std::endl;
        }
        
        for (int j=0; j<4*i; j++){
            outfile<<"R"<<std::to_string(1+j)<<" "<<"N00"<<std::to_string(2+j);
            if ((j+1) ==4*i){
                outfile<<" "<<"0"<<" "<<value/(4*i)<<std::endl;
            }
            else{
                outfile<<" "<<"N00"<<std::to_string(3+j)<<" "<<value/(4*i)<<std::endl; 
            }
             
        }
        outfile.close();
    }

}