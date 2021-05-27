#include <string>
#include <fstream>
#include <vector>
#include <utility> // std::pair

void write_csv(std::string filename, std::vector<std::pair<std::string, std::vector<float>>> transfer_func){
    // we meed both a string and an integer for parameters.
    std::ofstream myFile(filename);

    // Send column names to the stream
    for(int j = 0; j < transfer_func.size(); ++j)
    {
        myFile << transfer_func.at(j).first;// the string is for the columns
        if(j != transfer_func.size() - 1) myFile << ","; // No comma at end of line
    }
    myFile << "\n";

    // Send data to the stream
    for(int i = 0; i < transfer_func.at(0).second.size(); ++i)
    {
        for(int j = 0; j < transfer_func.size(); ++j)
        {
            myFile << transfer_func.at(j).second.at(i);
            if(j != transfer_func.size() - 1) myFile << ","; // No comma at end of line
        }
        myFile << "\n";
    }
    // Close the file
    myFile.close();
}

int main() {
    // Make three vectors, each of length 100 filled with 1s, 2s, and 3s
    std::vector<float> vec1;
    vec1.push_back(50);
    vec1.push_back(100);
    vec1.push_back(500);
    vec1.push_back(1000);
    vec1.push_back(5000);
    vec1.push_back(10000);
    vec1.push_back(50000);

   // std::vector<int> vec2(100, 2);
    std::vector<float> vec3;
    vec3.push_back(-12.8);vec3.push_back(-8.62);
    vec3.push_back(-2.62);vec3.push_back(-1.43);
    vec3.push_back(-0.323);vec3.push_back(-0.188);
    vec3.push_back(-0.100);


    // Wrap into a vector
    std::vector<std::pair<std::string, std::vector<float>>> vals = {{"Frequency", vec1}, {"Gain", vec3}};

    // Write the vector to CSV
    write_csv("transfunc.csv", vals);

    return 0;
}
