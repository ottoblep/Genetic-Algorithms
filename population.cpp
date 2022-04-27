#include <stdlib.h>
#include <vector>
#include <iostream>
#include <string>
#include <random>
#include <time.h>       
#include <chrono>
using namespace std;

class population {
    public:

    int popsize;
    int parametercount;
    vector<vector<double>> parameters;
    default_random_engine initialize_random(){
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator (seed);
        return generator;
    }
    default_random_engine rdev = initialize_random();

    population(int sizein,int parametercount): // Creates individuals with random parameters
        popsize(sizein), parametercount(parametercount){
        uniform_real_distribution<double> norm(-1.0,1.0);
        vector<double> tmp;
        for (int specimen=0; specimen<popsize;specimen++){
            tmp.clear();
            for (int parami=0; parami<parametercount; parami++) {
                double rnum = norm(rdev);
                tmp.push_back(rnum);
            }
            parameters.push_back(tmp);
        }
    }

    void generation_step(int n){
        for (int i=0; i<n; i++){
            vector<int> selected = selection();
            if (selected.size()==0 | selected.size()>95){ 
                cout << "Abort due to selection fault."; 
                return;
            }
            procreate(selected); 
        }
    }

    private:
    
    vector<int> selection(){ // Returns indices of individuals with better than average performance
       double b = 3;
       vector<int> selectedlist;
       double cumerrors = 0;
       double error;
       for (int specimen=0;specimen<popsize;specimen++){
            error = abs(3-parameters.at(specimen).at(0));
            cumerrors+=error; 
            //cout << "Error:" << error << " \n";
       } 
       double avgerror = cumerrors/popsize;
       for (int specimen=0;specimen<popsize;specimen++){
        
           error = abs (3-parameters.at(specimen).at(0));
           if (error<avgerror) {
                selectedlist.push_back(specimen); 
                //cout << "Selected " << specimen << "\n";
           }
       }
       cout<<"AvgError: "<<avgerror<<"\n";
       return selectedlist;
    }

    void procreate(vector<int> selected){ // Moves to a new generation
    
        vector<vector<double>> newparameters;
        for (int specimen=0; specimen<popsize; specimen++){
            vector<double> newgene = create_genes(selected);
            newparameters.push_back(newgene);
        }
        parameters = newparameters;
    }

    vector<double> create_genes(vector<int> selected){ // Creates a new individual from selected older individuals
    
        int numselected = selected.size();
        uniform_int_distribution<int> uni(0,numselected-1); // guaranteed unbiased
        normal_distribution<double> norm(0,0.1);

        //cout << selected.size() << " ";
        vector<double> mother = parameters.at(selected.at(uni(rdev) % numselected));
        vector<double> father = parameters.at(selected.at(uni(rdev) % numselected));
        vector<double> newgene;
        newgene = mother; // Clone mother
        for (int param=0;param<parametercount;param++){ 
            newgene.at(param) = (newgene.at(param) + father.at(param))/2; // Add father genome
            double rnum = norm(rdev);
            newgene.at(param) -= rnum; // Add random mutation;
        }
        //cout << "Gen " << newgene.at(0) << " \n";
        return newgene;
    }
};

int main(){
    // Solve a = 3
    population pop(100,1);
    pop.generation_step(100);
    return 0;
}