#include <iostream>
#include <cmath>
#include <random>
#include <chrono>

using namespace std;

struct assets{ // 
double* pMatrix; // First argument - initial matrix
double* pVector; // Second argument - initial vector
double* pResult; // Result vector for matrix-vector multiplication
int Size; // Sizes of initial matrix and vecto
};

// Generates random number modulary not greater
double GetRandom(double mod) {
    static std::mt19937 generator(
        std::chrono::steady_clock::now().time_since_epoch().count()
    );

    if(mod<=0){
        if(mod==0)mod++;
        else mod = mod*(-1);
    }

    std::uniform_real_distribution<double> distribution(-mod, mod);
    return distribution(generator);
}

assets Generator(int size, double mod){
    double mod;
    assets prod;
    prod.Size = size;

    if(size<=0){
        cout << "wron input";
        return prod;
    }

    prod.pMatrix = new double[size*size];
    prod.pVector = new double[size];
    prod.pResult = new double[size];

    for(int i=0; i<size; i++){
        prod.pVector[i]=GetRandom(mod);
        prod.pResult[i]=0;
    }
    
    for(int i=0; i<size*size; i++){
        prod.pMatrix[i] = GetRandom(mod);
    }

    return prod;
}

double D_Validator(){//makes sure your input is double
    double input;
    bool aux = true;
    while(aux){
    cin >> input;
    if(cin.fail()){
        cout<<"WRONG,try again"<<endl;
        cin.clear();
        cin.ignore(10000, '\n');
    }
    else aux = false;
}
return input;
}

assets Manual(){//manual input
    int size;
    bool aux=true;
    assets prod;

    while(aux){
        cout << "input size: " << endl;
        cin >> size;
        if(size<=0)cout << "wrong, try again " << endl;
        else aux = 0;
    }

    prod.Size = size;
    prod.pMatrix = new double[size*size];
    prod.pVector = new double[size];
    prod.pResult = new double[size];

    cout << "Enter first row of matrix" << endl;

    for(int i=0; i<size;i++){//first row
            cout <<"Enter element " << i <<endl;
            prod.pMatrix[i]=D_Validator();
    }
    
    for(int i=size;i<size*size;i++){//fills the rest
        prod.pMatrix[i]=prod.pMatrix[i%size]+1*(i/size);
    }

    cout << "Enter vector " << endl;

    for(int i=0;i<size;i++){
        cout <<"Enter element " << i <<endl;
        prod.pVector[i]=D_Validator();
        prod.pResult[i] = 0;
    }

    return prod;
}

int main() {
    cout << "Serial matrix-vector multiplication program" << endl;
    return 0;
}