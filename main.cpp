#include <iostream>
#include <random>
#include <chrono>
#include <limits>
using namespace std;

struct assets {
    double* pMatrix = nullptr;
    double* pVector = nullptr;
    double* pResult = nullptr;
    int Size = 0;
};

double GetRandom(double mod) {
    static mt19937 gen(chrono::steady_clock::now().time_since_epoch().count());
    if (mod <= 0) mod = 1.0;
    uniform_real_distribution<double> dist(-mod, mod);
    return dist(gen);
}

double D_Validator() {
    double x;
    for (;;) {
        if (cin >> x) return x;
        cout << "WRONG, try again\n";
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
    }
}

int I_Validator() {
    int x;
    for (;;) {
        if (cin >> x) return x;
        cout << "WRONG, try again\n";
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
    }
}

assets Generator(double mod) {
    assets prod;
    int size = 0;
    for (;;) {
        cout << "enter size\n";
        size = I_Validator();
        if (size > 0) break;
        cout << "wrong, try again\n";
    }
    prod.Size = size;
    prod.pMatrix = new double[size * size];
    prod.pVector = new double[size];
    prod.pResult = new double[size];

    for (int i = 0; i < size; ++i) {
        prod.pVector[i] = GetRandom(mod);
        prod.pResult[i] = 0.0;
    }
    for (int i = 0; i < size * size; ++i)
        prod.pMatrix[i] = GetRandom(mod);

    return prod;
}

assets Manual() {
    assets prod;
    int size = 0;
    for (;;) {
        cout << "enter size\n";
        size = I_Validator();
        if (size > 0) break;
        cout << "wrong, try again\n";
    }
    prod.Size = size;
    prod.pMatrix = new double[size * size];
    prod.pVector = new double[size];
    prod.pResult = new double[size];

    cout << "Enter first row of matrix\n";
    for (int i = 0; i < size; ++i) {
        cout << "Enter element " << i << "\n";
        prod.pMatrix[i] = D_Validator();
    }
    for (int i = size; i < size * size; ++i)
        prod.pMatrix[i] = prod.pMatrix[i % size] + 1.0 * (i / size);

    cout << "Enter vector\n";
    for (int i = 0; i < size; ++i) {
        cout << "Enter element " << i << "\n";
        prod.pVector[i] = D_Validator();
        prod.pResult[i] = 0.0;
    }
    return prod;
}

void calc(assets& a) {
    auto start = chrono::high_resolution_clock::now();
    for (int i = 0; i < a.Size; ++i)
        for (int j = 0; j < a.Size; ++j)
            a.pResult[i] += a.pMatrix[i * a.Size + j] * a.pVector[j];
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> dt = end - start;
    cout << "Calculation completed in " << dt.count() << " ms\n";
}

void Free(assets& a) {
    delete[] a.pMatrix;  a.pMatrix = nullptr;
    delete[] a.pVector;  a.pVector = nullptr;
    delete[] a.pResult;  a.pResult = nullptr;
    a.Size = 0;
}

int main() {
    cout << "Serial matrix-vector multiplication program\n";
    double mod = 10.0;
    int state = 1;

    while (state) {
        switch (state) {
            case 1: {
                cout << "choose option\n"
                     << "//////////////\n"
                     << "2. input matrix-vector manually and calculate\n"
                     << "3. generate matrix-vector randomly and calculate\n"
                     << "0. exit\n";
                state = I_Validator();
                break;
            }
            case 2: {
                assets vars = Manual();
                calc(vars);
                Free(vars);
                state = 1;
                break;
            }
            case 3: {
                assets vars = Generator(mod);
                calc(vars);
                Free(vars);
                state = 1;
                break;
            }
            case 0: {
                return 0;
            }
            default: {
                cout << "wrong option, try again\n";
                state = 1;
                break;
            }
        }
    }
    return 0;
}
