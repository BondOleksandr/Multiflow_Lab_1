#include <iostream>
#include <random>
#include <chrono>
#include <limits>
#include <mpi.h>
#include <vector>
#include <cstring>
#include <algorithm>

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

double D_Validator() {//makes sure your input is double
    double x;
    for (;;) {
        if (cin >> x) return x;
        cout << "WRONG, try again\n";
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
    }
}

int I_Validator() {//makes sure your input is int
    int x;
    for (;;) {
        if (cin >> x) return x;
        cout << "WRONG, try again\n";
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
    }
}

assets Generator(double mod) {//generates matrix and vector, asks for size
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

assets Manual() {//manual input of matrix and vector
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

void calc(assets& a) {//single-thread calculation, measures time
    auto start = chrono::high_resolution_clock::now();
    for (int i = 0; i < a.Size; ++i)
        for (int j = 0; j < a.Size; ++j)
            a.pResult[i] += a.pMatrix[i * a.Size + j] * a.pVector[j];
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> dt = end - start;
    cout << "Calculation completed in " << dt.count() << " ms\n";
}

void multiflow_calc(assets &a) {//MPI-CALCULATION
    // Start measuring total time
    auto start = std::chrono::high_resolution_clock::now();
    int argc = 0; char** argv = nullptr;
    MPI_Init(&argc, &argv);

    int rank = 0, procs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    int N = (rank == 0 ? a.Size : 0);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (N <= 10 || procs < 2 || procs > 8) {
        if (rank == 0)
            std::cerr << "[MPI] Invalid parameters (N=" << N
                      << ", procs=" << procs << ")\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Broadcast vector x
    std::vector<double> x(N);
    if (rank == 0) std::memcpy(x.data(), a.pVector, sizeof(double) * N);
    MPI_Bcast(x.data(), N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Distribute rows of matrix among processes
    std::vector<int> rows(procs, N / procs);
    for (int i = 0; i < N % procs; ++i) rows[i]++;
    std::vector<int> displ_rows(procs, 0);
    for (int i = 1; i < procs; ++i)
        displ_rows[i] = displ_rows[i - 1] + rows[i - 1];

    std::vector<int> sendcounts(procs), displs(procs);
    for (int i = 0; i < procs; ++i) {
        sendcounts[i] = rows[i] * N;
        displs[i] = displ_rows[i] * N;
    }

    std::vector<double> A_local(sendcounts[rank]);
    std::vector<double> y_local(rows[rank], 0.0);

    MPI_Scatterv(
        (rank == 0 ? a.pMatrix : nullptr), sendcounts.data(), displs.data(), MPI_DOUBLE,
        A_local.data(), sendcounts[rank], MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );

    // Local computation
    for (int r = 0; r < rows[rank]; ++r) {
        const double* row = &A_local[r * N];
        double sum = 0.0;
        for (int j = 0; j < N; ++j)
            sum += row[j] * x[j];
        y_local[r] = sum;
    }

    // Gather results
    std::vector<int> recvcounts(procs), recvdispl(procs);
    for (int i = 0; i < procs; ++i) {
        recvcounts[i] = rows[i];
        recvdispl[i] = displ_rows[i];
    }

    if (rank == 0) {
        if (!a.pResult) a.pResult = new double[N];
        a.Size = N;
    }

    MPI_Gatherv(
        y_local.data(), rows[rank], MPI_DOUBLE,
        (rank == 0 ? a.pResult : nullptr),
        recvcounts.data(), recvdispl.data(), MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );

    if (rank == 0)
        std::cout << "\n[MPI] Matrix " << N << "x" << N
                  << " | Processes: " << procs
                  << " | Total time: " << max_ms << " ms\n";

    MPI_Finalize();

    auto end = std::chrono::high_resolution_clock::now();
    chrono::duration<double, milli> dt = end - start;
    cout << "Calculation completed in " << dt.count() << " ms\n";
}

void Free(assets& a) {//cleaner
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
