#include <iostream>
#include <random>
#include <chrono>
#include <limits>
#include <mpi.h>
#include <vector>
#include <cstring>

using namespace std;

struct assets {
    double* pMatrix = nullptr;
    double* pVector = nullptr;
    double* pResult = nullptr;
    int Size = 0;
};

double GetRandom(double mod) {
    static std::mt19937 gen(
        static_cast<unsigned long>(
            std::chrono::steady_clock::now().time_since_epoch().count()
        )
    );
    if (mod <= 0) mod = 1.0;
    std::uniform_real_distribution<double> dist(-mod, mod);
    return dist(gen);
}

assets Generator(int size, double mod) {
    assets prod;
    prod.Size   = size;
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

void Free(assets& a) {
    delete[] a.pMatrix;  a.pMatrix = nullptr;
    delete[] a.pVector;  a.pVector = nullptr;
    delete[] a.pResult;  a.pResult = nullptr;
    a.Size = 0;
}

int main() {
    using clock = std::chrono::high_resolution_clock;

    int N = 0;
    const double mod = 10.0;

    int provided_argc = 0; char** provided_argv = nullptr;

    auto t0 = clock::now();
    MPI_Init(&provided_argc, &provided_argv);

    int rank = 0, procs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);

    if (rank == 0) {
        std::cout << "Enter matrix size N > 0: ";
        if (!(std::cin >> N) || N <= 0) {
            std::cerr << "Invalid N\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (N <= 10 || procs < 2 || procs > 8) {
        if (rank == 0) {
            std::cerr << "[MPI] Invalid parameters: N=" << N
                      << ", procs=" << procs << " (need N>10, 2..8 procs)\n";
        }
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    assets a;
    if (rank == 0) {
        a = Generator(N, mod);
    }

    std::vector<double> x(N);
    if (rank == 0) {
        std::memcpy(x.data(), a.pVector, sizeof(double) * N);
    }
    MPI_Bcast(x.data(), N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    std::vector<int> rows(procs, N / procs);
    for (int i = 0; i < N % procs; ++i) rows[i]++;

    std::vector<int> displ_rows(procs, 0);
    for (int i = 1; i < procs; ++i) displ_rows[i] = displ_rows[i - 1] + rows[i - 1];

    std::vector<int> sendcounts(procs), displs(procs);
    for (int i = 0; i < procs; ++i) {
        sendcounts[i] = rows[i] * N;
        displs[i]     = displ_rows[i] * N;
    }

    std::vector<double> A_local(sendcounts[rank]);
    std::vector<double> y_local(rows[rank], 0.0);

    MPI_Scatterv(
        (rank == 0 ? a.pMatrix : nullptr),
        sendcounts.data(), displs.data(), MPI_DOUBLE,
        A_local.data(), sendcounts[rank], MPI_DOUBLE,
        0, MPI_COMM_WORLD
    );

    for (int r = 0; r < rows[rank]; ++r) {
        const double* row = &A_local[r * N];
        double sum = 0.0;
        for (int j = 0; j < N; ++j) sum += row[j] * x[j];
        y_local[r] = sum;
    }

    std::vector<int> recvcounts(procs), recvdispl(procs);
    for (int i = 0; i < procs; ++i) {
        recvcounts[i] = rows[i];
        recvdispl[i]  = displ_rows[i];
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

    MPI_Barrier(MPI_COMM_WORLD);
    auto t1 = clock::now();
    double local_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    double global_ms = 0.0;
    MPI_Reduce(&local_ms, &global_ms, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    MPI_Finalize();

    if (rank == 0) {
        std::cout << "[MPI] Matrix " << N << "x" << N
                  << " | Processes: " << procs
                  << " | Total time (incl. init): " << global_ms << " ms\n";
        Free(a);
    }

    return 0;
}
