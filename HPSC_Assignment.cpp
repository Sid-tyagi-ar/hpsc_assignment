#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <omp.h>
#include <cstdlib>

void zero_forces(int N, std::vector<std::array<double,3>> &F) {
    for (int i = 0; i < N; i++)
        F[i] = {0.0, 0.0, 0.0};
}

void add_gravity(int N,
                 std::vector<std::array<double,3>> &F,
                 const std::vector<double> &m,
                 double g) {
    for (int i = 0; i < N; i++)
        F[i][2] -= m[i] * g;
}

void compute_particle_contacts(
    int N,
    const std::vector<std::array<double,3>> &x,
    const std::vector<std::array<double,3>> &v,
    std::vector<std::array<double,3>> &F,
    const std::vector<double> &R,
    double kn,
    double gamma_n
) {
    int num_threads = omp_get_max_threads();

    std::vector<std::vector<std::array<double,3>>> F_local(
        num_threads,
        std::vector<std::array<double,3>>(N, {0.0,0.0,0.0})
    );

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < N; i++) {
            for (int j = i+1; j < N; j++) {

                double rij[3] = {
                    x[j][0] - x[i][0],
                    x[j][1] - x[i][1],
                    x[j][2] - x[i][2]
                };

                double dij = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
                double delta = R[i] + R[j] - dij;

                if (delta > 0 && dij > 1e-12) {

                    double nij[3] = {
                        rij[0]/dij,
                        rij[1]/dij,
                        rij[2]/dij
                    };

                    double vij[3] = {
                        v[j][0] - v[i][0],
                        v[j][1] - v[i][1],
                        v[j][2] - v[i][2]
                    };

                    double vn = vij[0]*nij[0] + vij[1]*nij[1] + vij[2]*nij[2];

                    double Fn = kn*delta - gamma_n*vn;
                    if (Fn < 0) Fn = 0;

                    for(int k=0;k<3;k++){
                        double val = Fn * nij[k];
                        F_local[tid][i][k] += val;
                        F_local[tid][j][k] -= val;
                    }
                }
            }
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        for (int t = 0; t < num_threads; t++) {
            for(int k=0;k<3;k++)
                F[i][k] += F_local[t][i][k];
        }
    }
}

void compute_wall_contacts(int N,
    const std::vector<std::array<double,3>> &x,
    const std::vector<std::array<double,3>> &v,
    std::vector<std::array<double,3>> &F,
    const std::vector<double> &R,
    double kn,
    double gamma_n) {

    for(int i=0;i<N;i++){
        double delta = R[i] - x[i][2];
        if(delta>0){
            double vn = -v[i][2];
            double Fn = kn*delta - gamma_n*vn;
            if(Fn<0) Fn=0;
            F[i][2] += Fn;
        }
    }
}

void integrate(int N,
    std::vector<std::array<double,3>> &x,
    std::vector<std::array<double,3>> &v,
    const std::vector<std::array<double,3>> &F,
    const std::vector<double> &m,
    double dt){

    for(int i=0;i<N;i++){
        for(int k=0;k<3;k++){
            v[i][k] += (F[i][k]/m[i])*dt;
            x[i][k] += v[i][k]*dt;
        }
    }
}

int main(){

    double g = 9.81;

    {
        std::ofstream file("freefall.txt");

        int N = 1;
        double dt = 1e-4;
        int steps = 2000;

        std::vector<std::array<double,3>> x(N), v(N), F(N);
        std::vector<double> m(N,1.0), R(N,0.01);

        x[0] = {0,0,1.0};
        v[0] = {0,0,0};

        for(int step=0;step<steps;step++){
            double t = step*dt;

            zero_forces(N,F);
            add_gravity(N,F,m,g);
            integrate(N,x,v,F,m,dt);

            file << t << " " << x[0][2] << "\n";
        }
    }

    {
        std::ofstream errfile("error_dt.txt");

        for(double dt_test : {1e-3,5e-4,1e-4,5e-5}){

            int N = 1;
            std::vector<std::array<double,3>> x(N), v(N), F(N);
            std::vector<double> m(N,1.0), R(N,0.01);

            x[0] = {0,0,1.0};
            v[0] = {0,0,0};

            int steps = 1.0/dt_test;

            for(int step=0;step<steps;step++){
                zero_forces(N,F);
                add_gravity(N,F,m,g);
                integrate(N,x,v,F,m,dt_test);
            }

            double t = steps*dt_test;
            double z_exact = 1.0 - 0.5*g*t*t;
            double error = fabs(x[0][2] - z_exact);

            errfile << dt_test << " " << error << "\n";
        }
    }

    int N = 3000;
    double dt = 1e-4;
    int steps = 2000;

    std::vector<std::array<double,3>> x(N), v(N), F(N);
    std::vector<double> m(N,1.0), R(N,0.01);

    for(int i=0;i<N;i++){
        x[i] = {(double)rand()/RAND_MAX, (double)rand()/RAND_MAX, (double)rand()/RAND_MAX+0.5};
        v[i] = {0,0,0};
    }

    std::ofstream bounce("bounce.txt");
    std::ofstream energy("energy.txt");

    double start = omp_get_wtime();

    for(int step=0;step<steps;step++){

        zero_forces(N,F);
        add_gravity(N,F,m,g);

        compute_particle_contacts(N,x,v,F,R,1e5,10);
        compute_wall_contacts(N,x,v,F,R,1e5,10);

        integrate(N,x,v,F,m,dt);

        bounce << step*dt << " " << x[0][2] << "\n";

        double KE = 0;
        for(int i=0;i<N;i++){
            double v2 = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
            KE += 0.5*m[i]*v2;
        }
        energy << step*dt << " " << KE << "\n";

        if(step % 500 == 0){
            std::string filename = "snap_" + std::to_string(step) + ".txt";
            std::ofstream snap(filename);
            for(int i=0;i<N;i++){
                snap << x[i][0] << " " << x[i][2] << "\n";
            }
        }
    }

    double end = omp_get_wtime();
    std::cout << "Time: " << end-start << std::endl;

    return 0;
}
