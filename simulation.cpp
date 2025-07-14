#include <stdio.h>
#include <math.h>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <random>
#include <time.h>

#include "molecule.h"
#include "buffergas.h"
#include "flowfieldentries.h"


int main(int argc, char **argv)
{
    int run, i, numCollisions, inSphere;
    double sigmaBS = 1.2E-17; // unit 1/m^3 for naphthalene
    double collisionProbability = 0.1;
    double collisionRate, distance;
    double dt, t;
    std::fstream file, fileInfo, fileSteps;
    std::string fileName, intStr;
    unsigned int step = 0;
    BufferGas gasHe;
    Molecule naphthalene;
    clock_t cpuT;
    std::vector<double> v1;
    const int numPartners=10;
    double vbArray[numPartners][3];
    double relvelocities[numPartners];
    double vtmp[3];
    double norm,randomNumber;
    std::string name, outputname;
    int pos;
    double initcoord[3];

    cpuT = clock();

    name = argv[1];
    //name = "voxelized_D16mm_H8mm_d1.5mm_0deg_15Pa_0mW.txt";
    pos = name.find_last_of("/");
    //remove path and "voxelized"
    outputname = name.substr(pos+1+9);

    gasHe.load_Data(name);

    fileInfo.open("run_infos"+outputname, std::ios::out);
    fileInfo << "#runNr inSphere timesteps time numCollisions endx endy endz\n";

    run = 0;
    inSphere = 1;
    while (run < 100000){
        intStr = std::to_string(run);
        fileName = "naphthalene_steps_" + intStr + ".txt";
        fileSteps.open(fileName, std::ios::out);
        fileSteps << "#numCollisions t x y z vx vy vz\n";
        fileSteps << "#numCollisions ()  t (s)   x y z (m)    vx vy vz (m/s)\n";

        t = 0;
        numCollisions = 0;

        naphthalene.init_coordinates_velocity(500., gasHe.get_radius());
        initcoord[0] =naphthalene.coord[0];
        initcoord[1] =naphthalene.coord[1];
        initcoord[2] =naphthalene.coord[2];

        gasHe.update_coord(naphthalene.coord);

        fileSteps << numCollisions << " " <<  t << " ";
        fileSteps << naphthalene.coord[0] << " "  << naphthalene.coord[1] << " "  << naphthalene.coord[2] << " "  ;
        fileSteps << naphthalene.velocity[0] << " " << naphthalene.velocity[1] << " " << naphthalene.velocity[2] << "\n";


        for (step = 1; step <= 1000000000; ++step)
        {
            collisionRate = sigmaBS * gasHe.density * sqrt(gasHe.vbThermSquare 
                                                        + pow(naphthalene.velocity[0] - gasHe.vDrift[0], 2) 
                                                        + pow(naphthalene.velocity[1] - gasHe.vDrift[1], 2) 
                                                        + pow(naphthalene.velocity[2] - gasHe.vDrift[2], 2));
            dt = collisionProbability / collisionRate;

            for (i=0; i<3; ++i)
            {
                naphthalene.coord[i] += naphthalene.velocity[i] * dt;
            }
            t += dt;

            if (naphthalene.get_uniform_distibuted_variable() <= collisionProbability)
            {   
                
                //generate several possible collision partners with different velocities
                for (i=0; i<numPartners; ++i){
                    gasHe.get_thermal_velocity(vbArray[i]);
                    vbArray[i][0] += gasHe.vDrift[0];
                    vbArray[i][1] += gasHe.vDrift[1];
                    vbArray[i][2] += gasHe.vDrift[2];
                    // calculate relative velocitiy with possible collision partner i
                    vtmp[0] = vbArray[i][0] - naphthalene.velocity[0];
                    vtmp[1] = vbArray[i][1] - naphthalene.velocity[1];
                    vtmp[2] = vbArray[i][2] - naphthalene.velocity[2];
                    relvelocities[i] = sqrt(vtmp[0]*vtmp[0] + vtmp[1]*vtmp[1] + vtmp[2]*vtmp[2]);
                }
                //get total probability
                norm = 0;
                for (i=0; i<numPartners; ++i){
                    norm += relvelocities[i];
                }
                //choose collision partner: particle with higher relative velocity have higher collsion probability
                randomNumber = naphthalene.get_uniform_distibuted_variable()*norm;
                norm = 0;
                for (i=0; i<numPartners; ++i){
                    norm += relvelocities[i];
                    if (norm >= randomNumber){
                        break;
                    }
                }
                naphthalene.set_velocity_from_hard_shpere_collision_with_buffer_gas(vbArray[i], gasHe.mass);
                gasHe.update_coord(naphthalene.coord);
                     
                
                ++numCollisions;

                if ( (numCollisions<=200 && numCollisions % 1==0) || (numCollisions % 1000==0)){
                    fileSteps << numCollisions << " ";
                    fileSteps << t << " ";
                    fileSteps << naphthalene.coord[0] << " "  << naphthalene.coord[1] << " "  << naphthalene.coord[2] << " "  ;
                    fileSteps << naphthalene.velocity[0] << " " << naphthalene.velocity[1] << " " << naphthalene.velocity[2] << "\n";
                }
            }

            if (!gasHe.is_in_volume(naphthalene.coord))
            {
                if (gasHe.hits_exit(naphthalene.coord))
                {
                    inSphere = -1;
                }
                else
                {
                    inSphere = 0;
                }
                break;
            }
        }

        fileSteps.close();

        //r = sqrt(naphthalene.coord[0] * naphthalene.coord[0] + naphthalene.coord[1] * naphthalene.coord[1] + naphthalene.coord[2] * naphthalene.coord[2]);
        //theta = acos(naphthalene.coord[2] / r) * 180 / M_PI;
        //phi = atan2(naphthalene.coord[1], naphthalene.coord[0]) * 180 / M_PI;
        //temperature = naphthalene.get_temperature(vbDrift);

        //std::cout << inSphere << " tsteps " << step << " collisions " << numCollisions << " theta " << theta << " phi " << phi << " T " << temperature << '\n';

        fileInfo << run << ' ' << inSphere << ' ' << step << ' '  << t << ' ' << numCollisions << ' ';
        fileInfo << naphthalene.coord[0] << ' ' << naphthalene.coord[1] << ' ' << naphthalene.coord[2] << '\n';
        fileInfo.flush();

        //condition for saving the trajectory
        distance = pow(naphthalene.coord[0] - initcoord[0],2) + pow(naphthalene.coord[1] - initcoord[1],2) + pow(naphthalene.coord[2] - initcoord[2],2);
        if (distance > pow(0.003,2) ){
            ++run;
        } else {
            std::cout << distance << " is to small!\n";
        }
    }
    fileInfo.close();
    fileInfo.clear();

    cpuT = clock() - cpuT;
    std::cout << "CPU time: " << static_cast<double>(cpuT) / CLOCKS_PER_SEC << "s \n";

    return 0;
}