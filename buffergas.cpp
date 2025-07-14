#include "buffergas.h"

BufferGas::BufferGas()
{
    time_t timer;
    tm startTime;
    double seconds;

    mass = 6.64648E-27;             // mass for He, unit kg
    temperature = 4.5;              // unit K
    density = 1.50456E29 * 1.024e-6; // unit particle/m^3
    vbThermSquare = 8 * kB * 4.5 / (M_PI * mass);

    // get seconds since 01.01.2000
    startTime.tm_year = 2000;
    startTime.tm_mon = 0;
    startTime.tm_mday = 1;
    startTime.tm_hour = 0;
    startTime.tm_min = 0;
    startTime.tm_sec = 0;

    time(&timer);
    seconds = difftime(timer, mktime(&startTime));

    // initialize random number generator with time
    generator.seed(seconds);
}

bool BufferGas::is_in_volume(double const coord[3])
{
    double phi;
    double r[3];
    //is in sphere or cylinder
    if ((coord[1] >= -height && coord[1] < 0 && coord[0] * coord[0] + coord[2] * coord[2] <= halfDiameterSquare) ||
        coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2] <= halfDiameterSquare)
        return true;
    //is in exit region: first cylinder
    if (coord[2] > 0 && coord[2] <= halfDiameter + 0.0015
        && coord[0] * coord[0] + coord[1] * coord[1] <= pow(0.003, 2))
        return true;
    //is in exit region: conus
    if (coord[2] > halfDiameter + 0.0015 && coord[2] <= halfDiameter + 0.0015 + 0.003 - exitDiameter * 0.5
        && pow(0.003 - (coord[2] - halfDiameter - 0.0015), 2) >= coord[0] * coord[0] + coord[1] * coord[1])
        return true;
    //is in exit region: nozzle
    if (coord[2] > 0 && coord[2] <= halfDiameter + 0.0015 + 0.003 - exitDiameter * 0.5 + 0.0005
        && pow(exitDiameter * 0.5, 2) >= coord[0] * coord[0] + coord[1] * coord[1])
        return true;
    //entry region and other laser window
    phi = -50 * M_PI / 180;
    r[0] = cos(phi) * coord[0] + sin(phi) * coord[2];
    r[1] = coord[1];
    r[2] = -sin(phi) * coord[0] + cos(phi) * coord[2];
    if (abs(r[2]) <= halfDiameter + 0.0015 && r[0] * r[0] + r[1] * r[1] <= pow(0.002, 2))
        return true;
    return false;
}

bool BufferGas::hits_exit(double const coord[3]) {

    //if (coord[2] > halfDiameter+0.0015+0.003-exitDiameter*0.5+0.0005 && pow(exitDiameter*0.5,2) >= coord[0]*coord[0] + coord[1]*coord[1] )
    if (coord[2] > halfDiameter + 0.0005)
        return true;
    else
        return false;
}

void BufferGas::load_Data(const std::string filename)
{
    int ix, iy, iz;
    char c;
    std::fstream file;

    file.open(filename, std::ios::in);

    file >> c;//slot at the beginning of the header
    file >> diameter >> height >> exitDiameter >> delta >> nx >> ny >> nz;
    halfDiameter = diameter * 0.5;
    halfDiameterSquare = halfDiameter * halfDiameter;

    fField.resize(nx * ny * nz);

    for (ix = 0; ix < nx; ++ix)
    {
        for (iy = 0; iy < ny; ++iy)
        {
            for (iz = 0; iz < nz; ++iz)
            {
                file >> fField[ix + nx * (iy + ny * iz)];
                if (!file.good())
                {
                    std::cerr << "Error during load velocities from " << filename << "\n";
                    return;
                }
            }
        }
    }
    file.close();
}

double BufferGas::get_radius()
{
    return halfDiameter;
}

void BufferGas::update_coord(const double* coord)
{
    int ix, iy, iz;
    if (nx == 0)
    {
        std::cerr << "Error: No buffer gas data loaded!\n";
    }
    else
    {
        ix = static_cast<int>((coord[0] + halfDiameter + 0.001) / delta);
        iy = static_cast<int>((coord[1] + height) / delta);
        iz = static_cast<int>((coord[2] + halfDiameter) / delta);
        if (ix != ixold || iy != iyold || iz != izold)
        {
            if (ix >= 0 and ix < nx and iy >= 0 and iy < ny and iz >= 0 and iz < nz)
            {
                vDrift[0] = fField[ix + nx * (iy + ny * iz)].v[0];
                vDrift[1] = fField[ix + nx * (iy + ny * iz)].v[1];
                vDrift[2] = fField[ix + nx * (iy + ny * iz)].v[2];
                density = fField[ix + nx * (iy + ny * iz)].rho;
                temperature = fField[ix + nx * (iy + ny * iz)].T;
                if (density == 0)
                {
                    density = 1.50456E29 * 1.024e-6; // unit particle/m^3
                }
                if (temperature == 0)
                {
                    temperature = 4.5; //unit K
                }
                vbThermSquare = 8 * kB * temperature / (M_PI * mass);
                ixold = ix;
                iyold = iy;
                izold = iz;
            }
            else
            {
                vDrift[0] = 0;
                vDrift[1] = 0;
                vDrift[2] = 0;
                density = 1.50456E29 * 1.024e-6; // unit particle/m^3
                temperature = 4.5; // K
                vbThermSquare = 8 * kB * 4.5 / (M_PI * mass);
            }
        }
    }
}

void BufferGas::get_thermal_velocity(double* u)
{
    u[0] = sqrt(kB * temperature / mass) * normalDistribution(generator);
    u[1] = sqrt(kB * temperature / mass) * normalDistribution(generator);
    u[2] = sqrt(kB * temperature / mass) * normalDistribution(generator);
}
