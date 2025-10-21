#include "molecule.h"

Molecule::Molecule()
{
    time_t timer;
    tm startTime;
    double seconds;

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

    // set mass of naphthalene
    mass = 2.12881E-25; // unit kg
}

void Molecule::init_coordinates_velocity(const double temperature, double radius)
{
    //double radius = 0.00949; // unit m
    double theta, phi;      // angles of spheric cordinates
    double v, dummy;
    radius += 0.001;

    // angles in Deg, standard deviation of angles is 2 Deg
    theta = 130 + 4. * normalDistribution(generator);
    phi = 180 + 4. * normalDistribution(generator);
    // convert from Deg to Rad
    phi *= M_PI / 180;
    theta *= M_PI / 180;

    coord[0] = radius * sin(theta) * cos(phi);
    coord[1] = radius * sin(theta) * sin(phi);
    coord[2] = radius * cos(theta);

    // determine norm of velocity with help of normal distribution
    v = pow(normalDistribution(generator), 2.);
    v += pow(normalDistribution(generator), 2.);
    v += pow(normalDistribution(generator), 2.);
    // consider sigma = sqrt(kB*temperature/mass), here v^2
    v *= kB * temperature / mass;

    v = sqrt(v);
    // determine direction of v, angles now in Rad
    //theta only between 0 Deg and 90 Deg
    theta = inverse_cos_square_dist(uniformDistribution(generator));
    phi = 2 * M_PI * uniformDistribution(generator);

    velocity[0] = v * sin(theta) * cos(phi);
    velocity[1] = v * sin(theta) * sin(phi);
    velocity[2] = v * cos(theta);

    // rotate around y axis at 50 Deg
    //z nach vorne (von mir weg)
    //y nach oben
    //x nach links
    phi = 50 * M_PI / 180;
    dummy = velocity[0];
    velocity[0] = cos(phi) * dummy + sin(phi) * velocity[2];
    velocity[2] = -sin(phi) * dummy + cos(phi) * velocity[2];
}

double Molecule::inverse_cos_square_dist(const double u)
{
    // calculate the inverse of the distribution function
    //     F(x) = ( 2*x + sin(2*x) )/Pi
    //   of the distribution
    //     f(x) = 2/Pi * cos(x)^2
    //   with x in [0, Pi/2]
    //   using the regula falsi method, which searches the zero points of the function
    //     g(x) = ( 2*x + sin(2*x) )/Pi - u
    double ak = 0, bk = M_PI_2, ck;
    double fak = -u, fbk = 1. - u, fck;

    if (u == 0)
    {
        return ak;
    }
    else if (u == 1)
    {
        return bk;
    }
    else if (u < 0 or u > 1)
    {
        std::cerr << "inverse cos square distribtion of " << u << " is not defined!";
        return -10;
    }

    while (true)
    {
        ck = ak - (bk - ak) / (fbk - fak) * fak;
        fck = (2 * ck + sin(2 * ck)) / M_PI - u;
        if (abs(fck) < 1.0E-12)
        {
            break;
        }
        if (fck > 0)
        {
            bk = ck;
            fbk = fck;
        }
        else
        {
            ak = ck;
            fak = fck;
        }
    }

    return ck;
}

void Molecule::check_cos_square_dist()
{
    double d;
    int array[1000]{ 0 };
    int i;

    for (i = 0; i < 10000000; i++)
    {
        d = inverse_cos_square_dist(uniformDistribution(generator)) / M_PI_2;
        ++array[int(d * 1000.)];
    }
    std::ofstream file;
    file.open("cos_square_distribution.dat", std::ios::out);
    for (i = 0; i < 1000; i++)
    {
        file << i * M_PI_2 / 1000 << " " << array[i] << "\n";
    }
    file.close();
}

double Molecule::get_norm_vrel(double vb[3])
{
    return sqrt(pow(velocity[0] - vb[0], 2) + pow(velocity[1] - vb[1], 2) + pow(velocity[2] - vb[2], 2));
}

double Molecule::get_temperature(double vbDrift[3])
{
    return (pow(velocity[0] - vbDrift[0], 2) + pow(velocity[1] - vbDrift[1], 2) + pow(velocity[2] - vbDrift[2], 2)) * mass / (3 * kB);
}

double Molecule::get_uniform_distributed_variable()
{
    return uniformDistribution(generator);
}

void Molecule::set_velocity_from_hard_sphere_collision_with_buffer_gas(double vb[3], double bufferGasMass)
{
    double vcm[3]; // velocity of center of mass
    double g[3];
    double vrel; // relative velocity between buffer gas and molecule
    double totalMass = bufferGasMass + mass;
    int i;
    double r1 = uniformDistribution(generator);
    double r2 = uniformDistribution(generator);
    double phi, sintheta, costheta; // angles of spherical coordinates

    vrel = sqrt(pow(velocity[0] - vb[0], 2) + pow(velocity[1] - vb[1], 2) + pow(velocity[2] - vb[2], 2));

    for (i = 0; i < 3; i++)
    {
        vcm[i] = mass / totalMass * velocity[i] + bufferGasMass / totalMass * vb[i];
    }
    costheta = 2 * r1 - 1;
    sintheta = sqrt(1 - costheta * costheta);
    phi = 2 * M_PI * r2;

    g[0] = vrel * sintheta * cos(phi);
    g[1] = vrel * sintheta * sin(phi);
    g[2] = vrel * costheta;

    // update velocity of molecule after collision with buffer gas
    for (i = 0; i < 3; i++)
    {
        velocity[i] = vcm[i] + bufferGasMass / totalMass * g[i];
    }
}
