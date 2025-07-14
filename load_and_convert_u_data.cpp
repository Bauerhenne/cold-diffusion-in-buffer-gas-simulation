#include "vec3d.h"
#include "flowfieldentries.h"
#include<math.h>
#include<fstream>
#include<string>
#include<vector>
#include<iostream>


bool is_in_volume(double const coord[3], double const diameter, double const height, double const exitDiameter)
{
    double phi;
    double r[3];
    double halfDiameter = diameter * 0.5;
    double halfDiameterSquare = halfDiameter*halfDiameter;
    //is in sphere or cylinder
    if ((coord[1] >= -height && coord[1] < 0 && coord[0]*coord[0] + coord[2]*coord[2] <= halfDiameterSquare) ||
        coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2] <= halfDiameterSquare)
        return true;
    //is in exit region: first cylinder
    if (coord[2] > 0 && coord[2] <= halfDiameter+0.0015 
        && coord[0]*coord[0] + coord[1]*coord[1] <= pow(0.003,2))
        return true;
    //is in exit region: conus
    if (coord[2] >= halfDiameter+0.0015 && coord[2] <= halfDiameter+0.0015+0.003-exitDiameter*0.5 
        && pow(0.003 - (coord[2] - halfDiameter-0.0015),2) >= coord[0]*coord[0] + coord[1]*coord[1] ) 
        return true;
    //is in exit region: nozzle
    if (coord[2] > 0 && coord[2] <= halfDiameter+0.0015+0.003-exitDiameter*0.5+0.0005 
        && pow(exitDiameter*0.5,2) >= coord[0]*coord[0] + coord[1]*coord[1] )
        return true;
    //entry region and other laser window
    phi = -50 * M_PI / 180;
    r[0] = cos(phi) * coord[0] + sin(phi) * coord[2];
    r[1] = coord[1];
    r[2] = -sin(phi) * coord[0] + cos(phi) * coord[2];
    if (abs(r[2])<=halfDiameter+0.0015 && r[0]*r[0] + r[1]*r[1] <= pow(0.002,2))
        return true;
    return false;
}


int main(int argc, char **argv){
    std::ifstream filein;
    std::ofstream fileout;
    std::string filename;
    unsigned int i;
    std::vector<FlowFieldEntries> flowfielddata;
    std::vector<Vec3d> flowfieldcoord;
    Vec3d coorddummy;
    FlowFieldEntries flowfielddummy;
    std::vector<FlowFieldEntries> flowfield;
    std::vector<int> fieldcount;
    double dummy;
    double coord[3];
    double height;
    double diameter;
    double halfDiameter, exitDiameter;
    int ix,iy,iz, nx,ny,nz;
    int count;
    int countf=0;
    std::string name,s;
    unsigned int pos1, pos2;
    double delta;
    bool emptyvoxel=true;
    int tcount=0;

    name = argv[1];
    //name = "D16mm_H8mm_d1.5mm_0deg_15Pa_0mW.txt";

    pos1 = name.find(".txt");
    delta = 0.0005;
    
    //remove from name ".txt" and prepend "voxelixed_" and attach ".txt"
    // one may attach something else than .txt   
    filename = "voxelized_" + name.substr(0,name.length()-4) + ".txt";
    
    //extract the numbers from the filename
    pos1 = name.find("D");
    pos2 = name.find("mm");
    s = name.substr(pos1+1,pos2-pos1-1);
    diameter = std::stod(s);
    
    pos1 = name.find("H");
    pos2 = name.find("mm",pos1);
    s = name.substr(pos1+1,pos2-pos1-1);
    height = std::stod(s);

    pos1 = name.find("d");
    pos2 = name.find("mm",pos1);
    s = name.substr(pos1+1,pos2-pos1-1);
    exitDiameter = std::stod(s);

    diameter *= 0.001;//convert from mm to m
    exitDiameter *= 0.001; //convert from mm to m
    height *= 0.001; //convert from mm to m
    halfDiameter = diameter * 0.5; 
    
    nx = static_cast<int>((diameter+0.002) / delta) + 1;
    ny = static_cast<int>((halfDiameter + height) / delta) + 1;
    nz = static_cast<int>((diameter+0.005-exitDiameter*0.5) / delta) + 1;

    flowfield.resize(nx*ny*nz);
    fieldcount.resize(nx*ny*nz);

    std::fill(fieldcount.begin(), fieldcount.end(), 0);

    filein.open(name, std::ios::in);
    //do not take header into account
    for (i=0; i<9; ++i)
    {
        getline(filein,s);
    }
    while(filein.good())
    {        
        filein >> coorddummy;
        filein >> flowfielddummy;
        flowfielddummy.rho *= (6.022140e23/4.0026)*1e6; // convert g/cm^3 to parricles per m^3
        for (i=0; i<6; ++i)
        {
            filein >> dummy;
        }
        if (filein.good())
        {
            flowfielddata.emplace_back(flowfielddummy);
            flowfieldcoord.emplace_back(coorddummy);

        }
    }
    filein.close();

    std::cout << "Numdata: " << flowfielddata.size() << std::endl;

    if (flowfielddata.size() ==0 )
    {
        std::cout << "Error: No data loaded!\n";
        return 1;
    } 

    std::cout.precision(16);
    std::cout << "diameter:     " << diameter << " m\n";
    std::cout << "exitdiameter: " << exitDiameter << " m\n";
    std::cout << "height:       " << height << " m\n";
    std::cout << "halfdiameter: " << halfDiameter << " m\n";
    std::cout << "delta:        " << delta << " m\n";
    std::cout << "nx " << nx << "\n";
    std::cout << "ny " << ny << "\n";
    std::cout << "nz " << nz << "\n";

    for(i=0; i<flowfielddata.size(); ++i)
    {
        ix = static_cast<int>(round(( flowfieldcoord[i].r[0] + halfDiameter+0.001 ) / delta));
        iy = static_cast<int>(round(( flowfieldcoord[i].r[1] + height ) / delta));
        iz = static_cast<int>(round(( flowfieldcoord[i].r[2] + halfDiameter ) / delta));

        if (ix >=0 && ix < nx && iy >=0 && iy < ny && iz >=0 && iz < nz)
        {
            ++fieldcount[ix+nx*(iy+ny*iz)];
            flowfield[ix+nx*(iy+ny*iz)] += flowfielddata[i];
        } 
        else
        {
            //std::cout << "Error with index: " << ix << " " << iy << " " << iz << "\n";
        }
    }

    for (ix=0; ix<nx; ++ix)
    {
        coord[0] = ix * delta - halfDiameter-0.001;
        for(iy=0; iy<ny; ++iy)
        {
            coord[1] = iy * delta - height;
            for(iz=0; iz<nz; ++iz)
            {
                coord[2] = iz * delta - halfDiameter;
                //normalization
                if (fieldcount[ix+nx*(iy+ny*iz)] > 1)
                {
                    dummy = 1.0 / static_cast<double>(fieldcount[ix+nx*(iy+ny*iz)]);
                    flowfield[ix+nx*(iy+ny*iz)] *= dummy;
                }
                if (is_in_volume(coord, diameter, height, exitDiameter))
                {
                    tcount++;
                    if (fieldcount[ix+nx*(iy+ny*iz)]==0)
                    {
                        //std::cout << "No data point: " << ix << " " << iy << " " << iz << "\n";
                        ++countf;
                    } 
                }
            }
        }
    }

    std::cout << "Number of points without data in volume: " << countf  << " of " << tcount << "\n";

    emptyvoxel=true;
    while(emptyvoxel)
    {
        count=0;
        emptyvoxel = false;
        for (ix=0; ix<nx; ++ix)
        {
            for(iy=0; iy<ny; ++iy)
            {
                for(iz=0; iz<nz; ++iz)
                {
                    if (fieldcount[ix+nx*(iy+ny*iz)]==0)
                    {
                        emptyvoxel = true;
                        //edge cases near or in volume
                        if(ix>0 && fieldcount[ix-1+nx*(iy+ny*iz)]!=0){
                            flowfield[ix+nx*(iy+ny*iz)] += flowfield[ix-1+nx*(iy+ny*iz)];
                            ++fieldcount[ix+nx*(iy+ny*iz)];
                        }
                        if(iy>0 && fieldcount[ix+nx*(iy-1+ny*iz)]!=0){
                            flowfield[ix+nx*(iy+ny*iz)] += flowfield[ix+nx*(iy-1+ny*iz)];
                            ++fieldcount[ix+nx*(iy+ny*iz)];
                        }
                        if(iz>0  && fieldcount[ix+nx*(iy+ny*(iz-1))]!=0){
                            flowfield[ix+nx*(iy+ny*iz)] += flowfield[ix+nx*(iy+ny*(iz-1))];
                            ++fieldcount[ix+nx*(iy+ny*iz)];
                        }
                        if(ix<nx-1  && fieldcount[ix+1+nx*(iy+ny*iz)]!=0){
                            flowfield[ix+nx*(iy+ny*iz)] += flowfield[ix+1+nx*(iy+ny*iz)];
                            ++fieldcount[ix+nx*(iy+ny*iz)];
                        }
                        if(iy<ny-1  && fieldcount[ix+nx*(iy+1+ny*iz)]!=0){
                            flowfield[ix+nx*(iy+ny*iz)] += flowfield[ix+nx*(iy+1+ny*iz)];
                            ++fieldcount[ix+nx*(iy+ny*iz)];
                        }
                        if(iz<nz-1  && fieldcount[ix+nx*(iy+ny*(iz+1))]!=0){
                            flowfield[ix+nx*(iy+ny*iz)] += flowfield[ix+nx*(iy+ny*(iz+1))];
                            ++fieldcount[ix+nx*(iy+ny*iz)];
                        }

                        if(fieldcount[ix+nx*(iy+ny*iz)] > 0){
                            dummy = 1.0 / static_cast<double>(fieldcount[ix+nx*(iy+ny*iz)] );
                            flowfield[ix+nx*(iy+ny*iz)] *= dummy;
                            ++count;
                        }
                    } 
                }
            }
        }
        std::cout << "Number of points adapted from neighbors: " << count << "\n";
    }

    fileout.open(filename, std::ios::out);
    fileout.precision(7);
    fileout << "# " << diameter << " " << height << " " << " " << exitDiameter << " " << delta << " " << nx << " " << ny << " " << nz << "\n";
    for (ix=0; ix<nx; ++ix)
    {
        coord[0] = ix * delta - halfDiameter-0.001;
        for(iy=0; iy<ny; ++iy)
        {
            coord[1] = iy * delta - height;
            for(iz=0; iz<nz; ++iz)
            {
                coord[2] = iz * delta - halfDiameter;
                //fileout << ix <<' ' << iy << ' ' << iz << ' ' << " ";
                //fileout << coord[0] << ' ' << coord[1] << ' ' << coord[2] << "\n";
                fileout << flowfield[ix+nx*(iy+ny*iz)];
            }
        }
    }
    fileout.close();

    std::cout << "File " << filename << " written.\n";

    return 0;
}
