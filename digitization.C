#include <TMath.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TEllipse.h>
#include <TH1F.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <cmath>    
#include <utility>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <unordered_map>
#include <variant>
#include <TStopwatch.h>
#include <chrono>
#include <algorithm>
#include <set>
#include "LCClust.h"


using namespace std;

class Cluster
{
    public:
        std::map<std::pair<int,int>,std::vector<cpoint_t>> cluster; //
        std::map<std::pair<int,int>,std::vector<hit_point_t>> hits;    
        std::map<int,std::vector<cpoint_t>> mip_cluster;
        
void print_cpoint_t(int ii,int jj) // prints all the clusters
{

    std::pair<int,int> pr = {ii,jj};
    
    cout<<"(EventId,Layer): ("<< ii<<", "<<jj<<")"<<endl;

    std::vector<cpoint_t> clust_vec = cluster[pr];
    cout<<"Clusters: "<<endl;
    for(auto &vec : clust_vec)
        {
            cout<<"("<< vec.x <<", "<< vec.y <<", "<< vec.group <<" , "<<vec.e<<"), ";
        }
    cout<<"No of Clusters: "<< cluster[pr].size()<<endl;
}

void print_hit_point_t(int ii,int jj) // prints all the hit points
{

    std::pair<int,int> pr = {ii,jj};
    
    cout<<"(EventId,Layer): ("<< ii<<", "<<jj<<")"<<endl;

    std::vector<hit_point_t> clust_vec = hits[pr];
    cout<<"hits: " <<endl;
    for(auto &vec : clust_vec)
        {
            cout<<"("<< vec.x <<", "<< vec.y <<", "<< vec.group <<"), ";
        }
    cout<<"No of Points: "<< hits[pr].size()<<endl;
}

int cluster_number(int ii,int jj) // Number of clusters in a particular event, layer
{
    std::pair<int,int> pr={ii,jj};
    return cluster[pr].size();
}

std::vector<int> cluster_size(int ii, int jj) // A vector of cluster sizes in each event and layer is returned
{
    std::pair<int,int> pr = {ii,jj};
    std::vector<int> cluster_sizes;
    for(auto &cpt : cluster[pr])
        {
            int count=0;
             for(auto &hpt : hits[pr])
                 {
                     if(hpt.group == cpt.group)
                     {
                         count++;
                     }
                 }
            cluster_sizes.push_back(count);
        }
    return cluster_sizes;
}
int cluster_size_group(int ii, int jj, int g) // Size of a particular cluster
{
    std::pair<int,int> pr = {ii,jj};
    int count = 0;
    for(auto &hpt : hits[pr])
        {
            if(hpt.group==g)
            {
                count++;
            }
        }
    return count;
}

int max_cluster_size(int events, int layers) // Maximum cluster size in the experiment
{
    std::vector<int> cluster_sizes;

    for(int i=0; i<events; i++)
        {
            for(int j=1; j<= layers; j++)
                {
                    std::pair<int,int> pr = {i,j};

                    for(auto &cpt : cluster[pr])
                    {
                        int count=0;
                        for(auto &hpt : hits[pr])
                         {
                             if(hpt.group == cpt.group)
                             {
                                 count++;
                             }
                         }
                        cluster_sizes.push_back(count);
                    }
                }
        }
                
        
    auto max_it = std::max_element(cluster_sizes.begin(),cluster_sizes.end());

    if(max_it != cluster_sizes.end())
    {return *max_it;}
    else{return 0;}

}

int max_cluster_size_layer(int ii,int jj) // max cluster size in each layer of an event
{
    int cluster_s = 0;
    
    for(int i = 0; i< ii; i++)
        {
            std::vector<int> clusts = cluster_size(i,jj);
            auto max_it = std::max_element(clusts.begin(),clusts.end());
            if(max_it != clusts.end())
            {
                auto k = *max_it;
                if(cluster_s<k)
                {
                    cluster_s = k;
                }
            }

        }
    return cluster_s;
}

int fired_pads(int ii, int jj)
{
    std::pair<int,int> pr = {ii,jj};
    std::vector<hit_point_t> hvec = hits[pr];
    std::set<std::pair<int,int>> pads;
    for(auto &hit : hvec)
        {
            std::pair<int,int> xy = {hit.x,hit.y};
            pads.insert(xy);
        }
    return pads.size();
}

std::vector<int> mip_events(int entries,int max_layers) // Event selection for mip tracks for a given number of layers to analyze
{
    std::vector<int> events;

    for(int i=0; i<entries; i++)
        {
            for(int j=1; j<=max_layers; j++)
                {
                    if(j==1 && cluster_number(i,j)>2 )
                    {break;}
                    else if(j==1 && cluster_number(i,j)==0)
                    {break;}
                    else if(j>1 && cluster_number(i,j)!=1)
                    {break;}
                    else if(fired_pads(i,j)>2)
                    {break;}
                    

                    else if(j==max_layers && cluster_number(i,j) == 1 )
                    {
                        events.push_back(i);
                    }
                }
        }
    return events;
}
bool check_mip_event(int i, int mip_layers, int test_layer) // Checks if a particular event satisfies the required condition for the selection of mip tracks
{
    for(int j=1; j<=mip_layers; j++)
        {
            if(j!= test_layer)
            {
                if(fired_pads(i,j)>2 || fired_pads(i,j)==0)
                {break;}
                else if(j==mip_layers && fired_pads(i,j)<3)
                {return true;}
            }
        }
    return false; 
}
std::map<int,std::vector<cpoint_t>> mip_event_clusters(int entries,int max_layers) // Retruns all the clusters in the mip selected events in a map with event id on as the key
{
    auto events = mip_events(entries,max_layers);
    
    for(auto &element : events)
        {
            std::vector<cpoint_t> aux_;
            for(int i=0; i<=max_layers;i++)
                {
                    std::pair<int,int> pr = {element,i};
                    for(auto &clt : cluster[pr])
                        {
                            aux_.push_back(clt);
                        }
                }
            mip_cluster[element] = aux_;
            aux_.clear();
        }
    return mip_cluster;
}
std::pair<double,double> cluster_of_hit(int ii, int jj, int g) // the cluster to which the hit belongs
{
    std::pair<int,int> pr = {ii,jj};
    std::vector<cpoint_t> vec = cluster[pr];
    std::pair<double,double> xy = {0.0,0.0};

    for(auto &element : vec)
        {
            if(element.group==g)
            {
                xy = {element.x,element.y};
            }
        }
    return xy;
}

std::vector<int> cluster_size_pads(int i,int j)
{
    std::pair<int,int> pr = {i,j};
    std::vector<int> sizes ;

    for(auto &element : cluster[pr])
        {
            int g = element.group;
            std::set<std::pair<int,int>> set;
            for(auto &hits : hits[pr])
                {
                    if(hits.group == g)
                    {
                    std::pair<int,int> x = {element.x, element.y};
                    set.insert(x);}
                }
            sizes.push_back(set.size());
            set.clear();
        }
    return sizes;
}

};

class Track
{
    public:
            
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
        TVectorD direction;
        TVectorD centroid;

    Track() : direction(3), centroid(3) {}

    void fit_line(int test_layer)
{
    double xMean = 0, yMean = 0, zMean = 0;
    int nPoints = x.size();
    int cnt=0;
    std::vector<int> index;
    for (int i = 0; i < nPoints; ++i) 
    {
        if(z[i] != test_layer){
        xMean += x[i];
        yMean += y[i];
        zMean += z[i];
        cnt += 1;
        index.push_back(i);
        }
    }
    xMean /= cnt;
    yMean /= cnt;
    zMean /= cnt;

    centroid(0) = xMean;
    centroid(1) = yMean;
    centroid(2) = zMean;
    
    TMatrixD covariance(3, 3);
    for (int i = 0; i < 3; ++i) 
    {
        covariance(i, i) += 1e-8; 
    }

    for(auto &i : index)
    {
            double dx = x[i] - xMean;
            double dy = y[i] - yMean;
            double dz = z[i] - zMean;
    
            covariance(0, 0) += dx * dx;
            covariance(0, 1) += dx * dy;
            covariance(0, 2) += dx * dz;
            covariance(1, 1) += dy * dy;
            covariance(1, 2) += dy * dz;
            covariance(2, 2) += dz * dz;
    }
    
    covariance(1, 0) = covariance(0, 1);
    covariance(2, 0) = covariance(0, 2);
    covariance(2, 1) = covariance(1, 2);

    TVectorD eigenValues(3);
    TMatrixD eigenVectors = covariance.EigenVectors(eigenValues);

    direction(0) = eigenVectors(0, 2); 
    direction(1) = eigenVectors(1, 2);
    direction(2) = eigenVectors(2, 2);
}

std::pair<double,double> get_xy(int test_layer)
{
    double t = (test_layer - centroid(2))/direction(2);

    double x = direction(0) * t + centroid(0);
    double y = direction(1) * t + centroid(1);
    std::pair<double,double> xy = {x,y};
    
    return xy;
}
bool is_straight(int test_layer)
{
    int n = x.size();
    if (n < 3) 
    {
        return true;
    }
    
    int ref_1, ref_2;

    for(auto &el : z)
        {
            if(el != test_layer)
            {
                ref_1 = el;
            }
        }
    for(auto &el : z)
        {
            if(el != test_layer && ref_1 != el)
            {
                ref_2 = el;
            }
        }
    
    double dx_ref = x[ref_2] - x[ref_1];
    double dy_ref = y[ref_2] - y[ref_1];
    double dz_ref = z[ref_2] - z[ref_1];

    for (int i = 0; i < n; i++) 
    {
        if(z[i] != ref_1 && z[i] != ref_2)
        {
            double dx = x[i] - x[ref_1];
            double dy = y[i] - y[ref_1];
            double dz = z[i] - z[ref_1];
    
            
            double cx = dy_ref * dz - dz_ref * dy;
            double cy = dz_ref * dx - dx_ref * dz;
            double cz = dx_ref * dy - dy_ref * dx;
    
            
            const double epsilon = 1e-9; 
            if (std::abs(cx) > epsilon || std::abs(cy) > epsilon || std::abs(cz) > epsilon) 
            {return false; }
        }
    }
    return true;
}
};

void digitization()
{
    auto start = std::chrono::high_resolution_clock::now();

    Cluster clust;

    LCClust cluster_make;
    
    fstream config_file, mapping_file;

    std::string config_file_name = "config.txt";

    //Loading config file into a map

    double offset;
    double pad_size;
    std::string g4_file;
    std::string map_file;
    std::string line;
    double N_1, N_2, N_3, N_4, initial_efficiency;
    int test_layer;

    config_file.open(config_file_name,std::ios::in);

    if (config_file.is_open()) {
        std::cout << "Config File loaded successfully." << std::endl;

        while (std::getline(config_file, line)) {
        std::istringstream ss(line);
        std::string key;

        // Read the key and the associated value
        if (ss >> key) {
            if (key == "offset") {
                ss >> offset;
            } else if (key == "pad_size") {
                ss >> pad_size;
            } else if (key == "g4_file") {
                ss >> g4_file;
            } else if (key == "mapping_file") {
                ss >> map_file;
            } else if (key == "N_1") {
                ss >> N_1;
            } else if (key == "N_2") {
                ss >> N_2;
            } else if (key == "N_3") {
                ss >> N_3;
            } else if (key == "N_4") {
                ss >> N_4;
            } else if (key == "test_layer" ) {
                ss >> test_layer;
            } else if (key == "initial_efficiency" ) {
                ss >> initial_efficiency;
            }
        }
    }

    config_file.close();

    if(N_1 + N_2 + N_3 + N_4 != 1)
    {
        N_1 /= N_1 + N_2 + N_3 + N_4;
        N_2 /= N_1 + N_2 + N_3 + N_4;
        N_3 /= N_1 + N_2 + N_3 + N_4;
        N_4 /= N_1 + N_2 + N_3 + N_4;
    }

    } else {
        std::cerr << "Failed to open the Config file." << std::endl;
        exit(1);    
    }


    //Loading mapping file into a map
    std::map<std::pair<int,int>,std::pair<int,int>> mapping_map;

    mapping_file.open(map_file,std::ios::in);

    if (mapping_file.is_open()) {
        std::cout << "Mapping File loaded successfully." << std::endl;
            
        while(1)
        {
            int chip_id, channel_id, xpad, ypad;
            mapping_file>>chip_id>>channel_id>>xpad>>ypad;
            
            std::pair<int,int> value = {chip_id,channel_id};
            std::pair<int,int> key = {xpad,ypad}; //Loading x and y of mapping file

            mapping_map[key] = value;
            if(mapping_file.eof()) break;
        }
    
    mapping_file.close();
        
    } else {
        std::cerr << "Failed to open the Mapping file." << std::endl;
        exit(1);
    }

    
    TFile *input = new TFile(g4_file.c_str(),"READ");

    TTree *tree = (TTree*)input->Get("uMSDHCALtree");

    TFile *output = new TFile("digitization.root","RECREATE");

    TTree *tvec = new TTree("tvec","tvec");

    
    //setting max no of hits in an event to extract data from array

    
    int entries = tree->GetEntries();
    int nHit;
    tree->SetBranchAddress("nHit",&nHit);

    int max_hits = 0;

    for(int i=0; i<entries;i++)
        {
            tree->GetEntry(i);
            
            if(max_hits<nHit)
            {
                 max_hits=nHit+1;   
            }
        }

    int max_layers = 50; // max layers in the experiment 
    int mip_layers = 8; // layers we analyze to look for mip like tracks
    int pBeam__;
    double eBeam__;
    double xHit[max_hits];
    double yHit[max_hits];
    int lHit[max_hits];
    double eHit[max_hits];
    
    tree->SetBranchAddress("pBeam",&pBeam__);
    tree->SetBranchAddress("eBeam",&eBeam__);
    tree->SetBranchAddress("xHit",xHit);
    tree->SetBranchAddress("yHit",yHit);
    tree->SetBranchAddress("lHit",lHit);
    tree->SetBranchAddress("eHit",eHit);


    // Clustering


    long eventId;
    int layer; // long
    int xpos; // long
    int ypos; // long
    double eBeam;
    long pBeam;
    long digit;
    double cxpos;
    double cypos;
    int clusterId;

    tvec->Branch("eventId",&eventId,"eventId/L");
    tvec->Branch("layer",&layer,"layer/I");
    tvec->Branch("xpos",&xpos,"xpos/I");
    tvec->Branch("ypos",&ypos,"ypos/I");
    tvec->Branch("eBeam",&eBeam,"eBeam/D");
    tvec->Branch("pBeam",&pBeam,"pBeam/L");
    tvec->Branch("digit",&digit,"digit/L");
    tvec->Branch("cxpos",&cxpos,"cxpos/D");
    tvec->Branch("cypos",&cypos,"cypos/D");
    tvec->Branch("clusterId",&clusterId,"clusterId/I");
    //Digitization and Clustering

    //std::map<std::pair<int,int>,std::vector<cpoint_t>> cluster_map_linkn;
    //std::map<std::pair<int,int>,std::vector<hit_point_t>> cluster_map;

    std::map<std::pair<int,int>,std::vector<cpoint_t>> cluster_map_linkn_aux;
    std::map<std::pair<int,int>,std::vector<hit_point_t>> cluster_map_aux;

    std::vector<int> events; 
    int selected_tracks=0;
    // Adding Histograms

    std::vector<TH1F*> histograms;
    std::vector<TH1F*> MIP_histograms;
    int max_cl_sizes[8] ={0,0,0,0,0,0,0,0};
    for (int i = 1; i <= max_layers; ++i) {
        std::string histName = "H_layer_" + std::to_string(i) + "_Cl_size"; 
        std::string name = "Cluster Sizes in Layer " + std::to_string(i);
        TH1F* hist = new TH1F(histName.c_str(), name.c_str(), max_hits, 0, max_hits); 
        histograms.push_back(hist); 
    }
    for (int i = 1; i<=mip_layers; i++)
        {
            std::string histName = "MIP_layer_" + std::to_string(i) + "_Cl_size"; 
            std::string name = "Cluster Sizes in MIP Layer " + std::to_string(i);
            TH1F* hist = new TH1F(histName.c_str(), name.c_str(),max_hits, 0, max_hits); 
            MIP_histograms.push_back(hist);
        }
        
    TH1F *h2 = new TH1F("Custers_layers", "Number of Clusters in each Layer", max_layers,0,max_layers);
    TH1F *h3 = new TH1F("Cluster_Energies", "Number of Clusters with different Energies",entries,0,entries);

    TH1F *hx =new TH1F("Clust_Dist_x", "Distribution of Cluster Position in X",50,0,50);
    TH1F *hy =new TH1F("Clust_Dist_y", "Distribution of Cluster Position in Y",50,0,50);
    TH2F *hxy = new TH2F("Clust_Dist_xy", "2D distribution of clusters",50,0,50,50,0,50);
    for(int i=0; i<entries; i++)
        {
            tree->GetEntry(i);

            std::pair<int,int> event_layer;
            std::map<std::array<int,3>,long> mp;

            std::map<std::pair<int,int>,std::vector<hit_point_t>> cluster_map_init;

            for (int j=0; j<nHit; j++)
                {
                    int x_aux = trunc((xHit[j] + offset)/pad_size);
                    int y_aux = trunc((yHit[j] + offset)/pad_size);
                    std::pair<int,int> aux_pair = {x_aux,y_aux}; // long
                    
                    if(mapping_map.count(aux_pair) > 0)
                    {
                        hit_point_t clhit;
                        event_layer = {i,lHit[j]};
                        clhit.x = x_aux;
                        clhit.y = y_aux;
                        clhit.z = lHit[j];
                        clhit.e = eHit[j];
                        std::array<int,3> arr = {clhit.x,clhit.y,lHit[j]};
                        mp[arr]++;
                        cluster_map_init[event_layer].push_back(clhit);
                    }
                }
            for(auto &pair : cluster_map_init)
            {
                std::vector<hit_point_t> clhits = pair.second;
                std::vector<cpoint_t> clust_linkn;

                int n_cls = cluster_make.Clustering(clhits, clust_linkn, LCClust::Link_Neighbours);

                //cluster_map[pair.first] = clhits;
                //cluster_map_linkn[pair.first] = clust_linkn;

                cluster_map_aux[pair.first] = clhits;
                cluster_map_linkn_aux[pair.first] = clust_linkn;
            
            }

            cluster_map_init.clear();

            Cluster clust_aux;
            clust_aux.cluster = cluster_map_linkn_aux;
            clust_aux.hits = cluster_map_aux;

            for(auto &pair : clust_aux.hits)
            {
                std::pair<int,int> pr = pair.first;
                std::vector<hit_point_t> vec = pair.second;
                std::set<hit_point_t> unique_set(vec.begin(), vec.end());

                for(auto &clts : unique_set)
                    {
                        eventId = pr.first;
                        layer = pr.second;
                        xpos = clts.x;
                        ypos = clts.y;
                        int grp = clts.group;
                        digit = mp[{clts.x,clts.y,pr.second}];
                        eBeam = eBeam__;
                        pBeam = pBeam__;
                        std::pair<double,double> xy = clust_aux.cluster_of_hit(pr.first,pr.second,grp);
                        cxpos = xy.first;
                        cypos = xy.second;
                        clusterId = grp;
                        tvec->Fill();
                    }
                unique_set.clear();
                vec.clear();
            }

            if(clust_aux.check_mip_event(i,mip_layers,test_layer))
            {
                events.push_back(i);
                Track line;
                for(int j=1; j<=mip_layers; j++)
                    {
                        std::pair<int,int> pr = {i,j};
                        for(auto &element : clust_aux.cluster[pr])
                            {
                                line.x.push_back(element.x);
                                line.y.push_back(element.y);
                                line.z.push_back(element.z);
                            }
                    }
                line.fit_line(test_layer);
                std::pair<double,double> xy = line.get_xy(test_layer);

                std::pair<int,int> pr = {i,test_layer};
                
                for(auto &element : clust_aux.cluster[pr])
                    {
                        double k = sqrt((xy.first-element.x)*(xy.first-element.x) + (xy.second-element.y)*(xy.second-element.y));
                        if(k<1.5 && line.is_straight(test_layer))
                        {
                            selected_tracks+=1;
                            for(int j=1; j<=mip_layers; j++)
                                {
                                    for(auto &clsize : clust_aux.cluster_size(i,j))
                                        {
                                            
                                            if(clsize>max_cl_sizes[j-1])
                                            {max_cl_sizes[j-1] = clsize;}
                                            
                                            MIP_histograms[j-1]->Fill(clsize);
                                        }
                                    
                                }
                            break;
                        }
                    }
                
            }
            
            for(int j=0; j<max_layers; j++)
                {
                    std::pair<int,int> pr = {i,j};
                    for(auto &vec : clust_aux.cluster[pr])
                        {
                            hx->Fill(vec.x);
                            hy->Fill(vec.y);
                            hxy->Fill(vec.x,vec.y);
                            h3->Fill(vec.e);
                        }
                }
            
            for(auto &element : clust_aux.cluster)
                {
                    int lay = element.first.second;
                    std::vector<int> clsz = clust_aux.cluster_size(i,lay);
                    for(auto &sz : clsz)
                        {
                            histograms[lay-1]->Fill(sz);
                        }
                }
            
            cluster_map_aux.clear();
            cluster_map_linkn_aux.clear();
            mp.clear();
        }

    tvec->Write();
    int min_clsz = max_cl_sizes[0];
    for (auto &x : max_cl_sizes)
        {
            if(x<min_clsz)
            {
                min_clsz = x;
            }
        }
    //clust.cluster = cluster_map_linkn;
    //clust.hits = cluster_map;
    double nu_1 = 0;
    double nu_2 = 0;
    double nu_3 = 0;
    double nu_4 = 0;

    for(int i=0; i<max_layers; i++)
        {
            histograms[i]->Write();
            int count = histograms[i]->GetEntries();
            h2->SetBinContent(i,count);
        }
    for(int i=0; i<mip_layers; i++)
        {
            
            MIP_histograms[i]->GetXaxis()->SetRangeUser(0,max_cl_sizes[i]);
            std::string s = "cluster_sizes_for_layer_" + std::to_string(i+1);
            TH1F *hist = (TH1F *)MIP_histograms[i]->Rebin(max_cl_sizes[i]/4 + 1,s.c_str());
            //TH1F *hist = (TH1F *)MIP_histograms[i];
            hist->GetXaxis()->SetRangeUser(0,max_cl_sizes[i]);
            hist->Scale( 1./hist->Integral());
            /*
            MIP_histograms[i]->Scale(1. / MIP_histograms[i]->Integral());
            MIP_histograms[i]->Write();
            TH1F *hist = MIP_histograms[i];
            */
            nu_1 += hist->GetBinContent(1);
            nu_2 += hist->GetBinContent(2);
            nu_3 += hist->GetBinContent(3);
            nu_4 += hist->GetBinContent(4);
            hist->Write();
        }
    

    nu_1 /= mip_layers;
    nu_2 /= mip_layers;
    nu_3 /= mip_layers;
    nu_4 /= mip_layers;

    double epsilon_0, epsilon_1, epsilon_2, epsilon_3;

    epsilon_0 = N_1 / nu_1;
    epsilon_1 = N_2 / nu_1 - nu_2*epsilon_0*epsilon_0/nu_1;
    epsilon_2 = N_3 / nu_1 - 2*nu_2*epsilon_0*epsilon_1/nu_1 - nu_3 * epsilon_0 * epsilon_0 *epsilon_0 / nu_1;
    epsilon_3 = N_4 / nu_1 - (nu_2/nu_1) * (2 * epsilon_2 * epsilon_0 + epsilon_0*epsilon_0) - 3 * nu_3 * epsilon_1 * epsilon_0 * epsilon_0 / nu_1 - nu_4 * epsilon_0 * epsilon_0 * epsilon_0 * epsilon_0 / nu_1;
        
    h2->Write();
    h3->Write();
    hx->Write();
    hy->Write();
    hxy->Write();
    output->Close();
    
    int cnt= events.size();
    fstream file;
    file.open("MIP_like_Events", ios::out);
    cout<<"nu_1: "<<nu_1<<", nu_2: "<<nu_2<<", nu_3: "<<nu_3<<", nu_4: "<<nu_4<<endl;
    cout<<"epsilon_0: "<<epsilon_0 << ", epsilon_1 : "<<epsilon_1 << ", epsilon_2 : "<<epsilon_2 << ", epsilon_3 : "<<epsilon_3 << endl;
    
    for(auto &element : events) // printing out the selected event Ids
        {
            file<<element<<endl;
        }
    
    cout<<"No of mip like Events: "<<cnt<<endl;
    cout<<"no of selected tracks: "<<selected_tracks<<endl;
    file.close();
    // Time to run the code
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken: " << elapsed.count() << " seconds" << std::endl;

}