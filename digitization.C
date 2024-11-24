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
                    if(j==1 && cluster_number(i,j)>2)
                    {break;}
                    else if(j>1 && cluster_number(i,j)!=1)
                    {break;}
                    else if(fired_pads(i,j)>2)
                    {break;}
                    

                    else if(j==max_layers && cluster_number(i,j) == 1)
                    {
                        events.push_back(i);
                    }
                }
        }
    return events;
}
bool check_mip_event(int i, int mip_layers) // Checks if a particular event satisfies the required condition for the selection of mip tracks
{
    for(int j=1; j<=mip_layers; j++)
                {
                    if(j==1 && cluster_number(i,j)>2)
                    {break;}
                    else if(j>1 && cluster_number(i,j)!=1)
                    {break;}
                    else if(j>1 && fired_pads(i,j)>2)
                    {break;}

                    else if(j==mip_layers && cluster_number(i,j) == 1)
                    {
                        return true;
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
std::pair<int,int> cluster_of_hit(int ii, int jj, int g) // the cluster to which the hit belongs
{
    std::pair<int,int> pr = {ii,jj};
    std::vector<cpoint_t> vec = cluster[pr];
    std::pair<int,int> xy = {0,0};

    for(auto &element : vec)
        {
            if(element.group==g)
            {
                xy = {element.x,element.y};
            }
        }
    return xy;
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
    double
    pad_size;
    std::string g4_file;
    std::string map_file;
    std::string line;
    double N_1, N_2, N_3, N_4;


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
            }
        }
    }

    config_file.close();

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

    double initial_val = -1000.0;
    //
    for(int i=0; i<max_hits; i++)
        {
            xHit[i]=initial_val;
            yHit[i]=initial_val;
            lHit[i]=initial_val;
            eHit[i]=initial_val;
        }
    //
    
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
    int cxpos;
    int cypos;
    int clusterId;

    tvec->Branch("eventId",&eventId,"eventId/L");
    tvec->Branch("layer",&layer,"layer/I");
    tvec->Branch("xpos",&xpos,"xpos/I");
    tvec->Branch("ypos",&ypos,"ypos/I");
    tvec->Branch("eBeam",&eBeam,"eBeam/D");
    tvec->Branch("pBeam",&pBeam,"pBeam/L");
    tvec->Branch("digit",&digit,"digit/L");
    tvec->Branch("cxpos",&cxpos,"cxpos/I");
    tvec->Branch("cypos",&cypos,"cypos/I");
    tvec->Branch("clusterId",&clusterId,"clusterId/I");


    //Digitization and Clustering

    std::map<std::pair<int,int>,std::vector<cpoint_t>> cluster_map_linkn;
    std::map<std::pair<int,int>,std::vector<hit_point_t>> cluster_map;


    std::map<std::pair<int,int>,std::vector<cpoint_t>> cluster_map_linkn_aux;
    std::map<std::pair<int,int>,std::vector<hit_point_t>> cluster_map_aux;

    std::vector<int> events; 

    // Adding Histograms

    std::vector<TH1F*> histograms;

    for (int i = 1; i <= 50; ++i) {
        std::string histName = "H_layer_" + std::to_string(i) + "_Cl_size"; // Generate unique name
        std::string name = "Cluster Sizes in Layer " + std::to_string(i);
        TH1F* hist = new TH1F(histName.c_str(), name.c_str(), 100, 0, 100); // Define the histogram
        histograms.push_back(hist); // Store it in the vector
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

                cluster_map[pair.first] = clhits;
                cluster_map_linkn[pair.first] = clust_linkn;

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
                        std::pair<int,int> xy = clust_aux.cluster_of_hit(pr.first,pr.second,grp);
                        cxpos = xy.first;
                        cypos = xy.second;
                        clusterId = grp;
                        tvec->Fill();
                    }
                unique_set.clear();
                vec.clear();
            }

            if(clust_aux.check_mip_event(i,mip_layers))
            {events.push_back(i);}
            
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
    clust.cluster = cluster_map_linkn;
    clust.hits = cluster_map;

    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed1 = end1 - start;
    std::cout << "Time taken: " << elapsed1.count() << " seconds" << std::endl;
    
    int Num_Clusters_cal = 0;
    /*
    for(int i=1; i<=max_layers; i++)
        {
            std::string hist_name = "layer_" + std::to_string(i) + "cluster_size";
            std::string hist_title = "Cluster size in layer " + std::to_string(i);
            int x = clust.max_cluster_size_layer(entries,i);
            TH1F *h1 = new TH1F(hist_name.c_str(),hist_title.c_str(),x,0,x);
            
            for(int j=0; j<entries; j++)
                {
                    std::pair<int,int> pr = {j,i};
                    for(auto &vec : clust.cluster_size(j,i))
                        {
                            h1->Fill(vec);
                        }
                }
        
            int count = h1->GetEntries();

            h2->SetBinContent(i,count);
            
            h1->Write();
            h1->Reset();
        }
    */
    for(int i=0; i<max_layers; i++)
        {
            histograms[i]->Write();
            int count = histograms[i]->GetEntries();
            h2->SetBinContent(i,count);
            
        }
    h2->Write();
    h3->Write();
    hx->Write();
    hy->Write();
    hxy->Write();
    output->Close();
    
    
    //clust.mip_event_clusters(entries,mip_layers); // Selection of mip like events
    //std::vector<int> events = clust.mip_events(entries,mip_layers); // vector containing the required event Ids

    //cout<<clust.fired_pads(26,2,0)<<endl;
    
    int cnt= events.size();

    /*
    for(auto &element : events) // printing out the selected event Ids
        {
            
            cnt++;
            //cout<<element<<", ";
        }
    */
    cout<<"No of mip like Events: "<<cnt<<endl;
    //clust.print_cpoint_t(10,12);
    
    /*
    for(auto &i : events) // Checking for any anomalies in selected events
        {
            for(int j=0; j<mip_layers; j++)
                {
                    if(clust.cluster_number(i,j)>2)
        cout<<clust.cluster_number(i,j)<< " " << i<<","<<j <<" "<<endl;
                }
        }
    */
    // Time to run the code
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken: " << elapsed.count() << " seconds" << std::endl;

}