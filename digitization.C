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
#include "LCClust.h"

using namespace std;

class Cluster
{
    public:
        std::map<std::pair<int,int>,std::vector<cpoint_t>> cluster;
        std::map<std::pair<int,int>,std::vector<hit_point_t>> hits;    
        std::map<int,std::vector<cpoint_t>> mip_cluster;
        
void print_cpoint_t(int ii,int jj)
{

    std::pair<int,int> pr = {ii,jj};
    
    cout<<"(EventId,Layer): ("<< ii<<", "<<jj<<")"<<endl;

    std::vector<cpoint_t> clust_vec = cluster[pr];
    cout<<"Clusters: "<<endl;
    for(auto &vec : clust_vec)
        {
            cout<<"("<< vec.x <<", "<< vec.y <<", "<< vec.group <<"), ";
        }
    cout<<"No of Clusters: "<< cluster[pr].size()<<endl;
}

void print_hit_point_t(int ii,int jj)
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

int cluster_number(int ii,int jj)
{
    std::pair<int,int> pr={ii,jj};
    return cluster[pr].size();
}

std::vector<int> cluster_size(int ii, int jj)
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

int max_cluster_size(int events, int layers)
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

int max_cluster_size_layer(int ii,int jj)
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
std::vector<int> mip_event(int entries,int max_layers)
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

                    else if(j==max_layers && cluster_number(i,j) == 1)
                    {
                        events.push_back(i);
                    }
                }
        }
    return events;
}
std::map<int,std::vector<cpoint_t>> mip_event_clusters(int entries,int max_layers)
{
    auto events = mip_event(entries,max_layers);
    
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

    int max_layers = 50;
    int pBeam__;
    double eBeam__;
    double xHit[max_hits];
    double yHit[max_hits];
    int lHit[max_hits];
    double eHit[max_hits];

    double initial_val = -1000.0;
    
    for(int i=0; i<max_hits; i++)
        {
            xHit[i]=initial_val;
            yHit[i]=initial_val;
            lHit[i]=initial_val;
            eHit[i]=initial_val;
        }

    
    tree->SetBranchAddress("pBeam",&pBeam__);
    tree->SetBranchAddress("eBeam",&eBeam__);
    tree->SetBranchAddress("xHit",xHit);
    tree->SetBranchAddress("yHit",yHit);
    tree->SetBranchAddress("lHit",lHit);
    tree->SetBranchAddress("eHit",eHit);


    // Clustering


    long EventId;
    int layer; // long
    int xpos; // long
    int ypos; // long
    double eBeam;
    long pBeam;
    long digit;

    tvec->Branch("EventId",&EventId,"EventId/L");
    tvec->Branch("layer",&layer,"layer/I");
    tvec->Branch("xpos",&xpos,"xpos/I");
    tvec->Branch("ypos",&ypos,"ypos/I");
    tvec->Branch("eBeam",&eBeam,"eBeam/D");
    tvec->Branch("pBeam",&pBeam,"pBeam/L");
    tvec->Branch("digit",&digit,"digit/L");

    //Digitization and Clustering

    std::map<std::pair<int,int>,std::vector<cpoint_t>> cluster_map_linkn;
    std::map<std::pair<int,int>,std::vector<hit_point_t>> cluster_map;

    for(int i=0; i<entries; i++)
        {
            tree->GetEntry(i);

            std::vector<std::vector<int>> vec; // long
            std::pair<int,int> event_layer;

            std::map<std::pair<int,int>,std::vector<hit_point_t>> cluster_map_init;

            for (int j=0; j<nHit; j++)
                {
                    int x_aux = round((xHit[j] + offset)/pad_size);
                    int y_aux = round((yHit[j] + offset)/pad_size);
                    std::pair<int,int> aux_pair = {x_aux,y_aux}; // long
                    
                    if(mapping_map.count(aux_pair) > 0)
                    {
                        int l = lHit[j]; // long
                        int x = aux_pair.first; // long
                        int y = aux_pair.second; // long
                        
                        std::vector<int> hit = {l,x,y};
                        vec.push_back(hit);

                        hit_point_t clhit;
                        event_layer = {i,lHit[j]};
                        clhit.x = round((xHit[j] + offset)/pad_size);
                        clhit.y = round((yHit[j] + offset)/pad_size);
                        clhit.z = lHit[j];
                        clhit.e = eHit[j];

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
            
            }

            cluster_map_init.clear();
            
            std::map<std::vector<int>,long> m; // long

            for(auto element : vec)
                {
                    m[element]++;
                }
            
            for(auto &pair : m)
                {
                    EventId = i;
                    std::vector<int> v_t = pair.first; // long
                    layer = v_t[0];
                    xpos = v_t[1];
                    ypos = v_t[2];
                    digit = pair.second;
                    eBeam = eBeam__;
                    pBeam = pBeam__;
                    tvec->Fill();
                }
            vec.clear();
            m.clear();
        }

    tvec->Write();
    clust.cluster = cluster_map_linkn;
    clust.hits = cluster_map;

    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed1 = end1 - start;
    std::cout << "Time taken: " << elapsed1.count() << " seconds" << std::endl;
    // Adding Histograms

    TH1F *h2 = new TH1F("Custers_layers", "Number of Clusters in each Layer", max_layers,0,max_layers);

    TH1F *hx =new TH1F("Clust_Dist_x", "Distribution of Cluster Position in X",50,0,50);
    TH1F *hy =new TH1F("Clust_Dist_y", "Distribution of Cluster Position in Y",50,0,50);
    TH2F *hxy = new TH2F("Clust_Dist_xy", "2D distribution of clusters",50,0,50,50,0,50);

    TH1F *h3 = new TH1F("Cluster_Events", "Number of Clusters in each event",entries,0,entries);

    int Num_Clusters_cal = 0;
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
                    for(auto &vec : clust.cluster[pr])
                        {
                            hx->Fill(vec.x);
                            hy->Fill(vec.y);
                            hxy->Fill(vec.x,vec.y);
                        }
                }
            int count = h1->GetEntries();

            h2->SetBinContent(i,count);

            h1->Write();
            h1->Reset();
        }
    for(int i=0; i<entries; i++)
        {
            int count = 0;
            
            for (int j=1; j<=max_layers; j++)
                {
                    std::pair<int,int> pr = {i,j};
                    count += clust.cluster[pr].size();
                }
            h3->SetBinContent(i,count);
        }

    h2->Write();
    h3->Write();
    hx->Write();
    hy->Write();
    hxy->Write();
    output->Close();

    clust.mip_event_clusters(entries,max_layers);
    std::map<int,std::vector<cpoint_t>> mip_clust = clust.mip_cluster;
/*
    int cnt=0;
    for(auto &element : events)
        {
            
            cnt++;
            cout<<element<<", ";
        }
    cout<<"No of such Events: "<<cnt<<endl;
*/

    // Time to run the code
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken: " << elapsed.count() << " seconds" << std::endl;

}