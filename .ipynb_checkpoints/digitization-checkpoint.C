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
#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <chrono>
#include "LCClust.h"


void digitization()
{
    auto start = std::chrono::high_resolution_clock::now();

    LCClust cluster_make;
    
    fstream config_file, mapping_file;

    std::string config_file_name = "config.txt";

    //Loading config file into a map


    config_file.open(config_file_name,std::ios::in);

    double offset;
    double
    pad_size;
    std::string g4_file;
    std::string map_file;
    std::string line;

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
            }
        }
    }

    config_file.close();

    //Loading mapping file into a map

    mapping_file.open(map_file,std::ios::in);
    std::map<std::pair<int,int>,std::pair<int,int>> mapping_map;

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

    //Digitization


    for(int i=0; i<entries; i++)
        {
            tree->GetEntry(i);

            std::vector<std::vector<int>> vec; // long

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
                    }
                }
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

    output->Write();
    output->Close();

    // Time to sun the code
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken: " << elapsed.count() << " seconds" << std::endl;

}