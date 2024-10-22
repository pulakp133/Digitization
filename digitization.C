#include <iostream>
#include <fstream>
#include <cmath>    
#include <utility>
#include <string>
#include <vector>
#include <map>

void digitization()
{
    TFile *input = new TFile("uMSDHCAL_100ev.root","READ");

    TTree *tree = (TTree*)input->Get("uMSDHCALtree");

    TFile *output = new TFile("digitization.root","RECREATE");

    TTree *tvec = new TTree("tvec","tvec");

    int nHit;
    int pBeam__;
    double eBeam__;
    double xHit[5187];
    double yHit[5187];
    int lHit[5187];

    for(int i=0; i<5187; i++)
        {
            xHit[i]=1000.0;
            yHit[i]=1000.0;
            lHit[i]=1000.0;
        }

    tree->SetBranchAddress("nHit",&nHit);
    tree->SetBranchAddress("pBeam",&pBeam__);
    tree->SetBranchAddress("eBeam",&eBeam__);
    tree->SetBranchAddress("xHit",xHit);
    tree->SetBranchAddress("yHit",yHit);
    tree->SetBranchAddress("lHit",lHit);

    int val = 152235;

    int entries = tree->GetEntries();

    long EventId;
    long layer;
    long xpos;
    long ypos;
    double eBeam;
    long pBeam;
    long digit;

    tvec->Branch("EventId",&EventId,"EventId/L");
    tvec->Branch("layer",&layer,"layer/L");
    tvec->Branch("xpos",&xpos,"xpos/L");
    tvec->Branch("ypos",&ypos,"ypos/L");
    tvec->Branch("eBeam",&eBeam,"eBeam/D");
    tvec->Branch("pBeam",&pBeam,"pBeam/L");
    tvec->Branch("digit",&digit,"digit/L");
    

    for(int i=0; i<entries; i++)
        {
            tree->GetEntry(i);

            std::vector<std::vector<long>> vec;

            for (int j=0; j<nHit; j++)
                {
                    if(sqrt(pow(xHit[j],2) + pow(yHit[j],2)) <= 24)
                    {
                        long l = lHit[j];
                        long x = xHit[j];
                        long y = yHit[j];
                        
                        std::vector<long> hit = {l,x,y};
                        vec.push_back(hit);
                    }
                }
            std::map<std::vector<long>,long> m;

            int vec_size = vec.size();

            for(auto element : vec)
                {
                    m[element]++;
                }
            
            for(auto &pair : m)
                {
                    EventId = i;
                    std::vector<long> v_t = pair.first;
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
    
}