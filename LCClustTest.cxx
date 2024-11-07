///*
// Author: sqytau
// Date: 2024-11-02
//  Description:
// This is test program to check the performance of clustering algorithms with simulated clusters
// Usage:
// g++  -o lcclust_test LCClustTest.cxx LCClust.C
//

#include <cmath>
#include <iostream>
#include <vector>
#include <chrono>

#include "LCClust.h"


double randf(double m)
{
  return m * rand() / (RAND_MAX - 1.);
}



void gen_xy_lc(std::vector<hit_point_t> &pt, std::vector<cpoint_t> &clt, const int clcount)
{
  const int pad_range = 63;
  const int sec_range = 10;
  const double e_max = 50000.0;
  double e;
  int ang, r;
  cpoint_t p;

/* note: this is not a uniform 2-d distribution */
  for (int ii = 0; ii < clcount; ++ii) {

    r = randf(pad_range/clcount) + ii*pad_range/clcount;
    ang = randf(sec_range);
    e = randf(e_max);
    p.e = e;
    p.y = r;
    p.x = ang;
    p.z = 0;
    p.group = ii;
    clt.push_back(p);

    for (int pi = -14; pi < 15; ++pi) {
      double rp = 1.8 * pi;
//       double rp = 1.8 * pi + randf(1.0);
      for (int si = -1; si < 2; ++si) {
        double sp = si * 15.0;
//         double sp = si * 15.0 + randf(1.0);
        if ( (r+pi) < 0.0 || (r+pi) > pad_range || (ang+si) < 0.0  || (ang+si) > sec_range) continue;
        double ei = e * exp(-(sp*sp + rp*rp)/30.0);
        if (ei < 5.0) continue;
        hit_point_t ph;
        ph.e = ei;
        ph.y = r+pi;
        ph.x = ang+si;
        ph.z = 0;
        ph.group = ii;
        pt.push_back(ph);
      }
    }
  }
}



void gen_xy(std::vector<hit_point_t> &pt, std::vector<cpoint_t> &clt, const int clcount)
{
  const int pad_range = 50;
  const int sec_range = 50;
  const double e_max = 50000.0;
  double e;
  int ang, r;
  cpoint_t p;

/* note: this is not a uniform 2-d distribution */
  for (int ii = 0; ii < clcount; ++ii) {

    r = randf(pad_range/clcount) + ii*pad_range/clcount;
    ang = randf(sec_range);
    e = randf(e_max);
    p.e = e;
    p.y = r;
    p.x = ang;
    p.z = 0;
    p.group = ii;
    clt.push_back(p);

    for (int pi = -14; pi < 15; ++pi) {
      double rp = 1.8 * pi;
//       double rp = 1.8 * pi + randf(1.0);
      for (int si = -14; si < 15; ++si) {
        double sp = si * 1.8;
//         double sp = si * 15.0 + randf(1.0);
        if ( (r+pi) < 0.0 || (r+pi) > pad_range || (ang+si) < 0.0  || (ang+si) > sec_range) continue;
        double ei = e * exp(-(sp*sp + rp*rp)/10.0);
        if (ei < 5.0) continue;
        hit_point_t ph;
        ph.e = ei;
        ph.y = r+pi;
        ph.x = ang+si;
        ph.z = 0;
        ph.group = ii;
        pt.push_back(ph);

      }
    }
  }
}



void dump_points(const std::vector<hit_point_t> &pt)
{
  int n = 0;
  for (std::vector<hit_point_t>::const_iterator itr = pt.begin(); itr != pt.end(); ++itr) {
    std::cout << "Point " << n << ":  "  << itr->group << " " << itr->x << "  " << itr->y << "  " << itr->e << std::endl;
    ++n;
  }
}



int main(int argc, char **argv)
{
  std::vector<hit_point_t> hits;
  std::vector<cpoint_t>    clust, clust_reco, clust_reco_e, clust_reco_k;
  const int n_test = 1;

   srand (time(NULL));
//   srand (123456);
//   srand (123451);
//   srand (123411);

  auto timer0 = std::chrono::high_resolution_clock::now();
  LCClust   cluster_make;

  for (int ii = 0; ii < n_test; ++ii) {
    gen_xy(hits, clust, 5);
//     gen_xy_lc(hits, clust, 5);
//     dump_points(hits);
    cluster_make.print_eps(0, &hits.front(), hits.size(), &clust.front(), clust.size());

    int n_cls = cluster_make.Clustering(hits, clust_reco, LCClust::Link_Neighbours);
    cluster_make.print_eps(1, &hits.front(), hits.size(), &clust_reco.front(), clust_reco.size());

    int n_cls_e = cluster_make.Clustering(hits, clust_reco_e, LCClust::E_Clustering);
    cluster_make.print_eps(2, &hits.front(), hits.size(), &clust_reco_e.front(), clust_reco_e.size());

    int n_cls_k = cluster_make.Clustering(hits, clust_reco_k, LCClust::K_Means);
//     cluster_make.print_eps(3, &hits.front(), hits.size(), &clust_reco_k.front(), clust_reco_k.size());

//     std::cout << std::endl << "Clustering : " << std::endl;
//     dump_points(hits);

    hits.clear();
    clust.clear();
    clust_reco.clear();
    clust_reco_e.clear();
    clust_reco_k.clear();
  }

  auto timer1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = timer1 - timer0;
  std::cout << "Number of tests: " << n_test << ";  time: " <<  elapsed.count() << " s." << std::endl;

  return 0;
}


