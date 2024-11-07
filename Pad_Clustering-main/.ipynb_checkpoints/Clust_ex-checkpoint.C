#include "LCClust.h"

void Clust_ex(){

  LCClust   cluster_make;

  std::vector<hit_point_t>  clhits;
  std::vector<cpoint_t>     clust_linkn;
  std::map<int,  std::vector<int> > ev_clust_size;

  for (const auto &apixel : plane.hits) {
    int pxx, pxy;
    std::tie (pxx, pxy) = apixel;
    if (bad_pixels && bad_pixels->IsBad(planeid, pxx, pxy)) continue;
    lhist->FillHistW("tracking_planes_pixel_xy", planeid, pxx, pxy,
ev_weight);

    hit_point_t p1;
    p1.x = pxx;
    p1.y = pxy;
    p1.z = planeid;
    p1.e = 1.0;
    clhits.push_back(p1);
  }

  int n_cls = cluster_make.Clustering(clhits, clust_linkn,
LCClust::Link_Neighbours);  //clust_linkn are sorted


  for (const auto &clp : clust_linkn) {
 std::cout << "Cluster: " << clp.group << ",   Position:  " << clp.x
<< "  " << clp.y << std::endl;
    lhist->FillHistW("tracking_planes_cluster_xy", planeid, clp.x,
clp.y, ev_weight);
    std::vector<int> cluster_hits_ind;
    for (int hiti = 0; hiti < clhits.size(); ++hiti) {
      const auto &clht = clhits[hiti];
      if (clht.group == clp.group) {
//       std::cout << "Pixel: " << clht.x << ",  " << clht.y <<
std::endl;
         cluster_hits_ind.push_back(hiti);
      }
    }
    int clustsize = cluster_hits_ind.size();
    lhist->FillHistW("tracking_planes_cluster_size", planeid,
clustsize, ev_weight);
    if (clustsize==1) lhist->FillHistW("tracking_planes_one_pixel_cluster_xy", planeid, clp.x,
clp.y, ev_weight);

    //Fill histogram with big clusters
    if (clustsize > 25 && cluster_count < 20) {
      std::for_each(cluster_hits_ind.begin(), cluster_hits_ind.end(),
[&](const int jj) {
                    lhist->FillHistW("tracking_planes_cluster_pixels_xy", cluster_count,
clhits[jj].x, clhits[jj].y, ev_weight);});
      std::cout << "Event: " << ii << ":  cluster of " << clustsize <<" pixels in plane " << planeid << std::endl;
      ++cluster_count;
    }
    ev_clust_size[planeid].push_back(clustsize);
  }
}