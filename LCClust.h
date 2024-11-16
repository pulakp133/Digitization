// Author: sqytau
// Date: 2024-11-02
// This is implementation of clustering algorithms
//
// 1. Clustres are formed around pads with the maximum in the energy distribution. First, the pads are sorted,
//    then with the first pass each pad get a pointer on the neighbor with the highest energy if any, otherwise it points to itself.
//    In the second path the cluster are formed by following chains to pads with maximum energy.
//
// 2. Clustres are formed geometrically as islands.
//
// 3. K-Means
//

#if !defined LC_CLUST_H
#define LC_CLUST_H


#include <vector>
 

struct hit_point_t{
  hit_point_t(): x(0), y(0), z(0), e(0.0), group(0) {};
  hit_point_t& operator=(const hit_point_t &h) { if(&h == this) return *this; x=h.x; y=h.y, z=h.z; e=h.e; group=h.group; return *this; }

  int x, y, z; 
  double e; 
  int group; 
bool operator<(const hit_point_t& other) const {
        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
        if (z != other.z) return z < other.z;
        return group < other.group;
}
};

typedef hit_point_t *hpoint;


class cpoint_t {
  public:
  cpoint_t(): x(0.0), y(0.0), z(0.0), e(0.0), group(0) {};
//   cpoint_t(const cpoint_t &c): x(c.x), y(c.y), z(c.z), e(c.e), group(c.group) {};
  cpoint_t& operator=(const hit_point_t &h) {x=h.x; y=h.y, z=h.z; e=h.e; group=h.group; return *this;}
  double x, y, z, e; 
  int group;
};


typedef cpoint_t *cpoint;



template <class T> class point_t_cmp
{
public:  
  point_t_cmp(const std::vector<T> *ptr): p(ptr) {};
  bool operator() (size_t i, size_t j) const { if ((i < p->size()) && (j < p->size()) ) return ((p->at(i)).e > (p->at(j)).e); return false; }
  const std::vector<T> *p;
};


class LCClust {
public: 
  enum Algorithm {K_Means, E_Clustering, Link_Neighbours};
    
  LCClust(): fdebug(0) {};
  ~LCClust() {};
  
  int Clustering(std::vector<hit_point_t> &hits, std::vector<cpoint_t> &cls, const Algorithm alg = Link_Neighbours);
  int LinkClustering(std::vector<hit_point_t> &hits, std::vector<cpoint_t> &cls);
  int EClustering(std::vector<hit_point_t> &hits, std::vector<cpoint_t> &cls);
  int KMeansClustering(std::vector<hit_point_t> &hits, std::vector<cpoint_t> &cls, const int nclst = 2, const int ev_id=-1);

  int k_mean_testE(const int evid, std::vector<hit_point_t> &hits, std::vector<cpoint_t> &clust_reco);

  void print_eps(const int evid, hpoint pts, int len, cpoint cent, int n_cluster);

private:  
  double randf(double m);
  inline double dist2(cpoint a, hpoint b);
  inline int nearest(hpoint pt, cpoint cent, int n_cluster, double *d2);
  void kpp(hpoint pts, int len, cpoint cent, int n_cent);
  cpoint lloyd(hpoint pts, int len, int n_cluster);
  cpoint k_mean_test(const int evid, const int vsize, const int n_clust, hpoint v);
  
  int e_clustering(std::vector<hit_point_t> &hits, std::vector<cpoint_t> &cls, const std::vector<int> &pos);
  int link_clustering(std::vector<hit_point_t> &hits, std::vector<cpoint_t> &cls, const std::vector<int> &pos);
  
  void sort_index(const std::vector<hit_point_t> &hits, std::vector<int> &hit_i);
  
  int fdebug;
  
};


#endif
