//Author: sqytau
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

#include <stdio.h>
// #include <stdlib.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <list>
#include <map>
#include <algorithm>

#include "LCClust.h"



int LCClust::Clustering(std::vector<hit_point_t> &hits, std::vector<cpoint_t> &cls, const Algorithm alg)
{
  int n_clust = -1;
  std::vector<cpoint_t>  clsus;
  switch (alg) {
    case Link_Neighbours: n_clust = LinkClustering(hits, clsus);   break;
    case E_Clustering:    n_clust = EClustering(hits, clsus);      break;
    case K_Means:         n_clust = KMeansClustering(hits, clsus); break;
  }
  if (n_clust > 0) {
    std::vector<int> hit_i;
    for (size_t i = 0; i < clsus.size(); ++i) hit_i.push_back(i);
    point_t_cmp<cpoint_t> p_cmp(&clsus);
    std::sort(hit_i.begin(), hit_i.end(), p_cmp);
    for (std::vector<int>::const_iterator i = hit_i.begin(); i != hit_i.end(); ++i) { cls.push_back(clsus[*i]); }
  }
  return n_clust;
}


int LCClust::LinkClustering(std::vector<hit_point_t> &hits, std::vector<cpoint_t> &cls)
{
  std::vector<int> hit_i;
  sort_index(hits, hit_i);
  return link_clustering(hits, cls, hit_i);
}


int LCClust::EClustering(std::vector<hit_point_t> &hits, std::vector<cpoint_t> &cls)
{
  std::vector<int> hit_i;
  sort_index(hits, hit_i);
  return e_clustering(hits, cls, hit_i);
}


int LCClust::KMeansClustering(std::vector<hit_point_t> &hits, std::vector<cpoint_t> &cls, const int nclst, const int ev_id)
{
  cpoint c = lloyd(&hits[0], hits.size(), nclst);  
  for (int i = 0; i < nclst; ++i) cls.push_back(c[i]);
  free(c);
  
  return cls.size();
}


void LCClust::sort_index(const std::vector<hit_point_t> &hits, std::vector<int> &hit_i)
{
  for (unsigned int i = 0; i < hits.size(); ++i) hit_i.push_back(i);
  point_t_cmp<hit_point_t> p_cmp(&hits);
  std::sort(hit_i.begin(), hit_i.end(), p_cmp);
}


double LCClust::randf(double m)
{
	return m * rand() / (RAND_MAX - 1.);
}


double LCClust::dist2(cpoint a, hpoint b)
{
	double x = a->x - b->x, y = a->y - b->y;
	return (x*x + y*y)*(b->e)*(b->e);
}


int LCClust::nearest(hpoint pt, cpoint cent, int n_cluster, double *d2)
{
	int i, min_i;
	cpoint c;
	double d, min_d;
 
#	define for_n for (c = cent, i = 0; i < n_cluster; i++, c++)
	for_n {
		min_d = HUGE_VAL;
		min_i = pt->group;
		for_n {
			if (min_d > (d = dist2(c, pt))) {
				min_d = d; min_i = i;
			}
		}
	}
	if (d2) *d2 = min_d;
	return min_i;
}


void LCClust::kpp(hpoint pts, int len, cpoint cent, int n_cent)
{
#	define for_len for (j = 0, p = pts; j < len; j++, p++)
	int i, j;
	int n_cluster;
	double sum, *d = (double*)malloc(sizeof(double) * len);
 
	hpoint p;
	cpoint c;
	cent[0] = pts[ rand() % len ];
	for (n_cluster = 1; n_cluster < n_cent; n_cluster++) {
		sum = 0;
		for_len {
			nearest(p, cent, n_cluster, d + j);
			sum += d[j];
			
		}
		sum = randf(sum);
		for_len {
		
			if ((sum -= d[j]) > 0) continue;
			cent[n_cluster] = pts[j];
			break;
		}
	} 

	for_len p->group = nearest(p, cent, n_cluster, 0);
	free(d);
}


cpoint LCClust::lloyd(hpoint pts, int len, int n_cluster)
{
	int i, j, min_i;
	int changed;
 
	cpoint cent = (cpoint)malloc(sizeof(cpoint_t) * n_cluster), c;
        hpoint p;
 
	/* assign init grouping randomly */
	//for_len p->group = j % n_cluster;
 
	/* or call k++ init */
	kpp(pts, len, cent, n_cluster);
 
	do {
		/* group element for centroids are used as counters */
		for_n { c->group = 0; c->x = c->y = 0; }
		for_len {
			c = cent + p->group;
			c->group++;
			c->x += p->x; c->y += p->y;
		}
		for_n { c->x /= c->group; c->y /= c->group; }
 
		changed = 0;
		/* find closest centroid of each point */
		for_len {
			min_i = nearest(p, cent, n_cluster, 0);
			if (min_i != p->group) {
				changed++;
				p->group = min_i;
			}
		}
	} while (changed > (len >> 10)); /* stop when 99.9% of points are good */
 
	for_n { c->group = i; }
 
	return cent;
}


void LCClust::print_eps(const int evid, hpoint pts, int len, cpoint cent, int n_cluster)
{
#	define W 400
#	define H 400
#	define for_n for (c = cent, i = 0; i < n_cluster; i++, c++)
#	define for_len for (j = 0, p = pts; j < len; j++, p++)

        int i, j;
	hpoint p;
    cpoint c;
	double min_x, max_x, min_y, max_y, scale, cx, cy;
	double min_e, max_e, scale_e, ce;
	double *colors = (double*)malloc(sizeof(double) * n_cluster * 3);

	for_n {
		colors[3*i + 0] = (3 * (i + 1) % 11)/11.;
		colors[3*i + 1] = (7 * i % 11)/11.;
		colors[3*i + 2] = (9 * i % 11)/11.;
	}

	max_x = max_y = max_e = -(min_x = min_y = min_e = HUGE_VAL);
	for_len {
		if (max_x < p->x) max_x = p->x;
		if (min_x > p->x) min_x = p->x;
		if (max_y < p->y) max_y = p->y;
		if (min_y > p->y) min_y = p->y;
		if (max_e < p->e) max_e = p->e;
		if (min_e > p->e) min_e = p->e;
	}
	scale = W / (max_x - min_x);
	if (scale > H / (max_y - min_y)) scale = H / (max_y - min_y);
	cx = (max_x + min_x) / 2;
	cy = (max_y + min_y) / 2;

        if (min_e != max_e) scale_e = 1.0/log(max_e-min_e+1.0e-5); else scale_e = 1.0;

        FILE *fp;
        std::stringstream ss("");
        ss << "lcclustering_ev_" << evid << ".eps";
        fp = fopen(ss.str().c_str(), "w+");

	fprintf(fp, "%%!PS-Adobe-3.0\n%%%%BoundingBox: -5 -5 %d %d\n", W + 10, H + 10);
	fprintf(fp,  "/l {rlineto} def /m {rmoveto} def\n"
		"/c { .25 sub exch .25 sub exch 3.5 0 360 arc fill } def\n"
		"/s { moveto -2 0 m 2 2 l 2 -2 l -2 -2 l closepath "
		"	gsave 1 setgray fill grestore gsave 3 setlinewidth"
		" 1 setgray stroke grestore 0 setgray stroke }def\n"
	);
	for_n {
		for_len {
			if (p->group != i) continue;
//                         ce = scale_e*log(p->e - min_e+1.0e-5);
                        ce = 1.0;
	        	fprintf(fp, "%g %g %g setrgbcolor\n",
			colors[3*i]*ce, colors[3*i + 1]*ce, colors[3*i + 2]);
			fprintf(fp, "%.3f %.3f c\n",
// 				(p->x - cx) * scale * 5.0 + W / 2,
				(p->x - cx) * scale + W / 2,
				(p->y - cy) * scale + H / 2);
		}
		fprintf(fp, "\n0 setgray %g %g s\n",
// 			(c->x - cx) * scale * 5.0 + W / 2,
			(c->x - cx) * scale + W / 2,
			(c->y - cy) * scale + H / 2);
	}
	fprintf(fp, "\n%%%%EOF");
	free(colors);
    fclose(fp);
#	undef for_n
#	undef for_len
}


cpoint LCClust::k_mean_test(const int evid, const int vsize, const int n_clust, hpoint v)
{
  cpoint c = lloyd(v, vsize, n_clust);
  return c;
}


//************** E Clustering Algorithm***************
int LCClust::e_clustering(std::vector<hit_point_t> &hits, std::vector<cpoint_t> &cls, const std::vector<int> &pos)
{
 // std::cout << std::endl;
  for (std::vector<int>::const_reverse_iterator itr = pos.rbegin(); itr != pos.rend(); ++itr) {
    int i = *itr; 
   // std::cout << i << " ";
    double emax = hits[i].e;
    hits[i].group = i;
    for (std::vector<int>::const_reverse_iterator itr1 = itr+1; itr1 != pos.rend(); ++itr1) {
      int j = *itr1;
//      if ( abs(hits[i].x - hits[j].x) <= 1.0 && abs(hits[i].y - hits[j].y) <= 1.0 ) {
      if ( abs(hits[i].x - hits[j].x) < 1.1 && abs(hits[i].y - hits[j].y) < 1.1 ) {
        if (hits[j].e >= emax) { 
          hits[i].group = j;  
          emax = hits[j].e;  
        }
      }
    }
  }
  
  if (fdebug > 2) {
    for (size_t i = 0; i < hits.size(); ++i) 
       std::cout << "Cluster: " << pos[i] << "   " << hits[pos[i]].e 
                 << "   " << hits[pos[i]].x << "   " << hits[pos[i]].y 
                 << "   " << hits[pos[i]].group << std::endl;
  }

  int clid = 0;
  cpoint_t pp;
  for (std::vector<int>::const_iterator itr = pos.begin(); itr != pos.end(); ++itr) {
    int i = *itr;  
    if (hits[i].group == i) {
      hits[i].group = clid;
      pp.x = hits[i].x * hits[i].e;
      pp.y = hits[i].y * hits[i].e;
      pp.e = hits[i].e;
      pp.group = hits[i].group;
      cls.push_back(pp);
      ++clid;
// std::cout << "New cluster: " << (clid) << "  point: " << i << "  " << hits[i].e << std::endl;     
    } else {
      int j = hits[i].group;
      int gg = hits[j].group;
// if(gg >= cls.size()) std::cout << "Wrong cluster id: " << gg << " point: " << i << std::endl; 
      hits[i].group = gg;
      cls[gg].x += hits[i].x * hits[i].e;
      cls[gg].y += hits[i].y * hits[i].e;
      cls[gg].e += hits[i].e;
    }
  }
  
  for (std::vector<cpoint_t>::iterator itr = cls.begin(); itr != cls.end(); ++itr) {
    itr->x /= itr->e;
    itr->y /= itr->e;
// std::cout << "Cluster: " << pc.group << "  position: " << itr->x << "  " << itr->y << "  " << pc.x << "  " << pc.y << std::endl;     
  }
  return clid;
}


//************** Linking neighbors Clustering Algorithm ***************
int LCClust::link_clustering(std::vector<hit_point_t> &hits, std::vector<cpoint_t> &cls, const std::vector<int> &pos)
{
  std::map <int, std::list<int> > clusters;
  for (std::vector<hit_point_t>::iterator itr = hits.begin(); itr != hits.end(); ++itr) itr->group = -1;
  int cl = 0;
  for (std::vector<int>::const_iterator itr = pos.begin(); itr != pos.end(); ++itr) {
    if (hits[*itr].group > -1 ) continue;
    clusters[cl].push_back(*itr);
    hits[*itr].group = cl;
    std::list<int>::iterator cpitr = clusters[cl].begin();
    while (cpitr != clusters[cl].end()) {
      int i = *cpitr;
//      std::cout << "link_clustering:  i = " << i << std::endl;
      for (std::vector<int>::const_iterator nitr = itr+1; nitr != pos.end(); ++nitr) {
        int j = *nitr;  
        if (hits[j].group > -1 ) continue;
        if ( abs(hits[i].x - hits[j].x) <= 1.1 && abs(hits[i].y - hits[j].y) <= 1.1 ) {
          clusters[cl].push_back(j);
          hits[j].group = cl;
//          std::cout << "link_clustering:  j = " << j << std::endl;
        }
      }
      ++cpitr;
    }
    ++cl;
  }
  
  // create vector of cluster points (actually it could be doen also inside the selection loop)
  cpoint_t clp;
  double px, py, pe;
  for (std::map <int, std::list<int> >::iterator itr = clusters.begin(); itr != clusters.end(); ++itr) {
    px = py = pe = 0.0;
    for (std::list<int>::iterator nitr = itr->second.begin(); nitr != itr->second.end(); ++nitr) {
      px += hits[*nitr].x * hits[*nitr].e;
      py += hits[*nitr].y * hits[*nitr].e;
      pe += hits[*nitr].e;
    }
    clp.x = px/pe;
    clp.y = py/pe;
    clp.e = pe;
    clp.group = itr->first;
    cls.push_back(clp);
  }

  return cls.size();
}

