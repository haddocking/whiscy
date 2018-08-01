#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <map>
#include <algorithm>
#include <iostream>
using namespace std;
struct Residue {
  int nr;
  char code;
  float score;
};
struct Distance {
  short nr1;
  short nr2;
  float dis;
};

map<double, double> par;

inline double getweight(double val) {
  double weight;
  map<double,double>::iterator parpos = par.lower_bound(val);
  double vh = parpos->first;
  double wh = parpos->second;

  if (val == vh || parpos == par.begin()) return wh;
  else {
    parpos--;
    double vl = parpos->first;
    double wl = parpos->second;
    return wl + (wh - wl) * (val - vl)/(vh - vl);
  }
}

int compareres(Residue r1, Residue r2) {
  return (r1.score > r2.score);
}

int main(int argc, char *argv[]) {
  int n;
  char buf[500];
  bool err = 0;

  if (argc < 5) {
    fprintf(stderr, "ERROR: Too few arguments\n");
    fprintf(stderr, "Usage: parasmooth <surface conservation file> <low accessible conservation file> <residue distance matrix> <smoothing parameter file>\n");
    return 1;
  }
  FILE *surconsf = fopen(argv[1], "r");
  if (surconsf == NULL) {
    fprintf(stderr, "ERROR: Surface conservation file %s does not exist\n", argv[1]);
    fprintf(stderr, "Usage: parasmooth <surface conservation file> <low accessible conservation file> <residue distance matrix> <smoothing parameter file>\n");
    return(1);
  }
  Residue res_sur[10000];
  int res_surnr = 0;
  while (!feof(surconsf)) {
    Residue &r = res_sur[res_surnr];
    if (fgets(buf, 500, surconsf) != buf) break;
    if (sscanf(buf, "%f %c%d", &r.score, &r.code, &r.nr) != 3) {err = 1; break;}
    res_surnr++;
  }
  if (err) {
      fprintf(stderr, "ERROR: Reading error in surface conservation file %s\n", argv[1]);
    fprintf(stderr, "Usage: parasmooth <surface conservation file> <low accessible conservation file> <residue distance matrix> <smoothing parameter file>\n");
      return(1);
  }

  FILE *lacconsf = fopen(argv[2], "r");
  if (lacconsf == NULL) {
    fprintf(stderr, "ERROR: Low-accessible conservation file %s does not exist\n", argv[2]);
    fprintf(stderr, "Usage: parasmooth <surface conservation file> <low accessible conservation file> <residue distance matrix> <smoothing parameter file>\n");
    return(1);
  }
  Residue res_lac[10000];
  int res_lacnr = 0;
  while (!feof(lacconsf)) {
    Residue &r = res_lac[res_lacnr];
    if (fgets(buf, 500, lacconsf) != buf) break;
    if (sscanf(buf, "%f %c%d", &r.score, &r.code, &r.nr) != 3) {err = 1; break;}
    res_lacnr++;
  }
  if (err) {
      fprintf(stderr, "ERROR: Reading error in low-accessible conservation file %s\n", argv[2]);
    fprintf(stderr, "Usage: parasmooth <surface conservation file> <low accessible conservation file> <residue distance matrix> <smoothing parameter file>\n");
      return(1);
  }

  FILE *dismatf = fopen(argv[3], "r");
  if (dismatf == NULL) {
    fprintf(stderr, "ERROR: Residue distance matrix file %s does not exist\n", argv[3]);
    fprintf(stderr, "Usage: parasmooth <surface conservation file> <low accessible conservation file> <residue distance matrix> <smoothing parameter file>\n");
    return(1);
  }
  Distance *resdist = new Distance[2000000];
  int resdistnr = 0;
  while (!feof(dismatf)) {
    Distance &d = resdist[resdistnr];
    if (fgets(buf, 500, dismatf) != buf) break;
    if (sscanf(buf, "%hd %hd %f", &d.nr1, &d.nr2, &d.dis) != 3) {err = 1; break;}
    resdistnr++;
  }
  if (err) {
      fprintf(stderr, "ERROR: Reading error in distance matrix file %s\n", argv[3]);
    fprintf(stderr, "Usage: parasmooth <surface conservation file> <low accessible conservation file> <residue distance matrix> <smoothing parameter file>\n");
      return(1);
  }

  FILE *parfile = fopen(argv[4], "r");
  if (parfile == NULL) {
    fprintf(stderr, "ERROR: Smoothing parameter file %s does not exist\n", argv[4]);
    fprintf(stderr, "Usage: parasmooth <surface conservation file> <low accessible conservation file> <residue distance matrix> <smoothing parameter file>\n");
    return(1);
  }
  while (!feof(parfile)) {
    if (fgets(buf, 500, parfile) != buf) break;
    float v, w;
    if (sscanf(buf, "%f %f", &v, &w) != 2) {err = 1; break;}
    par[v] = w;
  }
  fclose(parfile);
  if (err) {
      fprintf(stderr, "ERROR: Reading error in smoothing parameter file %s\n", argv[5]);
    fprintf(stderr, "Usage: parasmooth <surface conservation file> <low accessible conservation file> <residue distance matrix> <smoothing parameter file>\n");
      return(1);
  }

  map<double,double>::iterator last = par.end(); last--;
  double maxdis = last->first;
  for (n = 0; n < res_surnr; n++) {
    double weight = 1;
    double score = res_sur[n].score;
    for (int i = 0; i < resdistnr; i++) {
      if (resdist[i].dis > maxdis) continue;
      int partner = -1;
      if (resdist[i].nr1 == res_sur[n].nr) partner = resdist[i].nr2;
      else if (resdist[i].nr2 == res_sur[n].nr) partner = resdist[i].nr1;
      if (!(partner + 1)) continue;
      bool found = 0;
      if (resdist[i].dis < maxdis) {
        for (int nn = 0; nn < res_surnr; nn++) {
      	  if (n == nn) continue;
      	  if (res_sur[nn].nr == partner) {
      	    found = 1;
      	    double currweight = getweight(resdist[i].dis);
      	    double currscore = currweight * res_sur[nn].score;
      	    weight += currweight;
      	    score += currscore;
      	    break;
      	  }
      	}
      }
      if (found) continue;
      if (resdist[i].dis < maxdis) {
        for (int nn = 0; nn < res_lacnr; nn++) {
      	  if (res_lac[nn].nr == partner) {
      	    double currweight = getweight(resdist[i].dis);
      	    double currscore = currweight * res_lac[nn].score;
      	    weight += currweight;
      	    score += currscore;
                  break;
      	  }
      	}
      }
    }
    res_sur[n].score = score/weight;
  }

  sort(res_sur, res_sur+n, compareres);

  for (n = 0; n < res_surnr; n++) {
    char id[5]; sprintf(id, "%c%d", res_sur[n].code, res_sur[n].nr);
    printf("%7.5f  %5s\n", res_sur[n].score, id);
  }

  return 0;
}
