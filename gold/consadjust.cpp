#include <cstdio>
#include <algorithm>
#include <functional>
#include <cmath>
#include <string>
#include <map>
using namespace std;

float z[25000];
map<char, double> resweight;

struct Residue {
  int nr;
  char res;
  string id;
  float score;
};

int residuecompare(Residue r1, Residue r2) {return (r1.score > r2.score);}

int main(int argc, char *argv[]) {
  int n;
  char buf[500];
  bool err = 0;
  if (argc < 4) {
    fprintf(stderr, "ERROR: Too few arguments\n");
    fprintf(stderr, "Usage: consadjust <conservation file> <residue weight file> <Z-table>\n");
    return 1;
  }
  char *fil = argv[1];
  FILE *f = fopen(fil, "r");
  if (f == NULL) {
    fprintf(stderr, "ERROR: Conservation file %s does not exist\n", fil);
    return(1);
  }

  char *wfil = argv[2];
  FILE *wf = fopen(wfil, "r");
  if (wf == NULL) {
    fprintf(stderr, "ERROR: Weight file %s does not exist\n", wfil);
    return(1);
  }
  for (n = 0; n < 20; n++) {
    if (fgets(buf, 500, wf) != buf) {err = 1; break;}
    char c; float weight;
    if (sscanf(buf, "%c %f", &c, &weight) != 2) {err = 1; break;}
    if ((c < 'A' || c > 'Y' || c == 'B' || c == 'J' || c == 'O' || c == 'U' || c == 'X' || c == 'Z')
     && (c < 'a' || c > 'y' || c == 'b' || c == 'j' || c == 'o' || c == 'u' || c == 'x' || c == 'z')) {
      err - 1;
      break;
    }

    if (c >='a' && c <= 'y') c -= 'a' - 'A';

    if (resweight.count(c)) {fprintf(stderr, "Warning: double weight specification\n");n--; continue;}
    resweight[c] = weight;
  }
  if (err) {
      fprintf(stderr, "ERROR: Reading error in weight file %s\n", wfil);
      return(1);
  }

  char *zfil = argv[3];
  FILE *zf = fopen(zfil, "r");
  if (zf == NULL) {
    fprintf(stderr, "ERROR: Z-table %s does not exist\n", zfil);
    return(1);
  }

  Residue *r = new Residue[100000];

  int consnr = 0;
  double sum = 0;
  double sumsq = 0;
  while (!feof(f)) {
    if (fgets(buf, 500, f) != buf) {break;}
    char id[500];
    if (sscanf(buf, "%f %s", &r[consnr].score, id) != 2) {err = 1; break;}
    r[consnr].id = id;
    r[consnr].res = id[0];
    sum += r[consnr].score;
    sumsq += r[consnr].score * r[consnr].score;
    consnr++;
  }
  if (err || !consnr) {
    fprintf(stderr, "ERROR: Reading error in conservation file %s\n", fil);
    return(1);
  }
  double stddev = sqrt((sumsq - sum*sum/consnr)/consnr);
  double mean = sum/consnr;
  if (!(stddev > 0)) {
    fprintf(stderr, "ERROR: All identical values in conservation file %s\n", fil);
    return(1);
  }

  double *zcalc = new double[consnr];
  for (n = 0; n < consnr; n++) {
    zcalc[n] = (r[n].score - mean) / stddev;
  }

  fgets(buf, 500, zf);
  err = 0;
  for (n = 0; n < 25000; n++) {
    if (!fgets(buf, 500, zf)) {err = 1; break;}
    if (sscanf(buf, "%f", &z[n]) != 1) {err = 1; break;}
  }
  fclose(zf);
  if (err) {
    fprintf(stderr, "ERROR: Reading error in Z-table %s\n", zfil);
    return(1);
  }

  for (n = 0; n < consnr; n++) {
    float pscore;

    bool neg = 0;
    float currz = zcalc[n];
    if (currz < 0) {currz *= -1; neg = 1;}
    float *upperz = upper_bound(z, z + 25000, currz, greater<float>());
    int index = upperz - z;

    do {
      if (index >= 24999) {pscore = .5; break;}
      if (index <= 0) {pscore = 0; break;}
      pscore = index - (currz - z[index]) / (z[index - 1] - z[index]);
    } while (0);

    bool newneg = 0;
    double pscore0 = pscore; if (neg) pscore0 = 50000 - pscore;
    double padj0 = 1000 * pscore0;
    if (resweight[r[n].res]) padj0 = pscore0 / resweight[r[n].res];
    double padj = padj0; if (padj > 25000) {padj = 50000 - padj; newneg = 1;}
    if (padj < 0) padj = 0;
    int pmin = (int) floor(padj);
    int pmax = (int) ceil(padj);
    double fac = padj - pmin;
    double znorm = fac * z[pmin] + (1-fac) * z[pmax];

    r[n].score = mean + znorm * (1-2*newneg) * stddev;
  }
  fprintf(stderr, "Subtracting average value ...\n");
  double scoresum = 0;
  for (n = 0; n < consnr; n++) scoresum += r[n].score;
  double scoreaverage = scoresum / consnr;
  for (n = 0; n < consnr; n++) r[n].score -= scoreaverage;

  fprintf(stderr, "Sorting scores...\n");
  sort(r, r + consnr, residuecompare);
  fprintf(stderr, "Writing scores...\n");
  for (n = 0; n < consnr; n++) {
    printf("%7.5f  %5s\n",
      r[n].score, r[n].id.c_str());
  }


  delete[] r, zcalc;
  fclose(f);
}
