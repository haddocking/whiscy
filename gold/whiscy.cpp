#include <cstdio>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
using namespace std;

extern void pamInit();
extern int pamLoadSequences(const char *seqfile, const char *distfile, int *retseqnr, int *retseqlen, char **refseq);
extern int pamCalcSimilarity(int pos, float *distances, float *scores);

void loadlist(FILE *f, vector<int> &list, const char *type, const char *filename) {
  while (!feof(f)) {
    char buf[2000];
    char *ret = fgets(buf, 2000, f);
    if (ret != buf) break;
    int num = atoi(buf);
    if (num <= 0) {
      fprintf(stderr, "ERROR: Reading error in %s list %s\n", type, filename);
      exit(1);
    }
    list.push_back(num);
  }
  fclose(f);
}

struct Residue {
  int nr;
  float score;
};

int residuecompare(Residue r1, Residue r2) {return (r1.score > r2.score);}

int main(int argc, char *argv[]) {
  int n;
  if (argc < 5) {
    fprintf(stderr, "ERROR: Too few arguments\n");
    fprintf(stderr, "Usage: whiscy <surface list> <conversion table> <alignment file> <distance file>\n");
    return 1;
  }
  fprintf(stderr, "Parsing surface list... \n");

  FILE *surfile = fopen(argv[1], "r");
  if (surfile == NULL) {
    fprintf(stderr,"ERROR: Surface file does not exist\n");
    fprintf(stderr, "Usage: whiscy <surface list> <conversion table> <alignment file> <distance file>\n");
    return 1;
  }
  vector<int> surlist;
  loadlist(surfile, surlist, "surface", argv[1]);
  
  #ifdef DEBUG
  for (vector<int>::const_iterator i = surlist.begin(); i != surlist.end(); ++i)
    cout << *i << ' ';
  cout << endl;
  #endif

  fprintf(stderr, "Loading conversion table...\n");

  FILE *convfile = fopen(argv[2], "r");
  if (convfile == NULL) {
    fprintf(stderr,"ERROR: Conversion table does not exist\n");
    fprintf(stderr, "Usage: whiscy <surface list> <conversion table> <alignment file> <distance file>\n");
    return 1;
  }
  map<int, int> conv;
  char buf[2000];
  if (!feof(convfile)) fgets(buf, 2000, convfile);
  while (!feof(convfile)) {
    fgets(buf, 2000, convfile);
    int from, to;
    if (sscanf(buf, "%d %d", &from, &to) != 2) {
      printf("%s", buf);
      fprintf(stderr, "ERROR: Reading error in conversion table %s\n", argv[2]);
      return 1;
    }
    conv[from] = to;
  }

  #if DEBUG
  for(map<int,int>::const_iterator it = conv.begin(); it != conv.end(); ++it)
    cout << it->first << ": " << it->second << ", ";
  cout << endl;
  #endif

  fprintf(stderr, "Converting...\n");

  for (n = 0; n < surlist.size(); n++) {
    if (!conv.count(surlist[n])) {
      fprintf(stderr, "WARNING: Surface residue number %d cannot be converted\n", surlist[n]);
      fprintf(stderr, "Continuing program...\n");
      surlist.erase(surlist.begin() + n);
      n--;
      continue;
    }
    surlist[n] = conv[surlist[n]];
  }

  #if DEBUG
  for (vector<int>::const_iterator i = surlist.begin(); i != surlist.end(); ++i)
    cout << *i << ' ';
  cout << endl;
  #endif

  int sum = surlist.size();
  if (!sum) {fprintf(stderr, "ERROR: No surface residues\n"); return 1;}

  Residue *totlist = new Residue[sum];
  fprintf(stderr, "Initializing score calculation...\n");

  pamInit();
  int seqnr, seqlen;
  char *refseq;
  int err = pamLoadSequences(argv[3], argv[4], &seqnr, &seqlen, &refseq);

  #if DEBUG
  cout << seqnr << endl;
  cout << seqlen << endl;
  cout << refseq << endl;
  #endif

  switch(err) {
    case 0: break;
    case 1: fprintf(stderr, "ERROR: Distance file %s does not exist\n", argv[4]); return 1;
    case 2: fprintf(stderr, "ERROR: Sequence file %s does not exist\n", argv[3]); return 1;
    case 3: fprintf(stderr, "ERROR: Invalid number of sequences\n"); return 1;
    case 4: fprintf(stderr, "ERROR: Invalid sequence length\n"); return 1;
    case 5: fprintf(stderr, "ERROR: Reading error in distance file %s\n", argv[4]); return 1;
    case 6: fprintf(stderr, "ERROR: Reading error in sequence file %s\n", argv[3]); return 1;
    default: fprintf(stderr, "ERROR: Unknown error\n"); return 1;
  }

  fprintf(stderr, "Calculating scores...\n");
  float *scores = new float[seqnr];
  float *distances = new float[seqnr];
  int realsum = 0;
  
  for (n = 0; n < sum; n++) {
    Residue &r = totlist[realsum];
    r.nr = surlist[n];
    r.score = - 1000;
    if (surlist[n] > seqlen || surlist[n] < 1) {
      fprintf(stderr, "ERROR: surface residue out of range\n");
      return 1;
    }
    int posnr = pamCalcSimilarity(surlist[n]-1, distances, scores);
    if (posnr <= 0 || posnr > seqnr) continue;
    realsum++;
    r.score = scores[posnr - 1];
  }

  if (!realsum) {fprintf(stderr, "ERROR: No sequence information for any surface residues\n"); return 1;}

  fprintf(stderr, "Subtracting average value ...\n");
  double scoresum = 0;
  for (n = 0; n < realsum; n++) scoresum += totlist[n].score;
  double scoreaverage = scoresum / realsum;
  for (n = 0; n < realsum; n++) totlist[n].score -= scoreaverage;

  fprintf(stderr, "Sorting scores...\n");
  sort(totlist, totlist+realsum, residuecompare);
  fprintf(stderr, "Writing scores...\n");
  for (n = 0; n < realsum; n++) {
    Residue &r = totlist[n];
    char id[5]; sprintf(id, "%c%d", refseq[r.nr-1], r.nr);
    printf("%7.5f  %5s\n",
      r.score, id);
  }

  delete [] totlist;
  return 0;
}

