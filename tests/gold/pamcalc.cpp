#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <ctime>
#include <cstring>
#include <iostream>
#include <iomanip>
using namespace std;

typedef double (matrix)[20][20];

int code[256];

struct Distance {
  int seq;
  float dist;
  matrix mat;
  double expect[20];
};

Distance *dis;
int *seqtodis;

int **seq;
int seqnr;
int seqlen;

int DistCompare(const void *pd1, const void *pd2) {
  const Distance *d1 = (const Distance *) pd1;
  const Distance *d2 = (const Distance *) pd2;

  if (d1->dist < d2->dist)
      return -1;
   else if (d1->dist > d2->dist)
      return 1;
   else
      return 0;

  /* Old implementation. Works fine in Linux (tested on Ubuntu),
   * but behavior is different on OSX
   */
  // return (d1->dist > d2->dist);
}

void pamInit() {
  for (int n = 0; n < 256; n++) code[n] = -1;
  int nr = -1;
  nr++ ;code['A'] = nr;
  nr++ ;code['R'] = nr;
  nr++ ;code['N'] = nr;
  nr++ ;code['D'] = nr;
  nr++ ;code['C'] = nr;
  nr++ ;code['Q'] = nr;
  nr++ ;code['E'] = nr;
  nr++ ;code['G'] = nr;
  nr++ ;code['H'] = nr;
  nr++ ;code['I'] = nr;
  nr++ ;code['L'] = nr;
  nr++ ;code['K'] = nr;
  nr++ ;code['M'] = nr;
  nr++ ;code['F'] = nr;
  nr++ ;code['P'] = nr;
  nr++ ;code['S'] = nr;
  nr++ ;code['T'] = nr;
  nr++ ;code['W'] = nr;
  nr++ ;code['Y'] = nr;
  nr++ ;code['V'] = nr;
  nr = -1;
  nr++ ;code['a'] = nr;
  nr++ ;code['r'] = nr;
  nr++ ;code['n'] = nr;
  nr++ ;code['d'] = nr;
  nr++ ;code['c'] = nr;
  nr++ ;code['q'] = nr;
  nr++ ;code['e'] = nr;
  nr++ ;code['g'] = nr;
  nr++ ;code['h'] = nr;
  nr++ ;code['i'] = nr;
  nr++ ;code['l'] = nr;
  nr++ ;code['k'] = nr;
  nr++ ;code['m'] = nr;
  nr++ ;code['f'] = nr;
  nr++ ;code['p'] = nr;
  nr++ ;code['s'] = nr;
  nr++ ;code['t'] = nr;
  nr++ ;code['w'] = nr;
  nr++ ;code['y'] = nr;
  nr++ ;code['v'] = nr;
}

void AssemblePam(double distance, matrix &m) {
 //Returns PAM(distance) matrix: m[i][j] is the probability of i being replaced by j
  static double logpameigval[20] = {-0.022091252,-0.019297602, 0.000004760,-0.017477817,
                           -0.016575549,-0.015504543,-0.002112213,-0.002685727,
                           -0.002976402,-0.013440755,-0.012926992,-0.004293227,
                           -0.005356688,-0.011064786,-0.010480731,-0.008760449,
                           -0.007142318,-0.007381851,-0.007806557,-0.008127024};
  static matrix pameigvec = {
  { -0.173697, -0.885673, -1.0049, 0.286697,
     1.88382, 1.83495, -0.198369, 0.197182,
    -0.328623, 0.46536, -0.29925, -0.0965166,
    -0.74739, -0.854738, 0.678401, -0.183266,
     0.192614, -0.482818, 0.146535, 0.0798611},
  { -0.184354, -0.175916, -1.0076, -0.142763,
    -0.0628839, 0.0139009, 0.141204, 0.331822,
    -0.539164, -0.708375, 0.0817457, -0.37228,
     2.4085, -2.80149, -2.21178, 1.28384, 1.42798,
     0.0188545, -0.4161, -0.332652},
  { -3.40863, -0.900981, -0.994619, 2.41327,
    -0.811826, -0.630945, -0.174313, 0.293746,
    -0.383093, -1.26984, 0.89058, -0.481288,
     0.12826, 0.0590369, 0.799897, -0.296714,
    -0.552095, 0.384358, 0.177192, 0.169113},
  {  2.46339, -1.69023, -0.992457, 0.277099,
    -0.415795, -0.0669855, -0.230705, 0.389199,
    -0.543538, -1.76778, 0.892756, -0.58089,
    -0.185477, -0.168277, -0.685745, -1.4301,
    -1.48203, 0.838871, 0.35372, 0.654423},
  { -0.0653222, -0.142554, -1.00875, -0.0956775,
     0.00982116, 0.0114652, -0.267974, -5.26286,
    -1.24758, 0.00651089, -0.0117717, 0.384433,
     0.233031, 0.151659, -0.112617, -0.07216,
    -0.130136, 0.277631, -0.271657, 0.0955676},
  {  0.411745, -1.0016, -0.993702, 1.55374,
    -0.744468, -0.843633, -0.175287, 0.385097,
    -0.406764, 3.44619, -1.97631, -0.351715,
     0.735297, 0.188087, -0.904437, 0.289089,
    -1.65422, 0.177167, -0.146884, 0.250162},
  { -1.46727, 2.10777, -0.9997, -1.83253,
     0.527169, 0.327848, -0.228532, 0.39569,
    -0.510345, 0.172614, -0.202514, -0.508889,
    -0.0845795, -0.16422, -1.27133, -1.46079,
    -1.68406, 0.685183, 0.331996, 0.684302},
  { -0.0358926, -0.101557, -1.0116, -0.332503,
    -0.196083, -0.33719, -0.250326, 0.316732,
    -0.566645, 0.055228, -0.0474809, -0.496718,
    -1.85709, 0.644994, -0.657107, 0.980698,
     0.993692, 1.38279, -0.619067, -0.747332},
  {  0.418693, 0.437099, -0.995008, -0.951753,
     0.506729, 0.416234, -0.116331, 0.255448,
    -0.181811, -0.851382, 0.440693, -0.67166,
     1.20247, 0.832078, 1.74873, 3.77201,
    -2.59574, -0.020501, -0.294246, -0.0926676},
  {  0.125546,-0.0372708, -0.995249, 0.281829,
     2.9408, -3.14634, -0.100311, 0.0550663,
     0.61025, -0.178163, 0.40751, 1.30446,
    -0.137214, 0.576156, -0.721845, 0.624154,
     0.424991, -0.423803, 1.73362, 1.19383},
  {  0.0180848, 0.0195127, -0.986272, -0.0422427,
    -0.0639326, 0.0913046, -0.0163454, 0.163283,
     1.15374, 0.263141, 0.678789, 2.26955,
     0.201714, -0.32496, 0.14775, -0.423257,
    -0.587479, 0.62248, -1.06261, -1.22256},
  {  0.254528, 0.0789705, -1.00771, -0.252216,
     0.231033, 0.338253, -0.122814, 0.41346,
    -0.510865, 0.460223, 0.425877, -0.347116,
     1.89864, 1.78614, 0.739422, -0.88218,
     1.38203, 0.548471, 0.111328, -0.142583},
  { -0.106978, -0.00225986, -1.01125, 0.0638037,
    -0.139527, 0.163789, -0.0895913, 0.233315,
     0.587038, -3.85451, -6.80573, 1.53798,
     0.593227, 1.02989, 0.45206, -0.415565,
     0.350157, 0.343607, 0.0730036, -0.369488},
  { -0.0168964, -0.0313974, -0.978517, -0.0355775,
    -0.186282, 0.218407, 0.470191, -0.473411,
     2.96084, 0.0480287, -0.0359644, -1.19118,
    -0.103219, 0.0888554, 0.0569158, 0.0670114,
     0.632513, 0.252839, -2.02657, 2.8489},
  { -0.0775297, -0.133835, -1.0047, -0.269255,
    -0.256457, -0.335363, -0.204311, 0.239236,
    -0.470702, -0.430193, 0.247265, -0.210275,
    -0.372086, 1.15706, -0.997508, -0.372601,
    -0.139173, -3.46403, -1.56201, -0.65},
  {  1.20632, 2.50065, -1.00745, 1.62844,
    -0.315781, -0.385228, -0.130661, 0.0778331,
    -0.422256, -0.0873078, 0.00251954, -0.279167,
    -0.424892, -0.931607, 1.05355, -0.221071,
     0.328534, -0.328235, -0.059254, -0.0121199},
  { -0.203699, -0.921196, -1.00239, -2.07912,
    -1.27451, -1.41652, -0.173309, 0.172953,
    -0.268621, 0.492529, -0.271765, -0.0140023,
    -0.313962, -1.39104, 1.95919, -0.381254,
     0.433291, -0.550867, 0.509705, 0.279679},
  { -0.0320838, -0.0845875, -1.09237, -0.0587499,
     0.0221192, 0.0126683, 10.1156, 0.135866,
    -1.72557, 0.0511482, -0.0169357, 0.554823,
    -0.563199, 0.343121, 0.15316, -0.115578,
    -0.285916, 0.0192524, 0.0998362, 0.102109},
  {  0.0694741, -0.0120275, -0.960742, -0.0175369,
     0.0378244, -0.0709755, 0.556068, -0.902598,
     2.98718, 0.0476938, -0.0730921, -2.46403,
     0.0035873, -0.032857, -0.324921, -0.394633,
    -0.27464, -0.374326, 2.21, -3.01709},
  {  0.00674906, 0.21934, -0.996274, 0.0661936,
    -1.71762, 1.4936, -0.151563, 0.0713457,
     0.355119, 0.140147, 0.366456, 1.27468,
    -0.40718, 0.563689, -0.771383, 0.858587,
     0.42959, -0.624028, 2.06336, 1.02797}
  };

  static matrix pameigvecinv =
  {
      {-0.01522976,-0.00746819,-0.13934468, 0.11755315,-0.00212101,
       0.01558456,-0.07408235,-0.00322387, 0.01375826, 0.00448826,
       0.00154174, 0.02013313,-0.00159183,-0.00069275,-0.00399898,
       0.08414055,-0.01188178,-0.00029870, 0.00220371, 0.00042546},
      {-0.07765582,-0.00712634,-0.03683209,-0.08065755,-0.00462872,
       -0.03791039,0.10642147,-0.00912185,0.01436308,-0.00133243,
        0.00166346,0.00624657,-0.00003363,-0.00128729,-0.00690319,
        0.17442028,-0.05373336,-0.00078751,-0.00038151,0.01382718},
      {-0.08810973,-0.04081786,-0.04066004,-0.04736004,-0.03275406,
      -0.03761164,-0.05047487, -0.09086213,-0.03269598,-0.03558015,
      -0.08407966,-0.07970977, -0.01504743,-0.04011920,-0.05182232,
      -0.07026991,-0.05846931, -0.01016998,-0.03047472,-0.06280511},
      { 0.02513756,-0.00578333, 0.09865453, 0.01322314,-0.00310665,
        0.05880899,-0.09252443,-0.02986539,-0.03127460, 0.01007539,
       -0.00360119,-0.01995024, 0.00094940,-0.00145868,-0.01388816,
        0.11358341,-0.12127513,-0.00054696,-0.00055627, 0.00417284},
      { 0.16517316,-0.00254742,-0.03318745,-0.01984173, 0.00031890,
       -0.02817810, 0.02661678,-0.01761215, 0.01665112, 0.10513343,
       -0.00545026, 0.01827470,-0.00207616,-0.00763758,-0.01322808,
       -0.02202576,-0.07434204, 0.00020593, 0.00119979,-0.10827873},
      { 0.16088826, 0.00056313,-0.02579303,-0.00319655, 0.00037228,
       -0.03193150, 0.01655305,-0.03028640, 0.01367746,-0.11248153,
        0.00778371, 0.02675579, 0.00243718, 0.00895470,-0.01729803,
       -0.02686964,-0.08262584, 0.00011794,-0.00225134, 0.09415650},
      {-0.01739295, 0.00572017,-0.00712592,-0.01100922,-0.00870113,
       -0.00663461,-0.01153857,-0.02248432,-0.00382264,-0.00358612,
       -0.00139345,-0.00971460,-0.00133312, 0.01927783,-0.01053838,
       -0.00911362,-0.01010908, 0.09417598, 0.01763850,-0.00955454},
      { 0.01728888, 0.01344211, 0.01200836, 0.01857259,-0.17088517,
        0.01457592, 0.01997839, 0.02844884, 0.00839403, 0.00196862,
        0.01391984, 0.03270465, 0.00347173,-0.01940984, 0.01233979,
        0.00542887, 0.01008836, 0.00126491,-0.02863042, 0.00449764},
      {-0.02881366,-0.02184155,-0.01566086,-0.02593764,-0.04050907,
       -0.01539603,-0.02576729,-0.05089606,-0.00597430, 0.02181643,
        0.09835597,-0.04040940, 0.00873512, 0.12139434,-0.02427882,
       -0.02945238,-0.01566867,-0.01606503, 0.09475319, 0.02238670},
      { 0.04080274,-0.02869626,-0.05191093,-0.08435843, 0.00021141,
        0.13043842, 0.00871530, 0.00496058,-0.02797641,-0.00636933,
        0.02243277, 0.03640362,-0.05735517, 0.00196918,-0.02218934,
       -0.00608972, 0.02872922, 0.00047619, 0.00151285, 0.00883489},
      {-0.02623824, 0.00331152, 0.03640692, 0.04260231,-0.00038223,
       -0.07480340,-0.01022492,-0.00426473, 0.01448116, 0.01456847,
        0.05786680, 0.03368691,-0.10126924,-0.00147454, 0.01275395,
        0.00017574,-0.01585206,-0.00015767,-0.00231848, 0.02310137},
      {-0.00846258,-0.01508106,-0.01967505,-0.02772004, 0.01248253,
       -0.01331243,-0.02569382,-0.04461524,-0.02207075, 0.04663443,
        0.19347923,-0.02745691, 0.02288515,-0.04883849,-0.01084597,
       -0.01947187,-0.00081675, 0.00516540,-0.07815919, 0.08035585},
      {-0.06553111, 0.09756831, 0.00524326,-0.00885098, 0.00756653,
        0.02783099,-0.00427042,-0.16680359, 0.03951331,-0.00490540,
        0.01719610, 0.15018204, 0.00882722,-0.00423197,-0.01919217,
       -0.02963619,-0.01831342,-0.00524338, 0.00011379,-0.02566864},
      {-0.07494341,-0.11348850, 0.00241343,-0.00803016, 0.00492438,
        0.00711909,-0.00829147, 0.05793337, 0.02734209, 0.02059759,
       -0.02770280, 0.14128338, 0.01532479, 0.00364307, 0.05968116,
       -0.06497960,-0.08113941, 0.00319445,-0.00104222, 0.03553497},
      { 0.05948223,-0.08959930, 0.03269977,-0.03272374,-0.00365667,
       -0.03423294,-0.06418925,-0.05902138, 0.05746317,-0.02580596,
        0.01259572, 0.05848832, 0.00672666, 0.00233355,-0.05145149,
        0.07348503, 0.11427955, 0.00142592,-0.01030651,-0.04862799},
      {-0.01606880, 0.05200845,-0.01212967,-0.06824429,-0.00234304,
        0.01094203,-0.07375538, 0.08808629, 0.12394822, 0.02231351,
       -0.03608265,-0.06978045,-0.00618360, 0.00274747,-0.01921876,
       -0.01541969,-0.02223856,-0.00107603,-0.01251777, 0.05412534},
      { 0.01688843, 0.05784728,-0.02256966,-0.07072251,-0.00422551,
       -0.06261233,-0.08502830, 0.08925346,-0.08529597, 0.01519343,
       -0.05008258, 0.10931873, 0.00521033, 0.02593305,-0.00717855,
        0.02291527, 0.02527388,-0.00266188,-0.00871160, 0.02708135},
      {-0.04233344, 0.00076379, 0.01571257, 0.04003092, 0.00901468,
        0.00670577, 0.03459487, 0.12420216,-0.00067366,-0.01515094,
        0.05306642, 0.04338407, 0.00511287, 0.01036639,-0.17867462,
       -0.02289440,-0.03213205, 0.00017924,-0.01187362,-0.03933874},
      { 0.01284817,-0.01685622, 0.00724363, 0.01687952,-0.00882070,
       -0.00555957, 0.01676246,-0.05560456,-0.00966893, 0.06197684,
       -0.09058758, 0.00880607, 0.00108629,-0.08308956,-0.08056832,
       -0.00413297, 0.02973107, 0.00092948, 0.07010111, 0.13007418},
      { 0.00700223,-0.01347574, 0.00691332, 0.03122905, 0.00310308,
        0.00946862, 0.03455040,-0.06712536,-0.00304506, 0.04267941,
       -0.10422292,-0.01127831,-0.00549798, 0.11680505,-0.03352701,
       -0.00084536, 0.01631369, 0.00095063,-0.09570217, 0.06480321}
    };
    static double disteigval[20];
    for (int n = 0; n < 20; n++) disteigval[n] = exp(distance * logpameigval[n]);
    for (int j = 0; j < 20; j++) {
      for (int i= 0; i < 20; i++) {
    double &v = m[i][j];
    v = 0;
    for (int n = 0; n < 20; n++) {
      v += pameigvec[i][n] * disteigval[n] * pameigvecinv[n][j];
    }
  }
    }
    return;
}

int pamLoadSequences(const char *seqfile, const char *distfile, int *retseqnr, int *retseqlen, char **refseq) {
  //refseq needs to be uninitialized
  /*Errors:
  0 none
  1 no distance file
  2 no sequence file
  3 invalid number of sequences
  4 invalid sequence length
  5 reading error in distance file
  6 reading error in sequence file
  */
  int n;

  const int buflen = 50000;
  FILE *f = fopen(distfile, "r");
  if (f == NULL) return 1;
  char *buf; buf = new char[buflen];
  fgets(buf, buflen, f);
  seqnr = atoi(buf);
  if (seqnr < 1 || seqnr > 10000) return 3;
  dis = new Distance[seqnr];
  if (feof(f)) return 5;
  fgets(buf, buflen, f);

  for (n = 0; n < seqnr; n++) {
    dis[n].seq = n;
    float val = atof(buf + 15 + 9 * n);
    if (val < 0 || val > 10) return 5;
    dis[n].dist = val;
    AssemblePam(100 * val, dis[n].mat);
    for (int i = 0; i < 20; i++) {
      dis[n].expect[i] = 0;
      for (int ii = 0; ii < 20; ii++) {
        dis[n].expect[i] += dis[n].mat[i][ii] * dis[n].mat[i][ii];
      }
    }
  }
  fclose(f);

  #if DEBUG
  cout << "***Dis***" << endl;
  for (n = 0; n < seqnr; n++) {
    cout << fixed << setprecision(6) << dis[n].dist << endl;
  }
  cout << "******" << endl;
  #endif

  f = fopen(seqfile, "r");
  if (f == NULL) return 2;
  seq = new int*[seqnr];
  fgets(buf, buflen, f);
  sscanf(buf,"%*d %d", &seqlen);
  if (seqlen < 1 || seqnr > 10000) return 4;
  for (n = 0; n < seqnr; n++) {
    seq[n] = new int[seqlen];
    fgets(buf, buflen, f);
    if (n == 0) {
      *refseq = new char[seqlen];
      memcpy(*refseq, buf + 10, seqlen);
    }
    for (int nn = 0; nn < seqlen; nn++) {
      seq[n][nn] = code[buf[10 + nn]];
    }
  }

  qsort(dis, seqnr, sizeof(Distance), DistCompare);

  #if DEBUG
  cout << "***SortedDis***" << endl;
  for (n = 0; n < seqnr; n++) {
    cout << fixed << setprecision(6) << dis[n].dist << endl;
  }
  cout << "******" << endl;
  #endif

  seqtodis = new int[seqnr];
  for (n = 0; n < seqnr; n++) {
    seqtodis[dis[n].seq] = n;
  }

  #if DEBUG
  cout << "***Seqtodis***" << endl;
  for (n = 0; n < seqnr; n++) {
    cout << seqtodis[n] << endl;
  }
  cout << "******" << endl;
  #endif

  fclose(f);
  delete[] buf;
  *retseqnr = seqnr;
  *retseqlen = seqlen;
  return 0;
}

int pamGetAllDistances(float *distances) {
  for (int n = 1; n < seqnr; n++) {
    distances[n] = dis[n].dist;
  }
  return seqnr;
}

int pamCalcRawScores(int pos, float *scores) {
//INPUT pos; OUTPUT scores
//always returns seqnr-1 values; writes -1000 for '-'
  int n;

  int currnr;

  int vref = seq[0][pos];
  for (currnr = 1; currnr < seqnr; currnr++) {
    if (seq[dis[currnr].seq][pos] == -1) {
      scores[currnr-1] = -1000;
      continue;
    }

    matrix &m = dis[currnr].mat;
    //Use a different formula, by calculating expect values
    int vcomp = seq[dis[currnr].seq][pos];
    double weight;
    if (n == seqnr -1) weight = dis[currnr].dist - dis[currnr-1].dist;
    else weight = .5 * (dis[currnr+1].dist - dis[currnr-1].dist);
    double sim = 2.4 * (m[vref][vcomp] - dis[currnr+1].expect[vref]);
    scores[currnr-1] = sim * weight;
  };
  return currnr-1;
}

int pamCalcSimilarity(int pos, float *distances, float *scores) {
  //INPUT pos; OUTPUT distances, scores
  int n;

  double lastdist;
  int currnr;
  double currdist = 0;
  double nextdist;
  int nextnr;

  double sim;
  double totsim = 0;
  double weight;
  double totweight = 0;

  for (n = 1; n < seqnr; n++) {
    if (seq[dis[n].seq][pos] >= 0) {
      nextnr = n;
      nextdist = dis[n].dist;
      break;
    }
  }
  if (n == seqnr) return 0;


  weight = .5 * nextdist; totweight += weight; sim = 0; totsim = 0;

  int vref = seq[0][pos];
  int counter = 0;
  do {
    lastdist = currdist;
    currnr = nextnr;
    currdist = nextdist;
    for (n = currnr + 1; n < seqnr; n++) {
      if (seq[dis[n].seq][pos] >= 0) {
        nextnr = n;
        nextdist = dis[n].dist;
        break;
      }
    }
    if (n == seqnr) break;
    if (currdist == lastdist) continue;
    matrix &m = dis[currnr].mat;
    int vcomp = seq[dis[currnr].seq][pos];
    weight = .5 * (nextdist - lastdist);
    sim = 2.4 * (m[vref][vcomp] - dis[currnr].expect[vref]);
    //This scaling factor of 2.4 is totally arbitrary, but gives a nice range of scores. Scaling does not affect the final ranking of scores whatsoever
    totsim += weight * sim;
    distances[counter] = currdist;
    scores[counter] = totsim;
    counter++;
  } while (1);
  return counter;
}

int pamCalcSimilarityForRate(double rate, int pos, float *scores) {
extern int **seq;
extern int seqnr;
extern int seqlen;

  int n;

  double lastdist;
  int currnr;
  double currdist = 0;
  double nextdist;
  int nextnr;

  double sim;
  double totsim = 0;
  double weight;
  double totweight = 0;

  static matrix m;

  for (n = 1; n < seqnr; n++) {
    if (seq[dis[n].seq][pos] >= 0) {
      nextnr = n;
      nextdist = dis[n].dist;
      break;
    }
  }
  if (n == seqnr) return 0;


  weight = .5 * nextdist; totweight += weight; sim = 0; totsim = 0;

  int vref = seq[0][pos];
  int counter = 0;
  do {
    lastdist = currdist;
    currnr = nextnr;
    currdist = nextdist;
    for (n = currnr + 1; n < seqnr; n++) {
      if (seq[dis[n].seq][pos] >= 0) {
        nextnr = n;
        nextdist = dis[n].dist;
        break;
      }
    }
    if (n == seqnr) break;
    if (currdist == lastdist) continue;
    AssemblePam(100 * rate*currdist, m);
    double expect = 0;
    for (int i = 0; i < 20; i++) {
      expect += m[vref][i] * m[vref][i];
    }
    int vcomp = seq[dis[currnr].seq][pos];
    weight = .5 * (nextdist - lastdist);
    sim = 2.4 * (m[vref][vcomp] - expect);
    totsim += weight * sim;
    scores[counter] = totsim/*/totweight*/;
    counter++;
  } while (1);
  return counter;
}


void TryRate(double rate, int pos, float &score, bool &ascending, double d, float *scores) {
    int n;
    int posnr = pamCalcSimilarityForRate(rate, pos, scores);
    score = scores[posnr -1];
    posnr = pamCalcSimilarityForRate(rate + d, pos, scores);
    float score0 = scores[posnr -1];
    if (score0 > score) ascending = 1; else ascending = 0;
    return;
}

double pamCalcRate(int pos, double d) {
  float *scores = new float[seqnr];
  double minrate = 0, maxrate = 100;

  float score; bool ascending;
  TryRate(minrate, pos, score, ascending, d, scores);
  if (!ascending || score >= 0) {
    delete[] scores;
    return minrate;
  }
  TryRate(maxrate, pos, score, ascending, d, scores);
  if (ascending && score < 0) {
     delete[] scores;
    return maxrate;
  }

  double rate;
  while (1) {
    rate = .5 * (minrate + maxrate);
    if (rate - minrate < d) break;
    TryRate(rate, pos, score, ascending, d, scores);
    if (ascending && score < 0) minrate = rate;
    else maxrate = rate;
  }
  delete [] scores;
  return rate;
}
double pamCalcRate(int pos) {
  const double d = .001;
  return pamCalcRate(pos, d);
}

int pamCalcSimilarity(int pos, const int *sequencenrs, int nrsequencenrs, float *distances, float *scores) {
//For this version you have to provide a list of sequence numbers to use
  int n;

  double lastdist;
  int currnr;
  double currdist = 0;
  double nextdist;
  int nextnr;

  double sim;
  double totsim = 0;
  double weight;
  double totweight = 0;

  static matrix m;

  int nrdistancenrs = nrsequencenrs + 1;
  int *distancenrs = new int[nrdistancenrs];
  distancenrs[0] = 0;
  for (n = 0; n < nrsequencenrs; n++) {
    distancenrs[n+1] = seqtodis[sequencenrs[n]];
  }
  sort(distancenrs, distancenrs+nrdistancenrs);

  for (n = 1; n < nrdistancenrs; n++) {
    if (seq[dis[distancenrs[n]].seq][pos] >= 0) {
      nextnr = n;
      nextdist = dis[distancenrs[n]].dist;
      break;
    }
  }
  if (n == nrdistancenrs) {
    delete[] distancenrs;
    return 0;
  }


  weight = .5 * nextdist; totweight += weight; sim = 0; totsim = 0;

  int vref = seq[0][pos];
  int counter = 0;
  do {
    lastdist = currdist;
    currnr = nextnr;
    currdist = nextdist;
    for (n = currnr + 1; n < nrdistancenrs; n++) {
      if (seq[dis[distancenrs[n]].seq][pos] >= 0) {
        nextnr = n;
        nextdist = dis[distancenrs[n]].dist;
        break;
      }
    }
    if (n == nrdistancenrs) break;
    if (currdist == lastdist) continue;
    matrix &m = dis[distancenrs[currnr]].mat;
    double expect = 0;
    for (int i = 0; i < 20; i++) {
      expect += m[vref][i] * m[vref][i];
    }
    int vcomp = seq[dis[distancenrs[currnr]].seq][pos];
    weight = .5 * (nextdist - lastdist);
    sim = 2.4 * (m[vref][vcomp] - expect);
    totsim += weight * sim;
    distances[counter] = currdist;
    scores[counter] = totsim/*/totweight*/;
    counter++;
  } while (1);
  return counter;
}

void pamAdjustExpectedRate(double rate) {
  //Adjusts the null conservation matrices and expected values, which are initially for a rate 1
    for (int n = 0; n < seqnr; n++) {
      AssemblePam(100 * rate * dis[n].dist, dis[n].mat);
      for (int i = 0; i < 20; i++) {
        dis[n].expect[i] = 0;
        for (int ii = 0; ii < 20; ii++) {
          dis[n].expect[i] += dis[n].mat[i][ii] * dis[n].mat[i][ii];
        }
      }
    }
}
