import os
from math import exp, isclose
from libwhiscy.pam_data import logpameigval, pameigvec, pameigvecinv, code


class Distance():
    def __init__(self, seq=0, dist=0., mat=None, expect=None):
        self.seq = seq
        self.dist = dist
        if mat is None:
            self.mat = [[0 for x in range(20)] for y in range(20)]
        else:
            self.mat = mat
        if expect is None:
            self.expect = [0. for x in range(20)]
        else:
            self.expect = expect


def get_pam_assemble(distance):
    m = [[0. for x in range(20)] for y in range(20)]
    disteigval = [0.0 for x in range(20)]

    for n in range(20):
        disteigval[n] = exp(distance * logpameigval[n])

    for j in range(20):
        for i in range(20):
            v = 0
            for n in range(20):
                v += pameigvec[i][n] * disteigval[n] * pameigvecinv[n][j]
            m[i][j] = v
    return m


def pam_load_sequences(alignment_file, distance_file):
    """Loads the aligment and distance files"""
    if not os.path.exists(distance_file):
        raise SystemExit("ERROR: Distance file {0} does not exist".format(distance_file))

    if not os.path.exists(alignment_file):
        raise SystemExit("ERROR: Sequence file {0} does not exist".format(alignment_file))

    refseq = ''
    seqtodis = []

    seqnr = 0
    distances = []
    with open(distance_file, "rU") as input_distances:
        first_line = input_distances.readline().rstrip(os.linesep)
        fields = first_line.split()
        seqnr = int(fields[0])
        if seqnr < 1 or seqnr > 10000:
            raise SystemExit("ERROR: Invalid number of sequences")

        raw = input_distances.readline().rstrip(os.linesep).split()

        if len(raw) != (seqnr + 1):
            raise SystemExit("ERROR: Reading error in distance file {0}".format(distance_file))

        for n in range(seqnr):
            seq = n
            dist = 0
            try:
                val = float(raw[1 + n])
                if val < 0 or val > 10:
                    raise ValueError()
                dist = val
            except ValueError:
                raise SystemExit("ERROR: Reading error in distance file {0}".format(distance_file))

            m = get_pam_assemble(100 * dist)

            expect = [0 for i in range(20)]
            for i in range(20):
                for ii in range(20):
                    expect[i] += m[i][ii] * m[i][ii]

            d = Distance(seq, dist, m, expect)
            distances.append(d)
    
    seqlen = 0
    sequences = [[] for _ in range(seqnr)]
    with open(alignment_file, "rU") as input_alignment:
        first_line = input_alignment.readline().rstrip(os.linesep)
        fields = first_line.split()
        seqlen = int(fields[1])
        if seqlen < 1 or seqnr > 10000:
            raise SystemExit("ERROR: Invalid sequence length")

        for n in range(seqnr):
            sequences[n] = []
            line = input_alignment.readline().rstrip(os.linesep)
            if line:
                name, sequence = line[:10], line[10:]
                if n == 0:
                    refseq = sequence
                for c in sequence:
                    sequences[n].append(code[ord(c)])

    # Sorted in ascending order as the C++ qsort
    sorted_distances = sorted(distances, key=lambda distance: distance.dist, reverse=False)
    seqtodis = [0 for _ in range(seqnr)]
    for n in range(seqnr):
        seqtodis[sorted_distances[n].seq] = n

    # print("***Dis***")
    # for x in distances:
    #     print("%.6f" % x.dist)
    #     print(' '.join([("%.6f" % i)  for i in x.expect]) + ' ')
    # print("******")

    # print("***SortedDis***")
    # for x in sorted_distances:
    #     print("%.6f" % x.dist)
    # print("******")

    # print("***Seqtodis***")
    # for x in seqtodis:
    #     print(x)
    # print("******")

    return seqnr, seqlen, refseq, sorted_distances, sequences, seqtodis


def pam_calc_similarity(pos, seqnr, seq, dis):
    nextnr = 0
    currnr = 0
    nextdist = 0.
    currdist = 0.
    lastdist = 0.
    scores = [0. for _ in range(seqnr)]
    distances = [0. for _ in range(seqnr)]
    for n in range(1, seqnr):
        if seq[dis[n].seq][pos] >= 0:
            nextnr = n
            nextdist = dis[n].dist
            break
    if n == seqnr: 
        return 0, distances, scores

    sim = 0.
    totsim = 0.
    weight = 0.5 * nextdist
    totweight = weight

    vref = seq[0][pos]
    counter = 0

    while True:
        lastdist = currdist
        currnr = nextnr
        currdist = nextdist
        for n in range(currnr + 1, seqnr):
            if seq[dis[n].seq][pos] >= 0:
                nextnr = n
                nextdist = dis[n].dist
                break
        if n == (seqnr - 1):
            break
        if isclose(currdist, lastdist): 
            continue
        
        m = dis[currnr].mat
        vcomp = seq[dis[currnr].seq][pos]
        weight = .5 * (nextdist - lastdist)
        # This scaling factor of 2.4 is totally arbitrary, but gives a nice range of scores. 
        # Scaling does not affect the final ranking of scores whatsoever
        sim = 2.4 * (m[vref][vcomp] - dis[currnr].expect[vref])
        
        totsim += weight * sim
        distances[counter] = currdist
        scores[counter] = totsim
        counter += 1

    return counter, distances, scores

