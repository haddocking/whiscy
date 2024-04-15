BEGIN {
  ok = 1;
  if (chainid == "") {
    print "USAGE: awk -f getsurface.awk -v chainid=<something> [cutoff=<something>]";
    ok = 0;
  }
  if (int(cutoff)<=0) cutoff = 15;
}

{ 
  chain = substr($0,9,1)
  nr = substr($0,10,4)
  reltot =  substr($0,23,6) + 0
  relside = substr($0,36,6) + 0
  relmain = substr($0,49,6) + 0
  #print $0
  #print chain, nr, reltot, relside, relmain
}
ok == 1 && $1 == "RES"  && chain == chainid {
  #g = 0; if (reltot > cutoff) g = 1
  #print "Test", reltot, cutoff, g
  if (reltot >= cutoff || relside >= cutoff || relmain >= cutoff) {
  
      r = 0;
      for (n = length(nr); n > 0; n--) {
        r = int(substr(nr,1,n));
        if (r > 0) break;
      }
      if (r > 0 && !(r in printed)) {
        printed[r] = 1;
	print r
      }
  }
}
