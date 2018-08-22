## our attempt to write a statistical accumulator.
BEGIN{firstfile = 0; nrec=0; nlines=0;}

# see if we're still dealing with the first file...
(! firstfile){if (NR != FNR) firstfile=1;}

# else, deal with it line by line, file by file
{
  if (! firstfile) {
    ++nlines;
    tag[FNR] = $1;
    if ( ($1+0) == $1) {
      data[FNR]=1;
      if ((NF-1)>nrec) nrec=(NF-1);
    }
    else data[FNR]=0;
  }
  if (data[FNR]) {
# there's data to accumulate
    for (j=0; j<nrec; ++j) {
      sum[FNR,j] += $(j+2);
      sum2[FNR,j] += $(j+2) * $(j+2);
    }
  }
}

END{nfiles = NR/FNR;
 for(i=1; i<=nlines; ++i) {
   if (! data[i])
     print tag[i];
   else {
     printf("%.8le", tag[i]);
     for (j=0; j<nrec; ++j) {
       ave = sum[i,j] / nfiles;
       dev = sum2[i,j]/nfiles - ave*ave;
       if (dev>0) dev = sqrt(dev); else dev=0;
       printf(" %.8le %.8le", ave, dev);
     }
     printf("\n");
   }
 }
}
