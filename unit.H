#ifndef __UNIT_H
#define __UNIT_H

#include "dcomp.H"

/*============================== FUNCTIONS ===========================*/
// Use these for making a list of atoms and sorting it.

typedef struct 
{
  int R[3];
  double Rmagn;
} unit_point_type;

int __COMPARE__UNIT__ (const void* a, const void* b) 
{
  double amag = ((unit_point_type*)a)->Rmagn;
  double bmag = ((unit_point_type*)b)->Rmagn;
  if (dcomp(amag,bmag)) return 0;
  if (amag<bmag) return -1;
  return 1;
}
#endif
