#pragma once

/*
  Arguments:
  - n: original degree
  - m: target degree
  - r: # of derivatives to retain at u = 0.0
  - s: # of derivatives to retain at u = 1.0
  - Q: an array of (m + 1) * (n + 1) doubles

  Q will be filled in row-major order:
    (0,0), (0,1), ..., (0,n), (1,0), (1,1), ..., (m,n)
 */
void reduction_matrix(int n, int m, int r, int s, double *Q);
