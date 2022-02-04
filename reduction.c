#include <stdlib.h>

#include <gmp.h>

/* See references (and a more readable code) in: bezier-reduction.jl */

static int min(int a, int b) {
  if (a <= b)
    return a;
  return b;
}

static int max(int a, int b) {
  if (a >= b)
    return a;
  return b;
}

static void transpose(int n, mpq_t *A) {
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < i; ++j) {
      int index1 = i * n + j;
      int index2 = j * n + i;
      mpq_swap(A[index1], A[index2]);
    }
}

static void add(int m, int n, mpq_t *A, mpq_t *B) {
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j) {
      int index = i * n + j;
      mpq_add(A[index], A[index], B[index]);
    }
}

static mpq_t *multiply(int m, int r, int n, mpq_t *A, mpq_t *B) {
  mpq_t *C = (mpq_t *)malloc(m * n * sizeof(mpq_t));
  mpq_t tmp;
  mpq_init(tmp);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j) {
      int index = i * n + j;
      mpq_init(C[index]);
      for (int k = 0; k < r; ++k) {
        int index_a = i * r + k, index_b = k * n + j;
        mpq_mul(tmp, A[index_a], B[index_b]);
        mpq_add(C[index], C[index], tmp);
      }
    }
  mpq_clear(tmp);
  return C;
}

static void clear(int size, mpq_t *A) {
  for (int i = 0; i < size; ++i)
    mpq_clear(A[i]);
  free(A);
}

static mpq_t *computeL(int n, int r, int s) {
  mpq_t *L = (mpq_t *)malloc((n + 1) * (n + 1) * sizeof(mpq_t));
  mpz_t b1z, b2z, b3z;
  mpq_t b1, b2, b3, tmp;
  mpq_t b, d;
  mpz_inits(b1z, b2z, b3z, NULL);
  mpq_inits(b1, b2, b3, tmp, b, d, NULL);
  for (int k = 0; k <= n; ++k)
    for (int j = 0; j <= n; ++j) {
      int index = k * (n + 1) + j;
      mpq_init(L[index]);
      for (int i = max(0, j + k - n); i <= min(j, k); ++i) {
        mpz_bin_uiui(b1z, k, i);
        mpz_bin_uiui(b2z, n - k, j - i);
        mpz_bin_uiui(b3z, n, j);
        mpq_set_z(b1, b1z);
        mpq_set_z(b2, b2z);
        mpq_set_z(b3, b3z);
        mpq_mul(b, b1, b2);
        mpq_div(b, b, b3);
        int sign = (k + i) % 2 == 0 ? 1 : -1;
        mpz_bin_uiui(b1z, k + r, i);
        mpz_bin_uiui(b2z, k + s, k - i);
        mpz_bin_uiui(b3z, k, i);
        mpq_set_z(b1, b1z);
        mpq_set_z(b2, b2z);
        mpq_set_z(b3, b3z);
        mpq_mul(d, b1, b2);
        mpq_div(d, d, b3);
        mpq_set_d(tmp, (double)sign);
        mpq_mul(tmp, tmp, b);
        mpq_mul(tmp, tmp, d);
        mpq_add(L[index], L[index], tmp);
      }
    }
  mpz_clears(b1z, b2z, b3z, NULL);
  mpq_clears(b1, b2, b3, tmp, b, d, NULL);
  return L;
}

static mpq_t *computeE(int n, int r, int s) {
  mpq_t *E = (mpq_t *)malloc((n + 1) * (n + 1) * sizeof(mpq_t));
  mpz_t b1z, b2z, b3z;
  mpq_t b1, b2, b3, tmp;
  mpq_t d;
  mpz_inits(b1z, b2z, b3z, NULL);
  mpq_inits(b1, b2, b3, tmp, d, NULL);
  for (int k = 0; k <= n; ++k)
    for (int j = 0; j <= n; ++j) {
      int index = k * (n + 1) + j;
      mpq_init(E[index]);
      for (int i = 0; i <= j; ++i) {
        int sign = (j + i) % 2 == 0 ? 1 : -1;
        mpz_bin_uiui(b1z, j + r, i);
        mpz_bin_uiui(b2z, j + s, j - i);
        mpz_bin_uiui(b3z, n + r + s + j, k + s + i);
        mpq_set_z(b1, b1z);
        mpq_set_z(b2, b2z);
        mpq_set_z(b3, b3z);
        mpq_mul(d, b1, b2);
        mpq_div(d, d, b3);
        mpq_set_d(tmp, (double)sign);
        mpq_mul(tmp, tmp, d);
        mpq_add(E[index], E[index], tmp);
      }
      mpz_bin_uiui(b1z, j + r + s, r);
      mpz_bin_uiui(b2z, n, k);
      mpz_bin_uiui(b3z, j + r, r);
      mpq_set_z(b1, b1z);
      mpq_set_z(b2, b2z);
      mpq_set_z(b3, b3z);
      mpq_mul(d, b1, b2);
      mpq_div(d, d, b3);
      mpq_set_d(tmp, (double)(2 * j + r + s + 1));
      mpq_mul(E[index], E[index], tmp);
      mpq_set_d(tmp, (double)(n + r + s + j + 1));
      mpq_div(E[index], E[index], tmp);
      mpq_mul(E[index], E[index], d);
    }
  mpz_clears(b1z, b2z, b3z, NULL);
  mpq_clears(b1, b2, b3, tmp, d, NULL);
  return E;
}

static mpq_t *computeA(int n, int m, int l) {
  mpq_t *A = (mpq_t *)malloc((l + 1) * sizeof(mpq_t));
  mpz_t bz;
  mpq_t b, tmp;
  mpz_init(bz);
  mpq_inits(A[0], b, tmp, NULL);
  mpq_set_d(A[0], 1.0);
  for (int k = 1; k <= l; ++k) {
    mpq_init(A[k]);
    for (int i = 0; i < k; ++i) {
      mpz_bin_uiui(bz, n - m, k - i);
      mpq_set_z(b, bz);
      mpq_mul(tmp, b, A[i]);
      mpq_sub(A[k], A[k], tmp);
    }
  }
  mpz_clear(bz);
  mpq_clears(b, tmp, NULL);
  return A;
}

static mpq_t *computeD(int N, int n, int r) {
  mpq_t *D = (mpq_t *)malloc((N + 1) * (N + 1) * sizeof(mpq_t));
  mpz_t b1z, b2z;
  mpq_t b1, b2;
  mpz_inits(b1z, b2z, NULL);
  mpq_inits(b1, b2, NULL);
  for (int i = 0; i <= N; ++i)
    for (int j = 0; j <= N; ++j) {
      int index = i * (N + 1) + j;
      mpq_init(D[index]);
      if (i != j)
        continue;
      mpz_bin_uiui(b1z, n, r + 1 + i);
      mpz_bin_uiui(b2z, N, i);
      mpq_set_z(b1, b1z);
      mpq_set_z(b2, b2z);
      mpq_div(D[index], b1, b2);
    }
  mpz_clears(b1z, b2z, NULL);
  mpq_clears(b1, b2, NULL);
  return D;
}

static mpq_t *computeC(int N, int n, int m, int r, int s, mpq_t *A) {
  mpq_t *C = (mpq_t *)malloc((N + 1) * (n + 1) * sizeof(mpq_t));
  mpz_t b1z, b2z;
  mpq_t b1, b2, d, tmp;
  mpz_inits(b1z, b2z, NULL);
  mpq_inits(b1, b2, d, tmp, NULL);
  for (int j = 0; j <= N; ++j)
    for (int k = 0; k <= n; ++k) {
      int index = j * (n + 1) + k;
      mpq_init(C[index]);
    }
  for (int j = r + 1; j <= n + r - m; ++j)
    for (int k = 0; k <= r; ++k) {
      int index = (j - r - 1) * (n + 1) + k;
      mpz_bin_uiui(b1z, n, k);
      mpz_bin_uiui(b2z, n, j);
      mpq_set_z(b1, b1z);
      mpq_set_z(b2, b2z);
      mpq_div(d, b1, b2);
      mpq_set_d(tmp, 0.0);
      for (int i = max(k, j - (n - m)); i <= r; ++i) {
        mpz_bin_uiui(b1z, n - m, j - i);
        mpq_set_z(b1, b1z);
        mpq_mul(b1, b1, A[i-k]);
        mpq_add(tmp, tmp, b1);
      }
      mpq_mul(tmp, tmp, d);
      mpq_sub(C[index], C[index], tmp);
    }
  for (int j = r + 1; j < n - s; ++j) {
      int index = (j - r - 1) * (n + 1) + j;
      mpq_set_d(C[index], 1.0);
  }
  for (int j = m - s; j < n - s; ++j)
    for (int k = n - s; k <= n; ++k) {
      int index = (j - r - 1) * (n + 1) + k;
      mpz_bin_uiui(b1z, n, k);
      mpz_bin_uiui(b2z, n, j);
      mpq_set_z(b1, b1z);
      mpq_set_z(b2, b2z);
      mpq_div(d, b1, b2);
      mpq_set_d(tmp, 0.0);
      for (int i = max(n - k, m - j); i <= s; ++i) {
        mpz_bin_uiui(b1z, n - m, j - m + i);
        mpq_set_z(b1, b1z);
        mpq_mul(b1, b1, A[i+k-n]);
        mpq_add(tmp, tmp, b1);
      }
      mpq_mul(tmp, tmp, d);
      mpq_sub(C[index], C[index], tmp);
    }
  mpz_clears(b1z, b2z, NULL);
  mpq_clears(b1, b2, d, tmp, NULL);
  return C;
}

static mpq_t *computeD1(int m, int M, int r, int s) {
  mpq_t *D1 = (mpq_t *)malloc((m + 1) * (M + 1) * sizeof(mpq_t));
  mpz_t b1z, b2z;
  mpq_t b1, b2;
  mpz_inits(b1z, b2z, NULL);
  mpq_inits(b1, b2, NULL);
  for (int i = 0; i <= m; ++i)
    for (int j = 0; j <= M; ++j) {
      int index = i * (M + 1) + j;
      mpq_init(D1[index]);
      if (i <= r || i >= m - s || j != i - r - 1)
        continue;
      mpz_bin_uiui(b1z, M, i - r - 1);
      mpz_bin_uiui(b2z, m, i);
      mpq_set_z(b1, b1z);
      mpq_set_z(b2, b2z);
      mpq_div(D1[index], b1, b2);
    }
  mpz_clears(b1z, b2z, NULL);
  mpq_clears(b1, b2, NULL);
  return D1;
}

static mpq_t *computeQ1(int m, int n, int r, int s, mpq_t *A) {
  mpq_t *Q1 = (mpq_t *)malloc((m + 1) * (n + 1) * sizeof(mpq_t));
  mpz_t b1z, b2z;
  mpq_t b1, b2;
  mpz_inits(b1z, b2z, NULL);
  mpq_inits(b1, b2, NULL);
  for (int j = 0; j <= m; ++j)
    for (int k = 0; k <= n; ++k) {
      int index = j * (n + 1) + k;
      mpq_init(Q1[index]);
      if (j <= r && k <= j) {
        mpz_bin_uiui(b1z, n, k);
        mpz_bin_uiui(b2z, m, j);
        mpq_set_z(b1, b1z);
        mpq_set_z(b2, b2z);
        mpq_div(Q1[index], b1, b2);
        mpq_mul(Q1[index], Q1[index], A[j-k]);
      } else if (j >= m - s && k >= j + (n - m)) {
        mpz_bin_uiui(b1z, n, k);
        mpz_bin_uiui(b2z, m, j);
        mpq_set_z(b1, b1z);
        mpq_set_z(b2, b2z);
        mpq_div(Q1[index], b1, b2);
        mpq_mul(Q1[index], Q1[index], A[k-j-(n-m)]);
      }
    }
  mpz_clears(b1z, b2z, NULL);
  mpq_clears(b1, b2, NULL);
  return Q1;
}

static mpq_t *computeQ2(int m, int n, int M, int N,
                        mpq_t *D1, mpq_t *L, mpq_t *E, mpq_t *D, mpq_t *C) {
  mpq_t *I = (mpq_t *)malloc((M + 1) * (N + 1) * sizeof(mpq_t));
  for (int i = 0; i <= M; ++i)
    for (int j = 0; j <= N; ++j) {
      int index = i * (N + 1) + j;
      mpq_init(I[index]);
      if (i == j)
        mpq_set_d(I[index], 1.0);
    }
  transpose(M + 1, L);
  transpose(N + 1, E);
  mpq_t *tmp1 = multiply(m + 1, M + 1, M + 1, D1, L);
  mpq_t *tmp2 = multiply(m + 1, M + 1, N + 1, tmp1, I);
  clear((m + 1) * (M + 1), tmp1);
  mpq_t *tmp3 = multiply(m + 1, N + 1, N + 1, tmp2, E);
  clear((m + 1) * (N + 1), tmp2);
  mpq_t *tmp4 = multiply(m + 1, N + 1, N + 1, tmp3, D);
  clear((m + 1) * (N + 1), tmp3);
  mpq_t *Q2 = multiply(m + 1, N + 1, n + 1, tmp4, C);
  clear((m + 1) * (N + 1), tmp4);
  return Q2;
}

void reduction_matrix(int n, int m, int r, int s, double *Q) {
  int N = n - (r + s + 2), M = m - (r + s + 2);
  mpq_t *L = computeL(M, 2 * r + 2, 2 * s + 2);
  mpq_t *E = computeE(N, 2 * r + 2, 2 * s + 2);
  mpq_t *A = computeA(n, m, max(r, s));
  mpq_t *D = computeD(N, n, r);
  mpq_t *C = computeC(N, n, m, r, s, A);
  mpq_t *D1 = computeD1(m, M, r, s);
  mpq_t *Q1 = computeQ1(m, n, r, s, A);
  mpq_t *Q2 = computeQ2(m, n, M, N, D1, L, E, D, C);
  add(m + 1, n + 1, Q1, Q2);
  for (int i = 0; i < (m + 1) * (n + 1); ++i)
    Q[i] = mpq_get_d(Q1[i]);
  clear((m + 1) * (n + 1), Q2);
  clear((m + 1) * (n + 1), Q1);
  clear((m + 1) * (M + 1), D1);
  clear((N + 1) * (n + 1), C);
  clear((N + 1) * (N + 1), D);
  clear(max(r, s) + 1, A);
  clear((N + 1) * (N + 1), E);
  clear((M + 1) * (M + 1), L);
}
