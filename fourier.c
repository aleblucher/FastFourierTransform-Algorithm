#include <math.h>
#include <stdio.h>

#include "fourier.h"

void nft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {
  for (int k = 0; k < n; k++) {
    t[k] = 0;

    for (int j = 0; j < n; j++) {
      t[k] += s[j] * cexp(sign * 2 * PI * k * j * I / n);
    }
  }
}

void nft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
  nft(s, t, n, -1);
}

void nft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
  nft(t, s, n, 1);

  for (int k = 0; k < n; k++) {
    s[k] /= n;
  }
}

void fft(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n, int sign) {
  // Base
  if (n == 1){
    t[0] = s[0];
    return;
  }
  
  // Recursão
  double complex s_p[n/2];
  double complex s_i[n/2];
  double complex t_p[n/2];
  double complex t_i[n/2];

  for (int i = 0; i < n/2; i++){
    s_p[i] = s[i*2];    //Criando vetor de índices pares
    s_i[i] = s[i*2+1];  //Criando vetor de índices ímpares
  } 
  
  fft(s_p,t_p,n/2,sign);
  fft(s_i,t_i,n/2,sign);

  // Passo
  for(int k = 0; k < n/2; k++){
    double complex cpx = cexp(sign * 2 * PI * k * I / n);
    t[k] = t_p[k] + t_i[k] * cpx;
    t[k + n/2] = t_p[k] - t_i[k] * cpx;
    //printf("%f", t[k]); 
  }

}

void fft_forward(double complex s[MAX_SIZE], double complex t[MAX_SIZE], int n) {
  fft(s, t, n, -1);
}

void fft_inverse(double complex t[MAX_SIZE], double complex s[MAX_SIZE], int n) {
  fft(t, s, n, 1);

  for (int k = 0; k < n; k++) {
    s[k] /= n;
  }
}

void fft_forward_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
  for (int i = 0; i <= height; i ++){
    for (int j = 0; j <= width; j ++){
      //double complex ind[width] = matrix[j];
      //fft_forward(ind, ind, 1);
    }
  }

  for (int j = 0; j <= height; j ++){
    for (int i = 0; i <= width; i ++){
      //fft_forward(matrix[j], matrix[j], 1);
    }
  }
}


void fft_inverse_2d(double complex matrix[MAX_SIZE][MAX_SIZE], int width, int height) {
  for (int i = 0; i <= width; i ++){
    for (int j = 0; j <= height; j ++){
      fft_inverse(matrix[i], matrix[i], 1);
    }
  }

  for (int j = 0; j <= height; j ++){
    for (int i = 0; i <= width; i ++){
      fft_inverse(matrix[j], matrix[j], 1);
    }
  }
}

void filter(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height, int flip) {
  int center_x = width / 2;
  int center_y = height / 2;

  double variance = -2 * SIGMA * SIGMA;

  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int dx = center_x - (x + center_x) % width;
      int dy = center_y - (y + center_y) % height;

      double d = dx * dx + dy * dy;

      double g = exp(d / variance);

      if (flip) {
        g = 1 - g;
      }

      output[y][x] = g * input[y][x];
    }
  }
}

void filter_lp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
  filter(input, output, width, height, 0);
}

void filter_hp(double complex input[MAX_SIZE][MAX_SIZE], double complex output[MAX_SIZE][MAX_SIZE], int width, int height) {
  filter(input, output, width, height, 1);
}
