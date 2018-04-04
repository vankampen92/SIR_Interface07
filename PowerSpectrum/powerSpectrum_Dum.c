/* Simultaneous calculation of the power spectra from two time series */
void powerSpectrum_K_Sets(float *x, float *y, int nP, int k, int m)
{
  float *px, *py;
  int i;

  if(k == 1){

    powerSpectrum(x,y,nP);

  }
  else{
    
    if(nP != 4*k*m){
      printf("The number of temporal points points is not correct!!\n");
      printf("Time Series Length is %d, but should be %d\n", nP, 4*k*m);
      exit(0);
    }
    
    px = vector(1,m+1); py = vector(1,m+1);
 
    spectral_Density(x, px, m, k, 0);
    spectral_Density(y, py, m, k, 0);
    
    for(i=0; i<=m; i++) {
      x[i] = px[i];
      y[i] = py[i];
    }

    free_vector(px, 1,m+1); free_vector(py, 1,m+1);
  }  
}

void powerSpectrum(float *x, float *y, int nP)
{
  /* Power Spectra of two series x and y at at time */
  int i,j;
  unsigned long n;
  float factor;
  float *fft1, *fft2;
  
  fft1 = vector(1,2*nP); fft2 = vector(1, 2*nP);
  n = nP;

  twofft(x,y,fft1,fft2,n);
  
  j=0;
  for(i=1; i<2*n; i+=2){
    j++;
    x[j] = fft1[i]*fft1[i] + fft1[i+1]*fft1[i+1];
    y[j] = fft2[i]*fft2[i] + fft2[i+1]*fft2[i+1];
  } 
 
  for(i=1; i<n/2; i++){
    x[i+1] = x[i+1] + x[n-i+1];
    y[i+1] = y[i+1] + y[n-i+1];
  }

  free_vector(fft1, 1,2*nP); free_vector(fft2, 1,2*nP);
}
