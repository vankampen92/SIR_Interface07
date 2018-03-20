#include <math.h>
 
#define JMAX 20 /* JMAX 20 */
#define JMAXP (JMAX+1)
#define K 5

float qromb_Accuracy(float (*func)(float), float a, float b, float EPS)
{
        /* Numerical Recipes Originaln Value EPS = 1.0e-6 */
	void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
	float trapzd(float (*func)(float), float a, float b, int n);
	void nrerror(char error_text[]);
	float ss,dss;
	float s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb_Accuracy");
	return 0.0;
}
#undef JMAX
#undef JMAXP
#undef K






