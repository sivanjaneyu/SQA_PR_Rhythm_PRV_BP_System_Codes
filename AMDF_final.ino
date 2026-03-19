#include <LiquidCrystal.h>
const int rs = 12, en = 11, d4 = 5, d5 = 4, d6 = 3, d7 = 2;
LiquidCrystal lcd(rs, en, d4, d5, d6, d7);

unsigned long t1,t2,t3,t4;
const int p=6; //length of a and b (order of filter+1)
float b[6]={0.864246524017246,-4.32123262008623,8.64246524017246,-8.64246524017246,4.32123262008623,-0.864246524017246};
float a[6]={1,-4.70835055262605,8.87551850311413,-8.37268603379321,3.95242919071962,-0.746904488298868};
const int Fs=125;
const int t=5;
const int M=Fs*t;
int i,PR,SQ=1;

float x1[M]={10464,9545,8364,6921,5412,3707,1936,296,-1279,-2919,-4559,-6133,-7642,-9020,-10332,-11512,-12496,-13415,-14071,-14530,-14858,-15120,-15317,-15448,-15580,-15645,-15645,-15711,-15711,-15776,-15908,-15973,-15973,-16039,-16039,-16170,-16301,-16498,-16629,-16826,-17023,-17220,-17416,-17679,-17876,-18072,-18335,-18663,-18925,-19187,-19450,-19712,-19975,-20237,-20499,-20631,-20565,-20303,-19778,-18859,-17548,-15776,-13677,-11250,-8560,-5740,-2853,-32,2657,4822,6593,8364,10004,11448,12628,13547,14072,14268,14137,13744,13088,12300,11316,10136,8758,7249,5609,3904,2198,427,-1213,-2853,-4296,-5674,-6855,-7773,-8692,-9544,-10266,-10922,-11447,-11840,-12103,-12234,-12300,-12300,-12300,-12234,-12103,-11906,-11840,-11709,-11644,-11578,-11578,-11512,-11512,-11512,-11512,-11512,-11512,-11512,-11512,-11578,-11644,-11775,-11972,-12234,-12562,-12890,-13218,-13480,-13743,-14005,-14202,-14268,-14202,-13808,-13087,-11972,-10397,-8429,-6133,-3640,-1410,755,3116,5609,8036,10332,12300,14072,15515,16630,17417,17876,18007,17876,17548,17089,16302,15318,14072,12628,10988,9348,7577,5806,4100,2526,1148,-360,-1869,-3312,-4624,-5805,-6855,-7773,-8495,-9085,-9610,-9938,-10200,-10463,-10725,-10856,-11053,-11250,-11381,-11447,-11512,-11119,-10463,-10135,-10004,-9938,-9938,-10004,-10135,-10266,-10397,-10528,-10725,-10922,-11119,-11381,-11578,-11775,-11972,-12168,-12496,-12824,-13152,-13480,-13743,-13874,-13743,-13349,-12628,-11775,-10660,-9216,-7380,-5346,-3181,-951,1345,3576,5609,7446,9086,10332,11251,12104,12825,13284,13481,13481,13219,12694,11907,10923,9742,8496,7249,6003,4560,3116,1608,99,-1344,-2722,-3968,-5149,-6199,-7183,-7904,-8495,-9020,-9544,-9938,-10200,-10463,-10660,-10791,-10922,-10988,-10988,-10988,-10988,-10988,-10988,-10988,-11053,-11184,-11316,-11512,-11644,-11840,-11972,-12103,-12168,-12234,-12300,-12496,-12693,-12890,-13152,-13415,-13743,-14071,-14333,-14596,-14727,-14792,-14858,-14727,-14464,-13874,-12956,-11775,-10200,-8364,-6330,-4100,-1869,33,1739,3641,5412,7118,8561,9873,10857,11644,12169,12366,12300,11972,11382,10660,9742,8692,7577,6265,4888,3379,1870,296,-1279,-2788,-4034,-5149,-6199,-7248,-8167,-9020,-9676,-10266,-10725,-11119,-11447,-11644,-11775,-11840,-11840,-11906,-11906,-11906,-11906,-11972,-11972,-12037,-12037,-12037,-12103,-12103,-12168,-12168,-12234,-12365,-12496,-12628,-12824,-13021,-13218,-13480,-13677,-13940,-14071,-14333,-14596,-14792,-14989,-15186,-15252,-15186,-14858,-14333,-13480,-12234,-10988,-9610,-7904,-5871,-3772,-1541,624,2657,4625,6396,7905,9217,10136,10726,11185,11513,11513,11316,10857,10201,9348,8233,6987,5675,4232,2920,1608,164,-1279,-2656,-4034,-5412,-6658,-7773,-8823,-9741,-10528,-11119,-11578,-11972,-12234,-12496,-12759,-12890,-13021,-13152,-13218,-13218,-13284,-13349,-13415,-13415,-13480,-13546,-13612,-13677,-13743,-13874,-14005,-14136,-14268,-14464,-14596,-14727,-14924,-15186,-15973,-17023,-17679,-18138,-18531,-18859,-19122,-19384,-19647,-19712,-19515,-19187,-18466,-17482,-16170,-14399,-12365,-10069,-7642,-5084,-2525,-426,1411,3313,5150,6790,8102,9217,10004,10464,10660,10464,10004,9283,8364,7315,6068,4560,2985,1280,-426,-2263,-3968,-5674,-7314,-8888,-10069,-11119,-12103,-13021,-13808,-14530,-15055,-15514,-15842,-16039,-16170,-16170,-16236,-16170,-16104,-16039,-16039,-16039,-16039,-16039,-16039,-16039,-16039,-16170,-16301,-16367,-16432,-16564,-16695,-16826,-17023,-17154,-17351,-17482,-17679,-17941,-18203,-18335,-18466,-18597,-18728,-18859,-18859,-18728,-18400,-17744,-16826,-15514,-13808,-11709,-9676,-7642,-5215,-2591,99,2657,5150,7446,9414,11120,12497,13547,14137,14465,14531,14465,14137,13547,12694,11644,10398,8955,7380,5675,3904,2395,1017,-360,-1804,-3181,-4493,-5674,-6724,-7642,-8298,-8888,-9348,-9676,-9872,-10004,-10069,-10069,-10004,-9938,-9872,-9741,-9676,-9544,-9544,-9544,-9544,-9610,-9676,-9741,-9872,-10004,-10200,-10332,-10463,-10660,-10791,-10988,-11119,-11250,-11512,-11775,-12103,-12431,-12759,-13021,-13284,-13480,-13612,-13612,-13480,-13152,-12562,-11512,-10069,-8298,-6133,-3706
}; //original signal
//float x1[M];  //real-time sensing
float x3[M]; //Normalised signal
float x4[M];  //Filtered signal
const int k=500; // length of partial AMDF
float th=-0.25; //Threshold for threshold crossing of derivative
int ZC[M];
int AMDF_zc[k]; //Threshold crossings of AMDF_ds
int num,N1,N2,Z;     //Number of th crossings
float AMDF_val[k];     //AMDF function
float AMDF_d[k];   //Derivative of AMDF
float AMDF_neg[k];   //Negative of AMDF
float AMDF_ds[k];  //Smoothed derivative of AMDF
float der[k],der1[k],der2[k];



void setup() {
  Serial.begin(9600);
  Serial1.begin(9600);
  
  lcd.begin(16, 4);
  // Print a message to the LCD.
  lcd.setCursor(0,1);
  lcd.print("Start!!");
  delay(1000);
}

void loop() {
  t1=micros();

// //Real-time sensing
//    i=0;
//    while(i<samples)
//    {
//      x1[i]=analogRead(A0);
//      delay(8);
//      i=i+1;
//    }


   normalize_amplitude(x1, M, x3);
   filtfilt(x3, M, a, b, p, x4);

  // Step 3: Derivative filtering and Zero Crossing Rate (ZCR)
    float xd[M - 1], xs[M - 1];
    derivative(x4, M, xd);
    MAF(xd, M - 1, 3, xs); // equivalent to rectwin(3)/3 in MATLAB

    int NPZCR = msm_zerocros(xs, 0.0, M - 1, ZC); // All zero crossings

   t2=micros();
   AMDF(x4, M, k, AMDF_val);
   t3=micros();
   mean_removal(AMDF_val, k, AMDF_val);
   normalize_amplitude(AMDF_val, k, AMDF_val);
   for (i=0; i<k; i++){
    AMDF_neg[i]=-AMDF_val[i];
   }
   derivative(AMDF_neg, k, AMDF_ds);
   //MAF(AMDF_d, k, 5, AMDF_ds);
   k_lag(AMDF_ds, k, 5, AMDF_ds);
   normalize_amplitude(AMDF_ds, k, AMDF_ds);
   num=neg_th_cross(AMDF_ds, th, k, AMDF_zc);
   minima_correction(AMDF_val, k, AMDF_zc, num, AMDF_zc);
   PR=ceil(Fs*60/(AMDF_zc[0]+1));
   t4=micros();

   N1=neg_th_cross(x4, 0.1, M, ZC);
   N2=neg_th_cross(x4, -0.1, M, ZC);
   Z = abs(N1-N2);


  if (NPZCR > 8 && NPZCR <= 38) SQ = 0;
  else if (NPZCR > 38 && NPZCR < 168) SQ = 1;
  else if (num < 2 || PR < 30 || PR > 300) SQ = 1;
  else if ((AMDF_val[AMDF_zc[0]] > -0.13 && AMDF_val[AMDF_zc[1]] > -0.35) && mag > 64) SQ = 1;
  else if ((AMDF_val[AMDF_zc[0]] < -0.37 && AMDF_val[AMDF_zc[1]] < -0.34) && mag <= 9) SQ = 0;
  else {
    for (int i = 0; i < num - 2; i++) {
      int d1 = AMDF_zc[i + 1] - AMDF_zc[i];
      int d2 = AMDF_zc[i + 2] - AMDF_zc[i + 1];
      int dist = abs(d2 - d1);
      if (dist > 8) {
        SQ = 1;
        break;
      } else {
        SQ = 0;
      }
    }
  }

//  for (int i=0; i< k;i++){
//    Serial.print(AMDF_val[i]-5);
//    Serial.print(',');
//    Serial.println(AMDF_ds[i]);
//  } 
 Serial.print("PR is");
 Serial.print(',');
Serial.println(PR);
Serial.print("Signal Quality is");
 Serial.print(',');
Serial.println(SQ);

// for (int i=0; i< 100;i++){
//    Serial.println(AMDF[i]-5);
//    Serial.print(',');
//    Serial.println(AMDF_der[i]+5);
//    Serial.print(',');
//   Serial.println(demo[i]);
//  }


//LCD display
lcd.setCursor(0,1);
lcd.print("Pulse rate:");
lcd.setCursor(12,1);
lcd.print(PR);      // print value on lcd display  
delay(1000);
lcd.setCursor(0,1);
lcd.print("Signal Quality:");
lcd.setCursor(12,1);
if (SQ==0) lcd.print("Good"); 
else  lcd.print("Bad"); 
delay(1000);
lcd.setCursor(0,1);
lcd.print("Time elapsed:");
lcd.setCursor(12,1);
lcd.print(t4-t1);      // print time on lcd display  
Serial1.print("Pulse rate:");
Serial1.println(PR);
Serial1.print(',');
Serial1.print("Time elapsed:");
Serial1.println(t4-t1);
Serial1.print(',');
Serial1.print("Signal Quality:");
if (SQ==0) Serial1.println("Good"); 
else  Serial1.println("Bad"); 
}



void normalize_amplitude(float x[], int N, float y[]){
      
//      Function to normalize the input signal to [-1,1] range
//      Input: signal x with length N
//      Output: Normalised signal y of length N
      
      int i;
      float max_x=x[0];

// Finding Max value
      for (int i=0; i<N; i++)
      {
        max_x=max(max_x, abs(x[i]));
      }
      
// Amplitude normalisation
      for(int i=0; i<N; i++)
      {
        y[i]=x[i]/max_x;  
      }
}


void AMDF(float x[], int N, int M, float y[]){

//  Function to find the partial AMDF of the given input signal
//  Input: signal x with length N
//  Output: AMDF y calculated for the signal upto a lag of M (length of AMDF)

  int i,j;
  float sum;
  
  for(int i=0; i<M; i++)
  {
       sum=0;
       for(int j=0; j<(N-i-1); j++)
       {
         sum=sum+abs(x[j]-x[j+i]);   
       }
       y[i]=sum/(N-i-1);
  }
}


void mean_removal(float x[], int N, float y[]){

//      Function for mean_removal of the input signal
//      Input: signal x with length N
//      Output: Zero-mean signal y of length N
      
      int i;
      float avg,sum=0;

// Mean Finding
      for (int i=0; i<N; i++)
      {
        sum=sum+x[i];
      }
      avg=sum/N;
      
// Mean removal
      for(int i=0; i<N; i++)
      {
        y[i]=x[i]-avg;  
      }
}


void k_lag(float x[], int N, int p, float y[]){

// Function to set the initial p values to (p+1)th
// Input: signal x with length N
// Output: signal y of length N with initial p values set to (p+1)th

  int i;
  
  for (int i=0; i<p; i++){
    y[i]=y[p];
  }

  for (int i=p; i<N; i++){
    y[i]=x[i];
  }
}

int neg_th_cross(float x[], float th, int N, int y[]){

//  Function to find the negative threshold-crossing points of the given input signal
//  Input: signal x with length N
//  Output: Array of 'th'-crossing locations along the entire length of the signal
//  Returns the number of the 'th-crossings' 

  int i,j=0;
  
  for(int i=0; i<N; i++)
      {  
        if (x[i]>th && x[i+1]<=th){
          y[j]=i;
          j++;
        }
      }

      if (j==0){
        y[j]=N-1;
        return 1;
      }
      return j;
}


void minima_correction(float x[], int N, int loc[], int M, int y[]){

//  Function to correct the minima value by finding the local minima in a small window of 8 samples
//  Input: signal x of length N whose M minima locations are given by array loc
//  Output: M corrected minima locations
  
  int i,j;
  float minima;
  
  for(int i=0; i<M;i++){
    minima=x[loc[i]];
    y[i]=loc[i];
    for (int j=loc[i]-3; j<loc[i]+4; j++){
      if (x[j]<minima){
          minima=x[j];
          y[i]=j;
        }
    }
  }    
}


void filtfilt(float x[], int N, float a[], float b[], int p, float y[]){

//  Function for filtfilt.
//  Implements 2 linear convolution: 1 in forward direction and again on the reversed signal. The signal is then reversed again to obtain filtfilt
//  Zero-phase forward and reverse digital IIR filtering.
//      After filtering in the forward direction, the filtered data is
//      reversed and run back through the filter; Y is the time reverse of
//      the output of the second filtering operation. The result has
//      precisely zero phase distortion, and magnitude modified by the
//      square of the filter's magnitude response. Startup and ending
//      transients are minimized by matching initial conditions.
//  Input: signal to be filtered of length N and filter coefficients a and b of length p.
//  NB: 1)Pad zeroes to the filter coefficients a and b and make it same length.
//      2)p is not the order of the filter. p can be order+1 for low and high pass filters whereas it can be 2*order+1 for band pass and bandstop filters
//  Output: filtered signal y of length N

  int i,j;
  float val;
  float z[N];
  
  filter(x, N, a, b, p, z);
  for (int i=0; i<N; i++){
    y[i]=z[N-i-1];  //Reversing
  }
  filter(y, N, a, b, p, z);
  for (int i=0; i<N; i++){
    y[i]=z[N-i-1];  //Reversing
  }
}


void filter(float x[], int N, float a[], float b[], int p, float y[]){
  
//  Function for applying linear filter in 1 direction
//  a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb) - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
//  Input: signal to be filtered of length N and filter coefficients a and b of length p.
//  NB: 1)Pad zeroes to the filter coefficients a and b and make it same length.
//      2)p is not the order of the filter. p can be order+1 for low and high pass filters whereas it can be 2*order+1 for band pass and bandstop filters
//  Output: filtered signal y of length N

  int i,j;
  float val;
  
  for (int i=0; i<N; i++){
    j=1;
    val=0;
    while(j<p && i>=j){
      val=val+b[j]*x[i-j]-a[j]*y[i-j];
      j++;
    }
    y[i]=b[0]*x[i]+ val;
  }
}



void MAF(float x[], int N, int p, float y[]){

//Function to apply a p-th order moving average filter on the given signal using filtfilt
//The filter coefficients for filtfilt are internally generated
//Input: signal to be filtered of length N and filter order p.
//Output: filtered signal y of length N

  int i;
  float a[p],b[p];

  b[0]=1/p;
  a[0]=1;
  for (int i=1; i<p; i++){
    b[i]=1/p;
    a[i]=0;
  }
  filtfilt(x, N, a, b, p, y);
}


void derivative(float x[], int N, float y[]){
  
//  Function to find the 1o' derivative of the given input signal y[n]=x[n+1]-x[n]
//  Input: signal to be differentiated of length N.
//  Output: 1 o' derivative y of signal x of length N

  int i;
  
  for(int i=0; i<N-1; i++){
    y[i]=x[i+1]-x[i];
  }
  y[N-1]=y[N-2];
}


int msm_zerocros(float* signal, int len, int* locations) {
  int count = 0;
  for (int i = 1; i < len; i++) {
    if (signal[i] * signal[i - 1] <= 0) {
      locations[count] = i;
      count++;
    }
  }
  return count;
}
