#include<opencv2/core/core.hpp>
#include<opencv2/highgui/highgui.hpp>
#include "opencv2/opencv.hpp"
#include<iostream>
#include<math.h>
#include<complex>
#include<math.h>
#define pi 3.14159
using namespace std;
using namespace cv;


Mat matSrc,matDst,matZero;
void init(String str){
	matSrc=imread(str,0);
	matDst=Mat(matSrc.size(),CV_64FC2);
	/*for(int i=0;i<matSrc.cols;i++)
		for(int j=0;j<matSrc.rows;j++){
			matDst.at<Vec2d>(j,i)[0]=matSrc.at<uchar>(j,i);
			cout<<matDst.at<Vec2d>(j,i)<<endl;
			cout<<matSrc.at<uchar>(j,i)<<endl;
		}*/
	matZero=Mat(matSrc.size(),CV_64FC1,0);
//	merge(matSrc,matZero,matDst);
	
	imshow("input",matSrc);
	waitKey(0);
	
}

void fourier_trans(){
	double xvalue=0,yvalue=0,value=0,sum_0=0,sum_1=0;
	int cols,rows;
	vector<Mat> channels(3);
	cols=matSrc.cols;
	rows=matSrc.rows;
	for(int u=0;u<matSrc.cols;u++)
		for(int v=0;v<matSrc.rows;v++){
			for(int i=0;i<matSrc.cols;i++)
				for(int j=0;j<matSrc.rows;j++){
					xvalue=2*pi*u*i/cols;
					yvalue=2*pi*v*j/rows;
					sum_0+=matSrc.at<uchar>(j,i)*cos(xvalue)*cos(yvalue)-sin(xvalue)*sin(yvalue);
					sum_1+=matSrc.at<uchar>(j,i)*-2*sin(xvalue)*sin(yvalue);
					//cout<<matDst.at<Vec2d>(v,u)<<endl;
				}
		matDst.at<Vec2d>(v,u)[0]=sum_0;
		matDst.at<Vec2d>(v,u)[1]=sum_1;
		cout<<matDst.at<Vec2d>(v,u)<<endl;
		sum_0=0;
		sum_1=0;

	}
	split(matDst,channels);
	imshow("Fourier result",channels[0]);
	waitKey('0');
	imshow("Fourier result",channels[1]);
	waitKey(0);
}

void fast_fourier(){
	int j=8;
	j=j>>1;
	cout<<j;


}
void swap(int *x,int *y){
	double tmp;
	tmp=*x;
	*x=*y;
	*y=tmp;
}
void bit_reverse(/*int *x,int *y*/){
	int nn,i,i1,j,k,i2;
	double tx,ty;
	int x[8]={1,2,3,4,5,6,7,8};
	int y[8]={8,7,6,5,4,3,2,1};
	int m=3;

  	nn = 1;
   	for (i=0;i<m;i++)
      		nn *= 2;

   	i2 = nn >> 1;
   	j = 0;
   	for (i=0;i<nn-1;i++) {
      		if (i < j) {
			swap(x[i],x[j]);
			swap(y[i],y[j]);
      		}
      	k = i2;
      	while (k <= j) {
        j -= k;
        k >>= 1;
        }
     	j += k;
     	cout<<j<<endl;
        }

}

int main(){
	init("test.tif");
	//fourier_trans();
/*			for(int i=0;i<matSrc.cols;i++)
				for(int j=0;j<matSrc.rows;j++){
					//matDst.at<Vec2d>(j,i)[0]=0;
				}*/
//	fourier_trans();
//	fast_fourier();
	bit_reverse();
	//imshow("Fourier result",matDst[0]);
//	waitKey(0);











}

