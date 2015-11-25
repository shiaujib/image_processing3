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


Mat matSrc,matDst;
void init(String str){
	matSrc=imread(str,1);
	matDst=Mat(matSrc.size(),CV_64FC2);
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
	imshow("Fourier result",channels[1]);
	waitKey(0);
}

int main(){
	init("test.tif");
	//fourier_trans();
/*			for(int i=0;i<matSrc.cols;i++)
				for(int j=0;j<matSrc.rows;j++){
					//matDst.at<Vec2d>(j,i)[0]=0;
				}*/
	fourier_trans();
	//imshow("Fourier result",matDst[0]);
	waitKey(0);











}

