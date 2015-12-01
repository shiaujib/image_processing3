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
	for(int i=0;i<matSrc.cols;i++)
		for(int j=0;j<matSrc.rows;j++){
			matDst.at<Vec2d>(j,i)[0]=matSrc.at<uchar>(j,i);
			cout<<matDst.at<Vec2d>(j,i)[1]<<"\t";
		}
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
void fft_2d(Mat matin,bool dir){
	int cols,rows;
	cols=matin.cols;
	rows=matin.rows;
	double realPart[100000],imagePart[100000];
	
	for(int i=0;i<rows;i++){
		for(int j=0;j<cols;j++){
			realPart[j]=matin.at<Vec2d>(i,j)[0];
			imagePart[j]=matin.at<Vec2b>(i,j)[1];
		}
				
	

void fft(int dir,int bits,int nn,double *x,double *y){
	int l1,l2=1,i1;
	double u1=1,u2=0,cof1=-1,cof2=0,t1,t2,value;
	for(int i=0;i<bits;i++){
		l1=l2;
		l2/=2;
		for(int i=0;i<l1;i++){
			for(int j=i;j<nn;j+=12){
				i1=j+l1;
				t1=u1*x[i1]-u2*y[i1];
				t2=u1*y[i1]+u2*x[i1];
				x[i1]=x[j]-t1;
				y[i1]=y[j]-t2;
				x[j]+=t1;
				y[j]+=t2;
			}
			value=u1*cof1-u2*cof2;
			u2=u1*cof2+u2*cof1;
			u1=value;
		}
		cof2=sqrt((1-cof1)/2);
		if(dir==1)
			cof2*=-1;
		cof1=sqrt((1+cof1)/2);


	}
	if(dir==1){
		for(int i=0;i<nn;i++){
			x[i]/=nn;
			y[i]/=nn;
		}
	}
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

void zero_padding(int m,int n){
	/*int num=5,index=0,bindex=0;
	int tmp;
	int a[10]={1,2,3,4,5};
	int b[10];
	//matSrc.cols=num;
	for(int i=0;i<10;i++){
		if(!num%2)
			tmp=bindex/2-1;
		else
			tmp=(bindex-1)/2;	
		if(index>tmp){
			num++;
			b[i]=0;
			bindex++;
		}
		else{
			b[i]=a[index++];
			bindex++;
		}
	}
	cout<<"num "<<num;
	for(int i=0;i<num;i++)
		cout<<"\t "<<b[i];*/
	int a[5]={1,2,3,4,5};
	int *b=(int *)malloc(10*sizeof(int));
	for(int i =0;i<m;i++){
		if(i<n)
			b[i]=a[i];
		else
			b[i]=0;
	}
	
	for(int i =0;i<m;i++)
		cout<<b[i]<<" \t";
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
//	bit_reverse();
//	zero_padding(10,5);

	//imshow("Fourier result",matDst[0]);
//	waitKey(0);











}

