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


Mat matSrc,matDst,matZero,matResult,matFFT;
void init(String str){
	matSrc=imread(str,0);
	matDst=Mat(matSrc.size(),CV_64FC2);
	matResult=Mat(matSrc.size(),CV_64FC2);
	matFFT=Mat(matSrc.size(),CV_64FC1);
	for(int i=0;i<matSrc.cols;i++)
		for(int j=0;j<matSrc.rows;j++){
			matDst.at<Vec2d>(j,i)[0]=matSrc.at<uchar>(j,i);
//			cout<<matDst.at<Vec2d>(j,i)[0]<<"----\t";
		}
	cout<<matSrc.cols<<endl;
	cout<<matSrc.rows<<endl;
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

void swap(int *x,int *y){
	double tmp;
	tmp=*x;
	*x=*y;
	*y=tmp;
}

void bit_reverse(double *real,double *image,int *bits,int *nn){
	int i,j,k,i2;
	double tx,ty;
	
	*nn=1;
   	for (i=0;i<*bits;i++)
      		(*nn) *= 2;

   	i2 = *nn >> 1;
   	j = 0;
   	for (i=0;i<*(nn)-1;i++) {
      		if (i < j) {
			swap(real[i],real[j]);
			swap(image[i],image[j]);
      		}
      	k = i2;
      	while (k <= j) {
        j -= k;
        k >>= 1;
        }
     	j += k;
     	//cout<<j<<endl;
        }

}



bool zero_padding(int length,int *bits,int *value){
	if(length<1){
		*bits=0;
		*value=1;
		return false;
	}
	*bits=1;
	*value=2;
	do{
		(*bits)++;
		(*value)*=2;
	}while(2*(*value)<=length);
	if((*value)!=length)
		return false;
	else
		return true;

}


void fft(int dir,int *bits,int *nn,double *real,double *image){
	int l1,l2=1,i1;
	double u1,u2,cof1=-1,cof2=0,t1,t2,value;
	for(int i=0;i<*bits;i++){
		u1=1;
		u2=0;
		l1=l2;
		l2<<=1;
		for(int i=0;i<l1;i++){
			for(int j=i;j<*nn;j+=l2){
				i1=j+l1;
				t1=u1*real[i1]-u2*image[i1];
				t2=u1*image[i1]+u2*real[i1];
				real[i1]=real[j]-t1;
				image[i1]=image[j]-t2;
				real[j]+=t1;
				image[j]+=t2;
			}
			value=u1*cof1-u2*cof2;
			u2=u1*cof2+u2*cof1;
			u1=value;
		}
		cof2=sqrt((1.0-cof1)/2.0);
		if(dir==1)
			cof2*=-1;
		cof1=sqrt((1.0+cof1)/2.0);


	}
	if(dir==1){
		for(int i=0;i<*nn;i++){
			real[i]/=(double)(*nn);
			image[i]/=(double)(*nn);
		}
	}
}
		
	
	

int  fft_2d(Mat matin,Mat matout,int dir){
	int cols,rows;
	int bits,value;
	int nn;
	double scale;
	nn=1;
	cols=matin.cols;
	rows=matin.rows;
	//double realPart[100000],imagePart[100000];
	double *realPart=(double *)malloc(cols*sizeof(double));
	double *imagePart=(double *)malloc(cols*sizeof(double));
	
	for(int i=0;i<cols;i++)
		for(int j=0;j<rows;j++)
			matin.at<Vec2d>(j,i)[0]=matin.at<Vec2d>(j,i)[0]*pow(-1,i+j);
	if(!zero_padding(cols,&bits,&value)){
		cout<<"image length is not power of two"<<endl;
		return 0;
	}
	
	for(int i=0;i<rows;i++){
		for(int j=0;j<cols;j++){
			realPart[j]=matin.at<Vec2d>(i,j)[0];
			imagePart[j]=matin.at<Vec2d>(i,j)[1];
		}
		
		bit_reverse(realPart,imagePart,&bits,&nn);
		fft(1,&bits,&nn,realPart,imagePart);
		for(int k=0;k<cols;k++){
			matout.at<Vec2d>(i,k)[0]=realPart[k];
			matout.at<Vec2d>(i,k)[1]=imagePart[k];
			
		}
	}
	free(realPart);
	free(imagePart);
	realPart=(double *)malloc(rows*sizeof(double));
	imagePart=(double *)malloc(rows*sizeof(double));
	if(!zero_padding(rows,&bits,&value)){
		cout<<"image length is not power of two"<<endl;
		return 0;
	}
	
	for(int i=0;i<cols;i++){
		for(int j=0;j<rows;j++){
			realPart[j]=matout.at<Vec2d>(j,i)[0];
			imagePart[j]=matout.at<Vec2d>(j,i)[1];
		}
		
		bit_reverse(realPart,imagePart,&bits,&nn);
		fft(1,&bits,&nn,realPart,imagePart);
		for(int k=0;k<rows;k++){
			matout.at<Vec2d>(k,i)[0]=realPart[k];
			matout.at<Vec2d>(k,i)[1]=imagePart[k];
			
		}
	}
	free(realPart);
	free(imagePart);
	for(int i=0;i<cols;i++)
		for(int j=0;j<rows;j++){
		//	if(matout.at<Vec2d>(j,i)[0])
		//		cout<<matout.at<Vec2d>(j,i)[1]<<"----"<<i<<"  "<<j<<endl;
			matFFT.at<Vec2d>(j,i)=sqrt(pow(matout.at<Vec2d>(j,i)[0],2)+pow(matout.at<Vec2d>(j,i)[1],2));
	//		matFFT.at<uchar>(j,i)=matout.at<Vec2d>(j,i)[1];
	}
	imshow("FFT result",matFFT);
	waitKey(0);
		
	
}
		


/*void zero_padding(int m,int n){
	int num=5,index=0,bindex=0;
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
		cout<<"\t "<<b[i];
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
}*/




int main(){
	init("figure1.tif");
	fft_2d(matDst,matResult,0);
	for(int i=0;i<matSrc.cols;i++){
		for(int j=0;j<matSrc.rows;j++){
		//	cout<<matResult.at<Vec2d>(j,i)[0]<<"\t"<<i<<"  "<<j<<endl;
		}
	//	cout<<endl;
	}
//	fourier_trans();
			for(int i=0;i<matSrc.cols;i++)
				for(int j=0;j<matSrc.rows;j++){
					matDst.at<Vec2d>(j,i)[1]=0;
				}
//	fourier_trans();
//	fast_fourier();
//	bit_reverse();
//	zero_padding(10,5);

	//imshow("Fourier result",matDst[0]);
//	waitKey(0);











}

