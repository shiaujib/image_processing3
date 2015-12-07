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
//typedef float real;
//typedef struct{real Re; real Im;} complex;

Mat matSrc,matSrc2,matDst,matZero,matResult,matFFT,matIFFT,matH,matHP,matIHP,matTmp,matRifft,matRifft2,matHFEF;

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
	if((*value)!=length){
		(*bits)++;
		(*value)*=2;
		
		return false;
	}
	else
		return true;


}
void init(String str){
	int bits,value;
	int newcols,newrows;
	matSrc=imread(str,0);
	//matResult=Mat(matSrc.size(),CV_64FC2);
	//matFFT=Mat(matSrc.size(),CV_64FC1);
	//matIFFT=Mat(matSrc.size(),CV_64FC1);
	zero_padding(matSrc.cols,&bits,&value);
	newcols=value;
	zero_padding(matSrc.rows,&bits,&value);
	newrows=value;
	matSrc2=Mat::zeros(Size(newcols,newcols),0);
	matResult=Mat(matSrc2.size(),CV_64FC2);
	matDst=Mat(matSrc2.size(),CV_64FC2);
	matFFT=Mat(matSrc2.size(),CV_64FC1);
	matIFFT=Mat(matSrc2.size(),CV_64FC1);
	matH=Mat(matSrc2.size(),CV_64FC1);
	matHP=Mat(matSrc2.size(),CV_64FC2);
	matIHP=Mat(matSrc2.size(),CV_64FC1);
	matTmp=Mat(matSrc2.size(),CV_64FC1);
	matRifft=Mat(matSrc.size(),CV_64FC1);
	matRifft2=Mat(matSrc.size(),CV_64FC1);
	matHFEF=Mat(matSrc.size(),CV_64FC1);
	for(int i=0;i<matSrc.cols;i++)
		for(int j=0;j<matSrc.rows;j++){
			matSrc2.at<uchar>(j,i)=matSrc.at<uchar>(j,i);
//			cout<<matDst.at<Vec2d>(j,i)[0]<<"----\t";
		}
	//imshow("new src",matSrc2);
	waitKey(0);
	for(int i=0;i<matSrc2.cols;i++)
		for(int j=0;j<matSrc2.rows;j++){
			matDst.at<Vec2d>(j,i)[0]=matSrc2.at<uchar>(j,i);
//			cout<<matDst.at<Vec2d>(j,i)[0]<<"----\t";
		}
	cout<<matSrc2.cols<<endl;
	cout<<matSrc2.rows<<endl;
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





void fft(int dir,int *bits,int *nn,double *real,double *image){
	
	int l1,l2=1,i1;
	double u1,u2,cof1=-1.0,cof2=0.0,t1,t2,value;
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
		cof2=(double)sqrt((1-cof1)/2);
		if(dir==1)
			cof2*=(double)-1;
		cof1=(double)sqrt((1+cof1)/2);


	}
	if(dir==1){
		for(int i=0;i<*nn;i++){
			real[i]/=(double)(*nn);
			image[i]/=(double)(*nn);
		}
	}
}

/*void fft_2(double *real,double *image,int *nn){
	for(int i=0;i<*nn;i++){
		for(int k=0;k<0.5*(*nn);k+=2){
			real[k]*
*/		

void center_trans(Mat matin){
	for(int i=0;i<matin.cols;i++)
		for(int j=0;j<matin.rows;j++)
			matin.at<Vec2d>(j,i)[0]=matin.at<Vec2d>(j,i)[0]*pow(-1,i+j);
}
	
void magnitude(Mat matin,Mat matout){
	int cols,rows;
	cols=matin.cols;
	rows=matin.rows;
	for(int i=0;i<cols;i++){
		for(int j=0;j<rows;j++){
			matout.at<Vec2d>(j,i*0.5)=sqrt(pow(matin.at<Vec2d>(j,i)[0],2)+pow(matin.at<Vec2d>(j,i)[1],2));
		}
	}

	//imshow("magnitude",matout);
	//waitKey(0);
}





int  fft_2d(Mat matin,Mat matout,int dir){
	int cols,rows;
	int bits,value;
	int nn;
	double scale;
	nn=1;
	cols=matin.cols;
	rows=matin.rows;
	double *realPart=(double *)malloc(cols*sizeof(double));
	double *imagePart=(double *)malloc(cols*sizeof(double));
	/*for(int i=0;i<cols;i++)
		for(int j=0;j<rows;j++)
			matin.at<Vec2d>(j,i)[0]=matin.at<Vec2d>(j,i)[0]*pow(-1,i+j);
*/
	center_trans(matin);	
	if(!zero_padding(cols,&bits,&value)){
		cout<<"image length is not power of two  "<<value<<endl;
		/*matSrc=Mat::zeros(Size(value,matSrc.rows),matSrc.type());
		imshow("zero_padding",matSrc);
		waitKey(0);*/
		return 0;
	}
	
	for(int i=0;i<cols;i++){
		for(int j=0;j<rows;j++){
			realPart[j]=matin.at<Vec2d>(j,i)[0];
			imagePart[j]=matin.at<Vec2d>(j,i)[1];
		}
		nn=1;
		bit_reverse(realPart,imagePart,&bits,&nn);
		//cout<<nn<<"---------"<<endl;
		//cout<<bits<<"---------"<<endl;
		fft(1,&bits,&nn,realPart,imagePart);
		for(int k=0;k<rows;k++){
			matout.at<Vec2d>(k,i)[0]=realPart[k];
			matout.at<Vec2d>(k,i)[1]=imagePart[k];
			
		}
	}
	free(realPart);
	free(imagePart);
	realPart=(double *)malloc(rows*sizeof(double));
	imagePart=(double *)malloc(rows*sizeof(double));
	if(!zero_padding(rows,&bits,&value)){
		cout<<"image length is not power of two  "<<value<<endl;
		return 0;
	}
	
	for(int i=0;i<rows;i++){
		for(int j=0;j<cols;j++){
			realPart[j]=matout.at<Vec2d>(i,j)[0];
			imagePart[j]=matout.at<Vec2d>(i,j)[1];
		}
		
		bit_reverse(realPart,imagePart,&bits,&nn);
		//cout<<"nn"<<nn<<endl<<"bits"<<bits<<endl;
		fft(1,&bits,&nn,realPart,imagePart);
		for(int k=0;k<cols;k++){
			matout.at<Vec2d>(i,k)[0]=realPart[k];
			matout.at<Vec2d>(i,k)[1]=imagePart[k];
			
		}
	}
//	cout<<"444444444444444"<<endl;
	free(realPart);
	free(imagePart);

//	cout<<"444444444444444"<<endl;

/*	for(int i=0;i<cols;i++)
		for(int j=0;j<rows;j++){
			matFFT.at<Vec2d>(j,i*0.5)=sqrt(pow(matout.at<Vec2d>(j,i)[0],2)+pow(matout.at<Vec2d>(j,i)[1],2));
			
			
			
			//cout<<matFFT.at<Vec3b>(j,i)<<endl;;	
			//		matFFT.at<uchar>(j,i)=matout.at<Vec2d>(j,i)[1];
		//	matFFT.at<Vec2d>(j,i)=log(matFFT.at<Vec2d>(j,i));
	}*/
	magnitude(matout,matFFT);
	imshow("FFT result",matFFT);
//	imwrite("fft_result.tif",matFFT);
	waitKey(0);
		
	
}
	


int  ifft_2d(Mat matin,Mat matresult,int dir,String str,int flag){
	int cols,rows;
	int bits,value;
	int nn;
	double scale;
	double result;
	nn=1;
	cols=matin.cols;
	rows=matin.rows;
	double *realPart=(double *)malloc(cols*sizeof(double));
	double *imagePart=(double *)malloc(cols*sizeof(double));
	Mat matout;
	matout=Mat(matSrc2.size(),CV_64FC2);
/*	for(int i=0;i<cols;i++)
		for(int j=0;j<rows;j++)
			matin.at<Vec2d>(j,i)[0]=matin.at<Vec2d>(j,i)[0]*pow(-1,i+j);
*/
	zero_padding(cols,&bits,&nn);	
	for(int i=0;i<cols;i++){
		for(int j=0;j<rows;j++){
			realPart[j]=matin.at<Vec2d>(j,i)[0];
			imagePart[j]=matin.at<Vec2d>(j,i)[1];
		}
		nn=1;
		bit_reverse(realPart,imagePart,&bits,&nn);
		fft(-1,&bits,&nn,realPart,imagePart);
		for(int k=0;k<rows;k++){
			matout.at<Vec2d>(k,i)[0]=realPart[k];
			matout.at<Vec2d>(k,i)[1]=imagePart[k];
			
		}
	}
	
	//cout<<"6666666666"<<endl;
	free(realPart);
	free(imagePart);
	//cout<<"6666666666"<<endl;
	realPart=(double *)malloc(rows*sizeof(double));
	imagePart=(double *)malloc(rows*sizeof(double));
//	zero_padding(rows,&bits,&nn);
	for(int i=0;i<rows;i++){
		for(int j=0;j<cols;j++){
			realPart[j]=matout.at<Vec2d>(i,j)[0];
			imagePart[j]=matout.at<Vec2d>(i,j)[1];
		}
		
		bit_reverse(realPart,imagePart,&bits,&nn);
		fft(-1,&bits,&nn,realPart,imagePart);
		for(int k=0;k<cols;k++){
			matout.at<Vec2d>(i,k)[0]=realPart[k];
			matout.at<Vec2d>(i,k)[1]=imagePart[k];
			
		}
	}
	free(realPart);
	free(imagePart);
	//center_trans(matout);

	double max,min;
	for(int i=0;i<cols;i++)
		for(int j=0;j<rows;j++){
			result=sqrt(pow(matout.at<Vec2d>(j,i)[0],2)+pow(matout.at<Vec2d>(j,i)[1],2))/(nn/4);
/*			if(result>255)
				result=255;
			if(result<0)
				result=0;*/
			
			matresult.at<double>(j,i)=result;
		/*	if(result>max)
				max=result;
			if(result<min);
				min=result;*/
		//	cout<<result<<endl;
			
			
			
	}
	
	//center_trans(matIFFT);
	for(int i=0;i<matSrc.cols;i++)
		for(int j=0;j<matSrc.rows;j++){
			matRifft.at<Vec2d>(j,i)=matresult.at<Vec2d>(j,i+0.5*(cols-matSrc.cols));
			if(flag==1)
				matRifft2.at<Vec2d>(j,i)=matresult.at<Vec2d>(j,i+0.5*(cols-matSrc.cols));
	}
	
//	imshow(str,matresult);
//	waitKey(0);
	imshow(str,matRifft);
//	imwrite("fft_result.tif",matFFT);
	waitKey(0);
		
	
}
	


void highPassFilter(Mat matin,Mat matout,double D0){
	int cols,rows;
	cols=matin.cols;
	rows=matin.rows;
	double D,value;
	for(int i=0;i<cols;i++)
		for(int j=0;j<rows;j++){
			D=sqrt(pow(i,2)+pow(j,2));
			/*if(D>D0)
				value=1;
			else if(D<=D0)
				value=0;*/
			value=1-exp((-pow(D,2)/(2*pow(D0,2))));
			//value=1-exp((-pow(D,2)/(2*pow(D0,2))));
			matout.at<Vec2d>(j,i)[0]=value*matin.at<Vec2d>(j,i)[0];
			matout.at<Vec2d>(j,i)[1]=value*matin.at<Vec2d>(j,i)[1];
		}	
	magnitude(matout,matTmp);
	imshow("test",matTmp);	
	//imshow("Gaussian high pass",matout);
	waitKey(0);
}


void HFEF(Mat matin ,Mat matout,int k){
	int cols,rows;
	cols=matin.cols;
	rows=matin.rows;
	double D,value;
	Mat gmask,tmp;
	double max,min;
	gmask=Mat(matSrc.size(),CV_64FC1);
	tmp=Mat(matSrc.size(),CV_8UC1);
	for(int i=0;i<cols;i++)
		for(int j=0;j<rows;j++){
			gmask.at<double>(j,i)=matRifft2.at<double>(j,i)-matin.at<double>(j,i);
			matout.at<double>(j,i)=matRifft2.at<double>(j,i)+k*gmask.at<double>(j,i);
			if(matout.at<double>(j,i)>max)
				max=matout.at<double>(j,i);
			if(matout.at<double>(j,i)<min)
				min=matout.at<double>(j,i);
		}
	
	for(int i=0;i<cols;i++)                         //normalize
		for(int j=0;j<rows;j++)
			tmp.at<uchar>(j,i)=(matout.at<double>(j,i)-min)*255/(max-min);
			
	imshow("HFEF RESULT",tmp);
	waitKey(0);
	equalizeHist( tmp,tmp );
	imshow("histogram equalize",tmp);
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
//	init("Fig0516(a)(applo17_boulder_noisy).tif");
//	init("figure1.tif");
//	init("test.tif");
//	init("Fig0516(c)(BW_banreject_order4).tif");
	init("456.tif");
	fft_2d(matDst,matResult,1);
	ifft_2d(matResult,matIFFT,-1,"ifft source image",1);
	highPassFilter(matResult,matHP,900);
	ifft_2d(matHP,matIHP,-1,"ifft hp image",0);
	HFEF(matRifft,matHFEF,1);
//	fourier_trans();
//	fourier_trans();
//	fast_fourier();
//	bit_reverse();
//	zero_padding(10,5);

	//imshow("Fourier result",matDst[0]);
//	waitKey(0);











}

