#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
//#include <opencv2/highgui/highgui_c.h>
#include <iostream>
#include <math.h>

using namespace std;
using namespace cv;
float const PI = 4*atan(1);

void makeControlPoints(Mat &,int);
void makeOffset(Mat &,bool);
void warp(Mat,int);
void warpTemp(Vec2f &,Mat,int,int);
float B(float,int,int);
float calColorWork(Mat,Mat,int,int);
float calLengthWork(Mat,Mat);
Mat takeQuad(Vec2f* , int, int);
Mat takeQuads(Vec2f*, int, int);
Vec2f pixel2grid (int, int);
Vec2i grid2pixel(Vec2f);
bool compareWarp(Mat,Mat,int);
void smooth(Mat,Mat&);
float calAngleWork (Mat, Mat);
float measureAngle(Vec2f, Vec2f, Vec2f);

Mat offsets_4,offsets_25;
Mat img_dst,img_src;

int main(int argc, const char** argv) {
	if (argc!=2)
		cout << "Usage ./test <picture>" << endl;
	char pictureName1[100];
	sprintf(pictureName1,"/Users/judy/Documents/senior/SpecialProject/morphing/Implementation/img/%s1.png",argv[1]);
	char pictureName2[100];
	sprintf(pictureName2,"/Users/judy/Documents/senior/SpecialProject/morphing/Implementation/img/%s2.png",argv[1]);
	img_src = imread(pictureName1,CV_LOAD_IMAGE_UNCHANGED);
	img_dst = imread(pictureName2,CV_LOAD_IMAGE_UNCHANGED);
	img_src.convertTo(img_src,CV_32FC3,1/255.0); //convert to 1~255 scale
	img_dst.convertTo(img_dst,CV_32FC3,1/255.0);
	cout << "src size: " <<  img_src.rows << "x" << img_src.cols << endl;
	// cout << img_src.at<Vec3b>(32,22) << endl;
	cout << "dst size: " <<  img_dst.rows << "x" << img_dst.cols << endl;
	//cout << img_dst << endl;
	if (img_src.rows != img_src.cols)
		cout << "please input square pictures" << endl;
	namedWindow("src",CV_WINDOW_AUTOSIZE);
	imshow("src",img_src);
	namedWindow("dst",CV_WINDOW_AUTOSIZE);
	imshow("dst",img_dst);

	makeOffset(offsets_4,0);
	makeOffset(offsets_25,1);

	// int i = 1;
	Mat ControlPoints;

	Mat tempControlPoints;
	bool stop=false;

	for ( int i=1; ; i++) {
		cout << " i: " << i << endl;
		Mat outputImg(img_dst.rows,img_dst.cols,CV_32FC3);
		Mat recordWarp = Mat::ones(img_dst.rows,img_dst.cols,CV_8U);

		if (i!=1) {
			tempControlPoints = ControlPoints.clone();
		}

		makeControlPoints(ControlPoints,i);


		// for (int j=0; j<ControlPoints.cols; j++) {
		// 	Vec2f temp = ControlPoints.ptr<Vec2f>(0)[j];
		// 	cout << "from: " << temp  ;
		// 	warpTemp(temp,ControlPoints,i);

		// 	cout << " to: " << temp <<endl;
		// }

		cout << "CP moving " << endl;
		// cout << ControlPoints << endl << endl;

		warp(ControlPoints,i);
		cout << "CP moved end" << endl;
		cout << ControlPoints << endl << endl;
		for (int x=0; x<img_src.rows; x++)
			for (int y=0; y<img_src.cols; y++){
				Vec2f out = pixel2grid(x,y);
				warpTemp(out,ControlPoints,i,0);
				Vec2i outi = grid2pixel(out);
				outputImg.at<Vec3f>(outi[1],outi[0]) = img_src.at<Vec3f>(y,x);
				recordWarp.at<uchar>(outi[1],outi[0]) = 0;
				// cout << "outi " << outi  << " " << x << " " << y << " " << outputImg.at<Vec3f>(outi[1],outi[0]) << endl;
			}
		cout << "without smoothing done" << endl;


		smooth(recordWarp,outputImg);
		cout << "smoothing done" << endl;


		char windowName[20];
		sprintf(windowName,"warped_%d",i);
		namedWindow(windowName,CV_WINDOW_AUTOSIZE);
		imshow(windowName,outputImg);
		if (i!=1)
			stop = compareWarp(ControlPoints,tempControlPoints,i-1);
		cout << "stop? " << stop << endl;	
		if (stop)
			break;

	}

	 
//char windowName[20];
	//for ( int i=0; i< img.rows ; i+=img.rows/2 )
			//for ( int j=0; j< img.cols; j+=img.cols/2 )
			//{
					//sprintf(windowName,"grid%d",i+j);
					//Mat tile = img(Range(i,i+img.rows/2),Range(j,j+img.cols/2));
					////namedWindow("windowName",
					//imshow(windowName,tile);
					//waitKey(0);
					//destroyWindow(windowName);
			//}
	////imshow("test",img);

	// namedWindow("warped",CV_WINDOW_AUTOSIZE);
	// imshow("warped",outputImg);

	waitKey(0);
	destroyWindow("src");
	destroyWindow("dst");
	destroyWindow("warped");
	return 0;
}
void smooth(Mat recordWarp, Mat& outputImg) {
	Mat nonZero;
	findNonZero(recordWarp,nonZero);
	// cout <<"nonzero: " << nonZero.total() << endl;

	for (int z=0; z<nonZero.total(); z++) {
		// cout << "z: " << z << endl;
		Point nz = nonZero.at<Point>(z);
		// cout << "nz: " << nz.y << " " << nz.x << endl;

		int count =0;
		Vec3f sum = Vec3f(0.0,0.0,0.0);
		if (nz.y-1>=0 && nz.x-1>=0 && recordWarp.at<uchar>(nz.y-1,nz.x-1) == 0) {
			// cout << "record: " << nz.y << " " << nz.x << " " <<recordWarp.at<uchar>(nz.y-1,nz.x-1) << endl;
			sum += outputImg.at<Vec3f>(nz.y-1,nz.x-1);
			// cout << outputImg.at<Vec3f>(nz.y-1,nz.x-1) << " " << sum << endl;
			count++;
		}
		if (nz.y-1 >= 0 && recordWarp.at<uchar>(nz.y-1,nz.x) == 0) {
			sum += outputImg.at<Vec3f>(nz.y-1,nz.x);
			count++;
		}
		if (nz.y-1>=0 && nz.x+1 < img_src.cols && recordWarp.at<uchar>(nz.y-1,nz.x+1) == 0) {
			sum += outputImg.at<Vec3f>(nz.y-1,nz.x+1);
			count++;
		}
		if (nz.x-1 >= 0 && recordWarp.at<uchar>(nz.y,nz.x-1) == 0) {
			sum += outputImg.at<Vec3f>(nz.y,nz.x-1);
			count++;
		}
		if (nz.x+1 < img_src.cols && recordWarp.at<uchar>(nz.y,nz.x+1) == 0) {
			sum += outputImg.at<Vec3f>(nz.y,nz.x+1);
			count++;
		}
		if (nz.y+1 < img_src.rows && nz.x-1 >= 0 && recordWarp.at<uchar>(nz.y+1,nz.x-1) == 0) {
			sum += outputImg.at<Vec3f>(nz.y+1,nz.x-1);
			count++;
		}					
		if (nz.y+1 < img_src.rows && recordWarp.at<uchar>(nz.y+1,nz.x) == 0) {
			sum += outputImg.at<Vec3f>(nz.y+1,nz.x);
			count++;
		}					
		if (nz.y+1 < img_src.rows && nz.x+1 < img_src.cols && recordWarp.at<uchar>(nz.y+1,nz.x+1) == 0) {
			sum += outputImg.at<Vec3f>(nz.y+1,nz.x+1);
			count++;
		}

		// cout << "sum: " << sum << " count " << count << " " << endl;

		if (count != 0) {
			outputImg.at<Vec3f>(nz.y,nz.x) = sum/count;
			// cout << "sum: " << sum << " count " << count << " " << outputImg.at<Vec3f>(nz.y,nz.x) << endl;
		}
		// else
			// cout << "still no value" << endl;
	}
}
bool compareWarp(Mat A, Mat B, int i) {
	int a_points = pow(2,i+1)+1;
	int b_points = pow(2,i)+1;
	Vec2f* b = B.ptr<Vec2f>(0);
	Vec2f* a = A.ptr<Vec2f>(0);

	for (int j=0; j<a_points; j++)
		for (int k=0; k<a_points; k++) {
			Vec2f sum = 0.25 * (b[(int)(floor(j/2)*b_points+floor(k/2))]
				+ b[(int)(floor((j+1)/2)*b_points+floor(k/2))]
				+ b[(int)(floor(j/2)*b_points+floor((k+1)/2))]
				+ b[(int)(floor((j+1)/2)*b_points+floor((k+1)/2))]);
			// cout << "sum: " << sum << " new: " << a[j*a_points+k] << endl;
			if (sum != a[j*a_points+k])
				return false;
		}
	return true;
}

void makeControlPoints(Mat &ControlPoints,int i) {
	int points = pow(2,i)+1;

	// if (i==1) {
		// cout << "points: " << points << endl;
		ControlPoints = Mat::zeros(2,points*points,CV_32FC2);
		//ControlPoints = (Mat_<double>(2,2) << 1,2,3,4);
		//cout << ControlPoints << endl << endl;
		Vec2f* ControlPoints0 = ControlPoints.ptr<Vec2f>(0);
		Vec2f* ControlPoints1 = ControlPoints.ptr<Vec2f>(1);
		//cout << "ControlPoints0: " << endl << ControlPoints0[1] << endl << endl;
		//cout << "ControlPoints1: " << endl << ControlPoints1[1] << endl << endl;
		for ( int j=0; j<=pow(2,i); j++ )
			for ( int k=0; k<=pow(2,i); k++ ) {
				Vec2f pt = Vec2f(j/pow(2,i),k/pow(2,i));
				//cout << "P(" << j << "," << k << ")=" << pt << endl;
				int idx = j*points + k ;
				ControlPoints0[idx] = pt;
				ControlPoints1[idx] = pt;
			}		
	// }
	// else {
	// 	Mat tempCP = Mat::zeros(2,points*points,CV_32FC2);
	// 	Vec2f* tempCP0 = tempCP.ptr<Vec2f>(0);
	// 	Vec2f* tempCP1 = tempCP.ptr<Vec2f>(1);
	// 	for ( int j=0; j<points; j++)
	// 		for (int k=0; k<points; k++) {
	// 			int idx = j*points + k ;
	// 			Vec2f pt;
	// 			if( j%2==0 && k%2==0 ) {
	// 				int tempidx = (j/2)*(pow(2,i-1)+1)+(k/2);
	// 				cout << "tempidx: " << tempidx << endl;
	// 				pt = ControlPoints.ptr<Vec2f>(0)[tempidx];
	// 			}
	// 			else
	// 				pt = Vec2f(j/pow(2,i),k/pow(2,i));
	// 			cout << "j: " << j << " k: " << k << " pt: " << pt << endl;
	// 			tempCP0[idx] = pt;
	// 			tempCP1[idx] = pt;
	// 		}
	// 	ControlPoints = tempCP.clone();
	// }


	// cout << ControlPoints << endl << endl;
}

void makeOffset(Mat &offsets, bool trymore) {
	if ( trymore ) {
		double theta = 2*PI/8; //45 deg
		// cout << "theta: " << theta << endl;
		float unit_rotate[] = {cos(theta),sin(theta),-sin(theta),cos(theta)};
		// cout << unit_rotate[0] << endl;
		Mat unit_rot = Mat(2,2,CV_32FC1,unit_rotate);
		// cout << unit_rot << endl << endl;
		//offsets = Mat::zeros(5,5,CV_32FC2);
		offsets = Mat::zeros(1,25,CV_32FC2);

		for (int i=0; i<5; i++)
			for (int j=0; j<5; j++) {
				float _temp[] = {(i-2)*0.14,(j-2)*0.14};
				Mat temp = Mat(2,1,CV_32FC1,_temp);
				temp = unit_rot*temp;
				//cout << "temp: " << temp << endl;
				int idx = i*5+j;
				offsets.at<Vec2f>(0,idx) = Vec2f(temp.at<float>(0,0),temp.at<float>(0,1));
				//offsets.at<Vec2f>(i,j) = Vec2f(temp.at<float>(0,0),temp.at<float>(0,1));
	//offsets = unit_rot*offsets;
			}
		// cout << offsets << endl << endl;
		//cout << offsets.len << endl;
	}
	else {
		offsets = Mat::zeros(1,4,CV_32FC2);
		offsets.at<Vec2f>(0,0) = Vec2f(-0.25,-0.25);
		offsets.at<Vec2f>(0,1) = Vec2f(-0.25,0.25);
		offsets.at<Vec2f>(0,2) = Vec2f(0.25,-0.25);
		offsets.at<Vec2f>(0,3) = Vec2f(0.25,0.25);
		cout << offsets << endl << endl;
	}
}
void warp( Mat ControlPoints, int i)  {
	//Mat colorQuad = Mat::zeros(pow(2,i),pow(2,i),CV_32FC1);
	//cout << colorQuad << endl << endl;
	//for (int j=0; j<pow(2,i); j++)
			//for (int k=0; k<pow(2,i); k++)
					//calColorWork(colorQuad.at<float>(j,k),ControlPoints,j,k);

	// cout << "CP : " << ControlPoints.cols << endl;
	Vec2f* ControlPoints0 = ControlPoints.ptr<Vec2f>(0);
	Vec2f* ControlPoints1 = ControlPoints.ptr<Vec2f>(1);
	Mat offsets;
	int n;
	switch(i) {
		case 1: n=5;
				break;
		case 2: n=3;
				break;
		case 3: n=2;
				break;
		default: n=1;
				break;
	}

	// n=1;

	for (int j=0; j<ControlPoints.cols; j++) {
		int type = 0; // 1:corner 2:edge_x 3:edge_y
		// cout << "original ControlPoints: " << ControlPoints0[j] << "  " << ControlPoints1[j] << " -> ";
		if (ControlPoints0[j] == Vec2f(0,0) ||ControlPoints0[j] == Vec2f(0,1) 
				||ControlPoints0[j] == Vec2f(1,0) ||ControlPoints0[j] == Vec2f(1,1) ) {
			// cout << " no move! " << endl;
			type = 1;
			continue;
		}
		for ( int nn=1; nn<=n; nn++ ) {
			if (nn==1 && i<4)// (i==2 || i==3))
				offsets = offsets_25;
			else
				offsets = offsets_4;
			//cout << "offsets: " << offsets.cols << " " << offsets << endl << endl;
			float s = 0.2/i + 0.3/nn;
			//cout << "s: " << s << endl;

			Vec2f minCP= ControlPoints0[j];
			Mat quads0 = takeQuads(ControlPoints0,j,i);
			float Color = calColorWork(quads0,ControlPoints,i,0);
			float Length = 0;//calLengthWork(quads0);
			float Angle = 0;//calLengthWork(quads0);
			float Work = pow((256.0/img_src.rows),2)*Color + (400.0/pow(256,2))*Length + (600.0/pow(256,2))*Angle;
				cout << "c: " << pow((256.0/img_src.rows),2)*Color << endl;
				cout << "l: " << (400.0/pow(256,2))*Length  << endl;
				cout << "a: " << (600.0/pow(256,2))*Angle << endl;			

			cout << "original Work: " << Work << endl;

			for (int o=0; o<offsets.cols; o++) {
				Vec2f moved = ControlPoints0[j];
				// cout << ControlPoints0[j] << " " << ControlPoints1[j] << endl;
				// moveby Offset
				Vec2f temp = ControlPoints0[j] + s* offsets.at<Vec2f>(0,o);
				// cout << "from: " << moved << "temp: " << temp ;
				warpTemp(temp,ControlPoints,i,0);
				// moved = temp;
				if ( moved[0] == 0 || moved[0]  == 1) {
					type = 2;
					if ( moved[0] == temp[0] )
						moved = temp;
				}
				else if (moved[1] == 0 || moved[1] == 1) {
					type = 3;
					if ( moved[1] == temp[1] )
						moved = temp;
				}
				else
					moved = temp;
			
				ControlPoints1[j] = moved;
				// cout << " to: " << ControlPoints1[j] <<endl;

				// Mat quads0 = takeQuads(ControlPoints0,j,i);
				Mat quads1 = takeQuads(ControlPoints1,j,i);
				// float Work = calColorWork(quads0,ControlPoints,i);
				float movedColor = calColorWork(quads1,ControlPoints,i,1);
				float movedLength = calLengthWork(quads0,quads1);
				float movedAngle = calAngleWork(quads0,quads1) / pow(2,i);
				float movedWork = pow((256.0/img_src.rows),2)*movedColor + (400.0/pow(256,2))*movedLength + (600.0/pow(256,2))*movedAngle;
				cout << "mc: " << pow((256.0/img_src.rows),2)*movedColor << endl;
				cout << "ml: " << (400.0/pow(256,2))*movedLength  << endl;
				cout << "ma: " << (600.0/pow(256,2))*movedAngle << endl;

				cout << "work0: " << Work << " movedWork: " << movedWork << endl;
				if (movedWork < Work) {
					// cout << "cal work move!!!! " << movedWork << endl;
					// ControlPoints0[j] = moved;
					minCP = moved;
					Work = movedWork;
				}
				// else
					ControlPoints1[j] = ControlPoints0[j];

			}
			ControlPoints0[j] = minCP;
			ControlPoints1[j] = minCP;

			// cout << endl;
		}
		// cout  << ControlPoints0[j] << "  " << ControlPoints1[j] << endl;
	}	
}

float calColorWork(Mat quads, Mat ControlPoints, int i, int flag) {
	double dist=0;
	//cout << "in cal: " << quads << endl << endl;
	for (int r=0; r<quads.rows; r++)
	{
		// cout << "quads: " << quads.row(r) << endl;
		//cout << quads.at<Vec2f>(x,0) << endl;
		vector<Point2f> contour(4);
		for (int v=0; v<4; v++) {
			Vec2f temp = Vec2f(quads.at<Vec2f>(r,v));
			contour[v] = Point(temp[0]*img_src.rows,temp[1]*img_src.cols);
			//cout << "temp.x " << temp[0] << " temp.y " << temp[1] << " rows: " << img_src.rows << endl;
			// cout << contour[v] << " " ;
		}

		vector<Mat> xy(2);
		split(quads.row(r),xy);
		//cout << xy[0] << endl; // quad x
		//cout << xy[1] << endl; // quad y
		double minx,miny,maxx,maxy;
		minMaxLoc(xy[0],&minx,&maxx,NULL,NULL);
		minMaxLoc(xy[1],&miny,&maxy,NULL,NULL);
		//cout << "minx: " <<minx *img_src.rows<< " maxx: " << maxx*img_src.rows << endl;
		//cout << "miny: " <<miny*img_src.cols << " maxy: " << maxy*img_src.cols << endl;
		for (int x=(minx*img_src.rows) ; x<(maxx*img_src.rows) ; x++ ) {
			for (int y=(miny*img_src.cols) ; y<(maxy*img_src.cols); y++) {
				
				int test = pointPolygonTest(Mat(contour), Point2f(x,y), false);
				if ( test >=0 ) {
					Vec3b colorc1 = img_src.at<Vec3b>(y,x);
					// Vec3f colorc1_f = (Vec3f)colorc1 / 255.0;
					Vec2f temp = pixel2grid(x,y);
					//Vec2f temp = Vec2f(x/img_src.rows,y/img_src.cols);
					warpTemp(temp,ControlPoints,i,flag);
					//cout << "compare : " << float(x)/float(img_src.rows) << "," << float(y)/float(img_src.cols) << " " << temp[0]*img_src.rows << " " << temp[1]*img_src.cols << endl;
					
					// cout << "temp" << temp << endl;
					Vec3b colorc2 = img_dst.at<Vec3b>(temp[1]*img_src.cols,temp[0]*img_src.rows);
					// Vec3f colorc2_f = (Vec3f)colorc / 255.0;

					// colorc2 /= 255.0;

					// cout << x << " " << y << " " << grid2pixel(temp) << endl;
					dist += norm(colorc1,colorc2)/255.0;
					// if ( norm(colorc1,colorc2)/255.0 >=1 ) {
					// 	cout << "color1: " << colorc1 << " color2: " << colorc2 << endl;
					// 	cout << "norm: " << norm(colorc1,colorc2)/255.0 << endl;						
					// }

				}
				//cout << "(x,y)" << x << " " << y << " test: " << test <<endl;
			}
		}
	}
	// cout << "colorWork: " << dist << endl;
	return dist;
}

float calLengthWork(Mat quads0, Mat quads1) {
	float dist=0;
	for (int r=0; r<quads0.rows; r++) {
		dist += pow(img_src.rows*(norm(Point2f(quads0.at<Vec2f>(r,0))-Point2f(quads0.at<Vec2f>(r,1)))-norm(Point2f(quads1.at<Vec2f>(r,0))-Point2f(quads1.at<Vec2f>(r,1)))),2);
		dist += pow(img_src.rows*(norm(Point2f(quads0.at<Vec2f>(r,0))-Point2f(quads0.at<Vec2f>(r,2)))-norm(Point2f(quads1.at<Vec2f>(r,0))-Point2f(quads1.at<Vec2f>(r,2)))),2);
		dist += pow(img_src.rows*(norm(Point2f(quads0.at<Vec2f>(r,3))-Point2f(quads0.at<Vec2f>(r,1)))-norm(Point2f(quads1.at<Vec2f>(r,3))-Point2f(quads1.at<Vec2f>(r,1)))),2);
		dist += pow(img_src.rows*(norm(Point2f(quads0.at<Vec2f>(r,2))-Point2f(quads0.at<Vec2f>(r,3)))-norm(Point2f(quads1.at<Vec2f>(r,2))-Point2f(quads1.at<Vec2f>(r,3)))),2);
	}
	// cout << "lengthWork: " << dist << endl;
	return dist;
}

float calAngleWork (Mat quads0, Mat quads1) {
	float angleDiff = 0;
	for (int r=0; r<quads0.rows; r++) {
		if (quads0.at<Vec2f>(r,0) != quads1.at<Vec2f>(r,0) || 
			quads0.at<Vec2f>(r,1) != quads1.at<Vec2f>(r,1) || 
			quads0.at<Vec2f>(r,2) != quads1.at<Vec2f>(r,2))
			angleDiff += pow((measureAngle(quads0.at<Vec2f>(r,0),quads0.at<Vec2f>(r,1),quads0.at<Vec2f>(r,2))- measureAngle(quads1.at<Vec2f>(r,0),quads1.at<Vec2f>(r,1),quads1.at<Vec2f>(r,2))),2);
		if (quads0.at<Vec2f>(r,1) != quads1.at<Vec2f>(r,1) || 
			quads0.at<Vec2f>(r,0) != quads1.at<Vec2f>(r,0) || 
			quads0.at<Vec2f>(r,3) != quads1.at<Vec2f>(r,3))
			angleDiff += pow((measureAngle(quads0.at<Vec2f>(r,1),quads0.at<Vec2f>(r,0),quads0.at<Vec2f>(r,3))- measureAngle(quads1.at<Vec2f>(r,1),quads1.at<Vec2f>(r,0),quads1.at<Vec2f>(r,3))),2);	
		if (quads0.at<Vec2f>(r,3) != quads1.at<Vec2f>(r,3) || 
			quads0.at<Vec2f>(r,1) != quads1.at<Vec2f>(r,1) || 
			quads0.at<Vec2f>(r,2) != quads1.at<Vec2f>(r,2))
			angleDiff += pow((measureAngle(quads0.at<Vec2f>(r,3),quads0.at<Vec2f>(r,1),quads0.at<Vec2f>(r,2))- measureAngle(quads1.at<Vec2f>(r,3),quads1.at<Vec2f>(r,1),quads1.at<Vec2f>(r,2))),2);		
		if (quads0.at<Vec2f>(r,2) != quads1.at<Vec2f>(r,2) || 
			quads0.at<Vec2f>(r,0) != quads1.at<Vec2f>(r,0) || 
			quads0.at<Vec2f>(r,3) != quads1.at<Vec2f>(r,3))
			angleDiff += pow((measureAngle(quads0.at<Vec2f>(r,2),quads0.at<Vec2f>(r,0),quads0.at<Vec2f>(r,3))- measureAngle(quads1.at<Vec2f>(r,2),quads1.at<Vec2f>(r,0),quads1.at<Vec2f>(r,3))),2);
	}
	// cout << "angleDiff:  " << angleDiff << endl;
	return angleDiff;
}
float measureAngle(Vec2f a, Vec2f b , Vec2f c) {
	// cout << "a: " << a << " b: " << b << " c: " << c << endl;
	// cout << b[0]-a[0] << " " << b[1]-a[1] << endl;
	Vec2f v1 = Vec2f(b[0]-a[0],b[1]-a[1]);
	Vec2f v2 = Vec2f(c[0]-a[0],c[1]-a[1]);
	// cout << "v1: " << v1 << " v2: " << v2 << endl;
	float len1 = sqrt(v1[0]*v1[0]+v1[1]*v1[1]);
	float len2 = sqrt(v2[0]*v2[0]+v2[1]*v2[1]);
	float dot = v1[0]*v2[0]+ v1[1]*v2[1] ;
	float temp = dot / (len1*len2);
	// cout << "cos: " << temp << endl;
	if (temp == 1)
		return 0.0;
	else if (temp == -1)
		return PI;
	else
		return acos(temp);
}
Vec2f pixel2grid ( int x, int y) {
	return Vec2f(float(x)/float(img_src.rows),float(y)/float(img_src.cols));
}
Vec2i grid2pixel ( Vec2f a ) {
	return Vec2i(a[0]*img_src.rows,a[1]*img_src.cols);
}


Mat takeQuad( Vec2f* ControlPoints0, int x, int points) {//int j, int k, int i) //left button point {
	Mat quad = Mat::zeros(1,4,CV_32FC2);
	quad.at<Vec2f>(0,0) = ControlPoints0[x];
	quad.at<Vec2f>(0,1) = ControlPoints0[x+1];
	quad.at<Vec2f>(0,2) = ControlPoints0[x+points];
	quad.at<Vec2f>(0,3) = ControlPoints0[x+1+points];
	return quad;
}
Mat takeQuads(Vec2f* ControlPoints0,int x,int i) {
	int points = pow(2,i)+1;
	Mat quads;

	//cout << ControlPoints0  << endl;
	//Vec2f* ControlPoints0 = ControlPoints.ptr<Vec2f>(0);
	if ( ControlPoints0[x](0) ==0 ) {
		quads = Mat(2,4,CV_32FC2);
		Mat quad1,quad2;
		quad1 = takeQuad(ControlPoints0,x-1,points);
		quad2 = takeQuad(ControlPoints0,x,points);
		quad1.copyTo(quads.row(0));
		quad2.copyTo(quads.row(1));
	}
	else if ( ControlPoints0[x](0) == 1) {
		quads = Mat(2,4,CV_32FC2);
		Mat quad1,quad2;
		quad1 = takeQuad(ControlPoints0,x-points-1,points);
		quad2 = takeQuad(ControlPoints0,x-points,points);
		quad1.copyTo(quads.row(0));
		quad2.copyTo(quads.row(1));
	}
	else if ( ControlPoints0[x](1) == 0) {
		quads = Mat(2,4,CV_32FC2);
		Mat quad1,quad2;
		quad1 = takeQuad(ControlPoints0,x-points,points);
		quad2 = takeQuad(ControlPoints0,x,points);
		quad1.copyTo(quads.row(0));
		quad2.copyTo(quads.row(1));
	}
	else if ( ControlPoints0[x](1) == 1) {
		quads = Mat(2,4,CV_32FC2);
		Mat quad1,quad2;
		quad1 = takeQuad(ControlPoints0,x-points-1,points);
		quad2 = takeQuad(ControlPoints0,x-1,points);
		quad1.copyTo(quads.row(0));
		quad2.copyTo(quads.row(1));
	}
	else {
		quads = Mat(4,4,CV_32FC2);
		Mat quad1,quad2,quad3,quad4;
		quad1 = takeQuad(ControlPoints0,x-points-1,points);
		quad2 = takeQuad(ControlPoints0,x-points,points);
		quad3 = takeQuad(ControlPoints0,x-1,points);
		quad4 = takeQuad(ControlPoints0,x,points);
		quad1.copyTo(quads.row(0));
		quad2.copyTo(quads.row(1));
		quad3.copyTo(quads.row(2));
		quad4.copyTo(quads.row(3));
	}
	// cout << "quads: " <<quads << endl << endl;
	return quads;
}

void warpTemp(Vec2f &v, Mat ControlPoints,int i,int flag) { 
	Vec2f r= Vec2f(0,0);
	for (int j=0; j<=pow(2,i); j++)
		for (int k=0; k<=pow(2,i); k++) {
			int idx = j*(pow(2,i)+1) + k;
			float Bj = B(v[0],j,i);
			float Bk = B(v[1],k,i);
			// cout << "Bj: " << Bj << " Bk: " << Bk << endl;
			// cout << "CP( " << j << "," << k << "): " << ControlPoints.ptr<Vec2f>(0)[idx] << endl;
			// r+= Vec2f( Bj * ControlPoints.ptr<Vec2f>(0)[idx][0] , Bk * ControlPoints.ptr<Vec2f>(0)[idx][1]) ;
			if (flag) 
				r+= Bj * Bk * ControlPoints.ptr<Vec2f>(1)[idx] ;
			else
				r+= Bj * Bk * ControlPoints.ptr<Vec2f>(0)[idx] ;
			// cout << "r: " << r << endl;
			// cout << "------------" << endl;
		}
	v = r;
}

float B(float x,int j,int i) {
	float basis = pow(2,i);
	//cout << "in B : " << x << " " << j<< " " << i << endl;
	if ( x >= (float(j-1)/basis) && x<= (float(j)/basis) ) {
		// cout << "left" << endl;
		return (basis*x-j+1);
	}
	else if ( x> (float(j)/basis) && x <= (float(j+1)/basis) ) {
		// cout << "right" << endl;
		return (j+1-basis*x);
	}
	else
		return 0;
}
