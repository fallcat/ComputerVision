// Source file for image class



// Include files

#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"
#include <algorithm>
#include <vector>
#include <cmath>



////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////


R2Image::
R2Image(void)
  : pixels(NULL),
    npixels(0),
    width(0),
    height(0)
{
}



R2Image::
R2Image(const char *filename)
  : pixels(NULL),
    npixels(0),
    width(0),
    height(0)
{
  // Read image
  Read(filename);
}



R2Image::
R2Image(int width, int height)
  : pixels(NULL),
    npixels(width * height),
    width(width),
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);
}



R2Image::
R2Image(int width, int height, const R2Pixel *p)
  : pixels(NULL),
    npixels(width * height),
    width(width),
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++)
    pixels[i] = p[i];
}



R2Image::
R2Image(const R2Image& image)
  : pixels(NULL),
    npixels(image.npixels),
    width(image.width),
    height(image.height)

{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++)
    pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
  // Free image pixels
  if (pixels) delete [] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
  // Delete previous pixels
  if (pixels) { delete [] pixels; pixels = NULL; }

  // Reset width and height
  npixels = image.npixels;
  width = image.width;
  height = image.height;

  // Allocate new pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels
  for (int i = 0; i < npixels; i++)
    pixels[i] = image.pixels[i];

  // Return image
  return *this;
}


void R2Image::
svdTest(void)
{
	// fit a 2D conic to five points
	R2Point p1(1.2,3.5);
	R2Point p2(2.1,2.2);
	R2Point p3(0.2,1.6);
	R2Point p4(0.0,0.5);
	R2Point p5(-0.2,4.2);

	// build the 5x6 matrix of equations
	double** linEquations = dmatrix(1,5,1,6);

	linEquations[1][1] = p1[0]*p1[0];
	linEquations[1][2] = p1[0]*p1[1];
	linEquations[1][3] = p1[1]*p1[1];
	linEquations[1][4] = p1[0];
	linEquations[1][5] = p1[1];
	linEquations[1][6] = 1.0;

	linEquations[2][1] = p2[0]*p2[0];
	linEquations[2][2] = p2[0]*p2[1];
	linEquations[2][3] = p2[1]*p2[1];
	linEquations[2][4] = p2[0];
	linEquations[2][5] = p2[1];
	linEquations[2][6] = 1.0;

	linEquations[3][1] = p3[0]*p3[0];
	linEquations[3][2] = p3[0]*p3[1];
	linEquations[3][3] = p3[1]*p3[1];
	linEquations[3][4] = p3[0];
	linEquations[3][5] = p3[1];
	linEquations[3][6] = 1.0;

	linEquations[4][1] = p4[0]*p4[0];
	linEquations[4][2] = p4[0]*p4[1];
	linEquations[4][3] = p4[1]*p4[1];
	linEquations[4][4] = p4[0];
	linEquations[4][5] = p4[1];
	linEquations[4][6] = 1.0;

	linEquations[5][1] = p5[0]*p5[0];
	linEquations[5][2] = p5[0]*p5[1];
	linEquations[5][3] = p5[1]*p5[1];
	linEquations[5][4] = p5[0];
	linEquations[5][5] = p5[1];
	linEquations[5][6] = 1.0;

	printf("\n Fitting a conic to five points:\n");
	printf("Point #1: %f,%f\n",p1[0],p1[1]);
	printf("Point #2: %f,%f\n",p2[0],p2[1]);
	printf("Point #3: %f,%f\n",p3[0],p3[1]);
	printf("Point #4: %f,%f\n",p4[0],p4[1]);
	printf("Point #5: %f,%f\n",p5[0],p5[1]);

	// compute the SVD
	double singularValues[7]; // 1..6
	double** nullspaceMatrix = dmatrix(1,6,1,6);
	svdcmp(linEquations, 5, 6, singularValues, nullspaceMatrix);

	// get the result
	printf("\n Singular values: %f, %f, %f, %f, %f, %f\n",singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6]);

	// find the smallest singular value:
	int smallestIndex = 1;
	for(int i=2;i<7;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

	// solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
	printf("Conic coefficients: %f, %f, %f, %f, %f, %f\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],nullspaceMatrix[3][smallestIndex],nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex]);

	// make sure the solution is correct:
	printf("Equation #1 result: %f\n",	p1[0]*p1[0]*nullspaceMatrix[1][smallestIndex] +
										p1[0]*p1[1]*nullspaceMatrix[2][smallestIndex] +
										p1[1]*p1[1]*nullspaceMatrix[3][smallestIndex] +
										p1[0]*nullspaceMatrix[4][smallestIndex] +
										p1[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #2 result: %f\n",	p2[0]*p2[0]*nullspaceMatrix[1][smallestIndex] +
										p2[0]*p2[1]*nullspaceMatrix[2][smallestIndex] +
										p2[1]*p2[1]*nullspaceMatrix[3][smallestIndex] +
										p2[0]*nullspaceMatrix[4][smallestIndex] +
										p2[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #3 result: %f\n",	p3[0]*p3[0]*nullspaceMatrix[1][smallestIndex] +
										p3[0]*p3[1]*nullspaceMatrix[2][smallestIndex] +
										p3[1]*p3[1]*nullspaceMatrix[3][smallestIndex] +
										p3[0]*nullspaceMatrix[4][smallestIndex] +
										p3[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #4 result: %f\n",	p4[0]*p4[0]*nullspaceMatrix[1][smallestIndex] +
										p4[0]*p4[1]*nullspaceMatrix[2][smallestIndex] +
										p4[1]*p4[1]*nullspaceMatrix[3][smallestIndex] +
										p4[0]*nullspaceMatrix[4][smallestIndex] +
										p4[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	printf("Equation #5 result: %f\n",	p5[0]*p5[0]*nullspaceMatrix[1][smallestIndex] +
										p5[0]*p5[1]*nullspaceMatrix[2][smallestIndex] +
										p5[1]*p5[1]*nullspaceMatrix[3][smallestIndex] +
										p5[0]*nullspaceMatrix[4][smallestIndex] +
										p5[1]*nullspaceMatrix[5][smallestIndex] +
										nullspaceMatrix[6][smallestIndex]);

	R2Point test_point(0.34,-2.8);

	printf("A point off the conic: %f\n",	test_point[0]*test_point[0]*nullspaceMatrix[1][smallestIndex] +
											test_point[0]*test_point[1]*nullspaceMatrix[2][smallestIndex] +
											test_point[1]*test_point[1]*nullspaceMatrix[3][smallestIndex] +
											test_point[0]*nullspaceMatrix[4][smallestIndex] +
											test_point[1]*nullspaceMatrix[5][smallestIndex] +
											nullspaceMatrix[6][smallestIndex]);

	return;
}



////////////////////////////////////////////////////////////////////////
// Image processing functions
// YOU IMPLEMENT THE FUNCTIONS IN THIS SECTION
////////////////////////////////////////////////////////////////////////

// Per-pixel Operations ////////////////////////////////////////////////

void R2Image::
Brighten(double factor)
{
  // Brighten the image by multiplying each pixel component by the factor.
  // This is implemented for you as an example of how to access and set pixels
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Pixel(i,j) *= factor;
      Pixel(i,j).Clamp();
    }
  }
}

void R2Image::
SobelX(void)
{

  // Apply the Sobel oprator to the image in X direction
  int s;
  s = 1;
  int weightSize;
  weightSize = 3;
  double weightArray[3][3] = {
    {1, 0, -1},
    {2, 0, -2},
    {1, 0, -1}
  };

  R2Pixel *tempPixels;
  tempPixels = new R2Pixel[width*height];

  for (int i = s; i < width-s; i++) {
    for (int j = s;  j < height-s; j++) {
      R2Pixel* val;
      val = new R2Pixel(0,0,0,0);
      for (int k = -s; k <= s; k ++) {
        for (int l = -s; l <= s; l ++) {
          *val = *val + Pixel(i+k,j+l)*weightArray[(s+l)][(s+k)];
        }
      }
      tempPixels[i*height+j] = *val;
    }
  }

  for (int i = s; i < width-s; i++) {
    for (int j = s;  j < height-s; j++) {
      Pixel(i,j) = tempPixels[i*height+j];
      //Pixel(i,j).Clamp();
    }
  }

  free(tempPixels);

}

void R2Image::
SobelY(void)
{
	// Apply the Sobel oprator to the image in Y direction
  int s;
  s = 1;
  int weightSize;
  weightSize = 3;
  double weightArray[3][3] = {
    {1, 2, 1},
    {0, 0, 0},
    {-1, -2, -1}
  };

  R2Pixel *tempPixels;
  tempPixels = new R2Pixel[width*height];

  for (int i = s; i < width-s; i++) {
    for (int j = s;  j < height-s; j++) {
      R2Pixel* val;
      val = new R2Pixel(0,0,0,0);
      for (int k = -s; k <= s; k ++) {
        for (int l = -s; l <= s; l ++) {
          *val = *val + Pixel(i+k,j+l)*weightArray[(s+l)][(s+k)];
        }
      }
      tempPixels[i*height+j] = *val;
    }
  }

  for (int i = s; i < width-s; i++) {
    for (int j = s;  j < height-s; j++) {
      Pixel(i,j) = tempPixels[i*height+j];
      //Pixel(i,j).Clamp();
    }
  }

  free(tempPixels);
}

void R2Image::
LoG(void)
{
  // Apply the LoG oprator to the image

  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "LoG() not implemented\n");
}



void R2Image::
ChangeSaturation(double factor)
{
  // Changes the saturation of an image
  // Find a formula that changes the saturation without affecting the image brightness

  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "ChangeSaturation(%g) not implemented\n", factor);
}


// Linear filtering ////////////////////////////////////////////////
void R2Image::
Blur(double sigma)
{
  // Gaussian blur of the image. Separable solution is preferred

  int s;
  s = (int)sigma*3;
  int weightSize;
  weightSize = s*2+1;
  double * weightArray = (double*)malloc(weightSize*sizeof(double));
  double weights = 0;
  for (int i = -s; i <= s; i++) {
    weightArray[s+i] = 1/sqrt(2*M_PI*pow(sigma,2))*exp(-pow(i,2)/(2*pow(sigma,2)));
    weights += weightArray[s+i];
  }

  R2Pixel *tempPixels1;
  tempPixels1 = new R2Pixel[width*height];

  for (int i = 0; i < width; i++) {
    double tempWeights = 0;
    if (i <= s || i >= width - s) {
      for (int lx = -s; lx <= s; lx ++) {
        if (i+lx >= 0 && i+lx < width) {
          tempWeights += weightArray[s+lx];
        }
      }
    }
    else {
      tempWeights = weights;
    }
    for (int j = 0;  j < height; j++) {
      R2Pixel* val;
      val = new R2Pixel(0,0,0,0);
      for (int lx = -s; lx <= s; lx ++) {
        if (i+lx >= 0 && i+lx < width) {
          *val += Pixel(i+lx,j)*weightArray[s+lx];
        }
      }
      tempPixels1[i*height+j] = *val/tempWeights;
    }
  }


  for (int j = 0;  j < height; j++) {
    double tempWeights = 0;
    if (j <= s || j >= height - s) {
      for (int ly = -s; ly <= s; ly ++) {
        if (j+ly >= 0 && j+ly < height) {
          tempWeights += weightArray[s+ly];
        }
      }
    }
    else {
      tempWeights = weights;
    }
    for (int i = 0; i < width; i++) {
      R2Pixel* val;
      val = new R2Pixel(0,0,0,0);
      for (int ly = -s; ly <= s; ly ++) {
        if (j+ly >= 0 && j+ly < height) {
          *val = *val + tempPixels1[i*height+j+ly]*weightArray[s+ly];
        }
      }
      Pixel(i,j) = *val / tempWeights;
      //Pixel(i,j).Clamp();
    }
  }

  free(weightArray);
  free(tempPixels1);
}

void R2Image::
HighPassSharpen(double sigma)
{
  // Sharpening filter using high pass filter

  R2Image blurredPicture(*this);
  blurredPicture.Blur(sigma);

  R2Pixel *highPassPicture;
  highPassPicture = new R2Pixel[width*height];

  R2Pixel *grayPixel;
  grayPixel = new R2Pixel(0.5,0.5,0.5,0.5);

  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      highPassPicture[i*height+j] = (Pixel(i,j) - blurredPicture.Pixel(i,j))+*grayPixel;
    }
  }

  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      double temprgba[4];
      for (int k = 0; k < 4; k ++) {
        // Overlay
        if (highPassPicture[i*height+j].Component(k) < 0.5) {
          temprgba[k] = 2*highPassPicture[i*height+j].Component(k)*Pixel(i,j).Component(k);
        }
        else {
          temprgba[k] = 1 - 2*(1-highPassPicture[i*height+j].Component(k))*(1-Pixel(i,j).Component(k));
        }
      }

      Pixel(i,j).Reset(temprgba[0], temprgba[1], temprgba[2], temprgba[3]);

      Pixel(i,j).Clamp();
    }
  }

  free(highPassPicture);
}


void R2Image::
Harris(double sigma)
{
    // Harris corner detector. Make use of the previously developed filters, such as the Gaussian blur filter
	// Output should be 50% grey at flat regions, white at corners and black/dark near edges

  double alpha = 0.04;

  R2Image Img1(*this);
  R2Image Img2(*this);
  R2Image Img3(*this);

  Img1.SobelX();
  Img2.SobelY();
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Img3.Pixel(i,j) = Img1.Pixel(i,j) * Img2.Pixel(i,j);
      Img1.Pixel(i,j) = Img1.Pixel(i,j) * Img1.Pixel(i,j);
      Img2.Pixel(i,j) = Img2.Pixel(i,j) * Img2.Pixel(i,j);
    }
  }

  Img1.Blur(sigma);
  Img2.Blur(sigma);
  Img3.Blur(sigma);

  R2Pixel *halfPixel = new R2Pixel(0.5,0.5,0.5,0.5);

  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Pixel(i,j) = Img1.Pixel(i,j)*Img2.Pixel(i,j) - Img3.Pixel(i,j)*Img3.Pixel(i,j) - alpha*((Img1.Pixel(i,j)+Img2.Pixel(i,j))*(Img1.Pixel(i,j)+Img2.Pixel(i,j)));
      Pixel(i,j) += *halfPixel;
      Pixel(i,j).Clamp();
    }
  }

  free(halfPixel);

}

void R2Image::FeatureDetection(double sigma)
{
  // Harris

  R2Image harris(*this);

  double alpha = 0.04;

  R2Image Img1(*this);
  R2Image Img2(*this);
  R2Image Img3(*this);

  Img1.SobelX();
  Img2.SobelY();
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Img3.Pixel(i,j) = Img1.Pixel(i,j) * Img2.Pixel(i,j);
      Img1.Pixel(i,j) *= Img1.Pixel(i,j);
      Img2.Pixel(i,j) *= Img2.Pixel(i,j);
    }
  }

  Img1.Blur(sigma);
  Img2.Blur(sigma);
  Img3.Blur(sigma);

  //R2Pixel *halfPixel = new R2Pixel(0.5,0.5,0.5,0.5);

  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      harris.Pixel(i,j) = Img1.Pixel(i,j)*Img2.Pixel(i,j) - Img3.Pixel(i,j)*Img3.Pixel(i,j) - alpha*((Img1.Pixel(i,j)+Img2.Pixel(i,j))*(Img1.Pixel(i,j)+Img2.Pixel(i,j)));
      harris.Pixel(i,j).Reset(harris.Pixel(i,j)[0],harris.Pixel(i,j)[1],harris.Pixel(i,j)[2],harris.Pixel(i,j)[3]);
      // don't clamp here!
    }
  }

  // feature detection: detect the highest n features
  bool *features;
  features = new bool[width*height];

  std::vector<Feature> featureVec;
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      features[i*height+j] = false;
      featureVec.push_back(Feature(i, j, harris.Pixel(i,j)));
    }
  }

  std::sort (featureVec.begin(), featureVec.end());

  int total = 150;
  int count = total;
  int minDistance = 10;
  int *chosenFeatures;
  chosenFeatures = new int[2*count];
  while (count > 0 || featureVec.empty()) {
    bool flag = false;
    Feature feature = featureVec.back();
    featureVec.pop_back();
    for (int i = feature.centerX-minDistance; i <= feature.centerX+minDistance; i++){
      for (int j = feature.centerY-minDistance; j <= feature.centerY+minDistance; j++){
        if (i >= 0 && i < width && j >= 0 && j < height) {
          if (features[i*height+j]) {
            flag = true;
          }
        }
        if (flag) {
          break;
        }
      }
      if (flag) {
        break;
      }
    }

    if (!flag) {
      count --;
      features[feature.centerX*height+feature.centerY] = true;
      chosenFeatures[count*2] = feature.centerX;
      chosenFeatures[count*2+1] = feature.centerY;
      //printf("x: %d, y: %d, count: %d", chosenFeatures[count*2], chosenFeatures[count*2+1], count);

    }
  }

  // mark on the original graph
  // int markOffsets[45][2] = {
  //   {-1,4},{0,4},{1,4},
  //   {-1,3},{0,3},{1,3},
  //   {-1,2},{0,2},{1,2},
  //   {-4,1},{-3,1},{-2,1},{-1,1},{0,1},{1,1},{2,1},{3,1},{4,1},
  //   {-4,0},{-3,0},{-2,0},{-1,0},{0,0},{1,0},{2,0},{3,0},{4,0},
  //   {-4,-1},{-3,-1},{-2,-1},{-1,-1},{0,-1},{1,-1},{2,-1},{3,-1},{4,-1},
  //   {-1,-2},{0,-2},{1,-2},
  //   {-1,-3},{0,-3},{1,-3},
  //   {-1,-4},{0,-4},{1,-4}
  // };
  //int marks =  sizeof markOffsets / sizeof markOffsets[0];
  for (int i = 0; i < total; i++) {
    int x = chosenFeatures[2*i];
    int y = chosenFeatures[2*i+1];
    //printf("x: %d, y: %d", x, y);
    //for (int pair = 0; pair < marks ; pair ++) {
    for (int m = -5; m <= 5; m++) {
      for (int n = -5; n <= 5; n++) {
        if (x+m >= 0 && x+m < width && y+n >= 0 && y+n < width && (m<-2 || m>2 || n<-2 || n>2)) {
          Pixel(x+m,y+n).Reset(1,0,0,0.5);
          Pixel(x+m,y+n).Clamp();
        }
      }
    }
      // if (x+markOffsets[pair][0] >= 0 && x+markOffsets[pair][0] < width && y+markOffsets[pair][1] >= 0 && y+markOffsets[pair][1] < width) {
      //   Pixel(x+markOffsets[pair][0],y+markOffsets[pair][1]).Reset(1,0,0,0);
      //   Pixel(x+markOffsets[pair][0],y+markOffsets[pair][1]).Clamp();
      // }
    //}
  }

  free(chosenFeatures);
  free(features);

}

void R2Image::ScaledFeatureDetection()
{
  // // Harris
  //
  // R2Image harris(*this);
  // double *sigmas;
  // sigmas = new double[width*height];
  // double initialSigma = sqrt(width*height)/10;
  //
  // double alpha = 0.04;
  //
  // for (double sigma = initialSigma; sigma > initialSigma/10; sigma -= initialSigma/10) {
  //   R2Image Img1(*this);
  //   R2Image Img2(*this);
  //   R2Image Img3(*this);
  //
  //   Img1.SobelX();
  //   Img2.SobelY();
  //   for (int i = 0; i < width; i++) {
  //     for (int j = 0;  j < height; j++) {
  //       Img3.Pixel(i,j) = Img1.Pixel(i,j) * Img2.Pixel(i,j);
  //       Img1.Pixel(i,j) *= Img1.Pixel(i,j);
  //       Img2.Pixel(i,j) *= Img2.Pixel(i,j);
  //     }
  //   }
  //
  //   Img1.Blur(sigma);
  //   Img2.Blur(sigma);
  //   Img3.Blur(sigma);
  //
  //   //R2Pixel *halfPixel = new R2Pixel(0.5,0.5,0.5,0.5);
  //
  //   for (int i = 0; i < width; i++) {
  //     for (int j = 0;  j < height; j++) {
  //       if (sigma == initialSigma) {
  //         harris.Pixel(i,j) = Img1.Pixel(i,j)*Img2.Pixel(i,j) - Img3.Pixel(i,j)*Img3.Pixel(i,j) - alpha*((Img1.Pixel(i,j)+Img2.Pixel(i,j))*(Img1.Pixel(i,j)+Img2.Pixel(i,j)));
  //       }
  //       else {
  //         R2Pixel temp = Img1.Pixel(i,j)*Img2.Pixel(i,j) - Img3.Pixel(i,j)*Img3.Pixel(i,j) - alpha*((Img1.Pixel(i,j)+Img2.Pixel(i,j))*(Img1.Pixel(i,j)+Img2.Pixel(i,j)));
  //         if
  //       }
  //
  //       harris.Pixel(i,j).Reset(harris.Pixel(i,j)[0],harris.Pixel(i,j)[1],harris.Pixel(i,j)[2],harris.Pixel(i,j)[3]);
  //       // don't clamp here!
  //     }
  //   }
  // }
  //
  //
  //
  // // feature detection: detect the highest n features
  // bool *features;
  // features = new bool[width*height];
  //
  // std::vector<Feature> featureVec;
  // for (int i = 0; i < width; i++) {
  //   for (int j = 0;  j < height; j++) {
  //     features[i*height+j] = false;
  //     featureVec.push_back(Feature(i, j, harris.Pixel(i,j)));
  //   }
  // }
  //
  // std::sort (featureVec.begin(), featureVec.end());
  //
  // int total = 150;
  // int count = total;
  // int minDistance = 10;
  // int *chosenFeatures;
  // chosenFeatures = new int[2*count];
  // while (count > 0 || featureVec.empty()) {
  //   bool flag = false;
  //   Feature feature = featureVec.back();
  //   featureVec.pop_back();
  //   for (int i = feature.centerX-minDistance; i <= feature.centerX+minDistance; i++){
  //     for (int j = feature.centerY-minDistance; j <= feature.centerY+minDistance; j++){
  //       if (i >= 0 && i < width && j >= 0 && j < height) {
  //         if (features[i*height+j]) {
  //           flag = true;
  //         }
  //       }
  //       if (flag) {
  //         break;
  //       }
  //     }
  //     if (flag) {
  //       break;
  //     }
  //   }
  //
  //   if (!flag) {
  //     count --;
  //     features[feature.centerX*height+feature.centerY] = true;
  //     chosenFeatures[count*2] = feature.centerX;
  //     chosenFeatures[count*2+1] = feature.centerY;
  //     //printf("x: %d, y: %d, count: %d", chosenFeatures[count*2], chosenFeatures[count*2+1], count);
  //
  //   }
  // }
  //
  // for (int i = 0; i < total; i++) {
  //   int x = chosenFeatures[2*i];
  //   int y = chosenFeatures[2*i+1];
  //   //printf("x: %d, y: %d", x, y);
  //   //for (int pair = 0; pair < marks ; pair ++) {
  //   for (int m = -5; m <= 5; m++) {
  //     for (int n = -5; n <= 5; n++) {
  //       if (x+m >= 0 && x+m < width && y+n >= 0 && y+n < width && (m<-2 || m>2 || n<-2 || n>2)) {
  //         Pixel(x+m,y+n).Reset(1,0,0,0.5);
  //         Pixel(x+m,y+n).Clamp();
  //       }
  //     }
  //   }
  // }
  //
  // free(chosenFeatures);
  // free(features);

}


void R2Image::
Sharpen()
{
  // Sharpen an image using a linear filter. Use a kernel of your choosing.
  int s;
  s = 1;
  int weightSize;
  weightSize = 3;
  double weightArray[3][3] = {
    {-1, -1, -1},
    {-1, 9, -1},
    {-1, -1, -1}
  };

  R2Pixel *tempPixels;
  tempPixels = new R2Pixel[width*height];

  for (int i = s; i < width-s; i++) {
    for (int j = s;  j < height-s; j++) {
      R2Pixel* val;
      val = new R2Pixel(0,0,0,0);
      for (int k = -s; k <= s; k ++) {
        for (int l = -s; l <= s; l ++) {
          *val = *val + Pixel(i+k,j+l)*weightArray[(s+k)][(s+l)];
        }
      }
      tempPixels[i*height+j] = *val;
    }
  }

  for (int i = s; i < width-s; i++) {
    for (int j = s;  j < height-s; j++) {
      Pixel(i,j) = tempPixels[i*height+j];
      Pixel(i,j).Clamp();
    }
  }

  free(tempPixels);
}

void R2Image::
LensDistortion()
{
  // lensDistortion

  R2Pixel *tempPixels;
  tempPixels = new R2Pixel[width*height];
  double scale = 2;

  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      int x = i - width/2;
      int y = j - height/2;
      double rd = sqrt(pow((double)x,2)+pow((double)y,2));
      double percent = sqrt(1-(pow(rd/width,1.5)));
      //1-(rd/width + 0.1*pow(rd/width,3));
      //1+((0.5-pow(rd/width,1.5))/0.5);
      int xold = (int)x/percent/scale;
      int yold = (int)y/percent/scale;
      int iold = width/2+xold;
      int jold = height/2+yold;
      if (iold < 0) {
        iold = 0;
      }
      if (iold >= width) {
        iold = width-1;
      }
      if (jold < 0) {
        jold = 0;
      }
      if (jold >= height) {
        jold = height-1;
      }
      tempPixels[i*height+j] = Pixel(iold,jold);
    }
  }

  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Pixel(i,j) = tempPixels[i*height+j];
      Pixel(i,j).Clamp();
    }
  }

  free(tempPixels);
}

void R2Image::
Median(int radius)
{
  // Median filter

  R2Image tempPixels(*this);

  for (int i = radius; i < width-radius; i++) {
    for (int j = radius;  j < height-radius; j++) {
      double temprgba[4];
      for (int color = 0; color < 4; color ++) {
        double sortArray[(2*radius+1)*(2*radius+1)];
        for (int k = -radius; k <= radius; k ++) {
          for (int l = -radius; l <= radius; l ++) {
            sortArray[(radius+k)*(2*radius+1)+radius+l] = tempPixels.Pixel(i+k,j+l).Component(color);
            // for (int s = (radius+k)*(2*radius+1)+radius+l; s > 0; s--) {
            //   if (sortArray[s] > sortArray[s-1]) {
            //     double tempColor = sortArray[s-1];
            //     sortArray[s-1] = sortArray[s];
            //     sortArray[s] = tempColor;
            //   }
            // }
          }
        }
        std::vector<double> sortVector (sortArray, sortArray+(2*radius+1)*(2*radius+1));
        std::sort (sortVector.begin(), sortVector.begin()+(2*radius+1)*(2*radius+1));
        temprgba[color] = sortVector.at(((2*radius+1)*(2*radius+1)+1)/2);
      }
      Pixel(i,j).Reset(temprgba[0],temprgba[1],temprgba[2],temprgba[3]);
      Pixel(i,j).Clamp();
    }
  }

}

void R2Image::
Bilateral(double sigmaS, double sigmaR)
{
  // Bilateral blur

  int s;
  s = (int)sigmaS*3;
  int weightSize;
  weightSize = s*2+1;
  double * spaceWeightArray = (double*)malloc(weightSize*sizeof(double));
  for (int i = -s; i <= s; i++) {
    spaceWeightArray[s+i] = 1/sqrt(2*M_PI*pow(sigmaS,2))*exp(-pow(i,2)/(2*pow(sigmaS,2)));
  }
  double * rangeWeightArray = (double*)malloc(201*sizeof(double));
  for (int i = 0; i < 201; i++) {
    rangeWeightArray[i] = 1/sqrt(2*M_PI*pow(sigmaR,2))*exp(-pow(((double)i-100)/100,2)/(2*pow(sigmaR,2)));
  }

  R2Pixel *tempPixels1;
  tempPixels1 = new R2Pixel[width*height];

  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      double temprgba[4] = {0, 0, 0, 0};
      for (int color = 0; color < 4; color ++) {
        double totalWeights = 0;
        for (int lx = -s; lx <= s; lx ++) {
          if (i+lx >= 0 && i+lx < width) {
            double newWeight = spaceWeightArray[s+lx]*rangeWeightArray[(int)((Pixel(i+lx,j).Component(color)-Pixel(i,j).Component(color))*100+100)];
            //printf("x weight: %f", newWeight);
            temprgba[color] += Pixel(i+lx,j).Component(color)*newWeight;
            totalWeights += newWeight;
          }
        }
        if (totalWeights != 0) {
          temprgba[color] /= totalWeights;
        }
      }
      tempPixels1[i*height+j].Reset(temprgba[0],temprgba[1],temprgba[2],temprgba[3]);
    }
  }


  for (int j = 0;  j < height; j++) {
    for (int i = 0; i < width; i++) {
      double temprgba[4] = {0, 0, 0, 0};
      for (int color = 0; color < 4; color ++) {
        double totalWeights = 0;
        for (int ly = -s; ly <= s; ly ++) {
          if (j+ly >= 0 && j+ly < height) {
            double newWeight = spaceWeightArray[s+ly]*rangeWeightArray[(int)((Pixel(i,j+ly).Component(color)-Pixel(i,j).Component(color))*100+100)];
            //printf("y weight: %f", newWeight);
            temprgba[color] += tempPixels1[i*height+j+ly].Component(color)*newWeight;
            totalWeights += newWeight;
          }
        }
        if (totalWeights != 0) {
          temprgba[color] /= totalWeights;
        }
      }
      Pixel(i,j).Reset(temprgba[0],temprgba[1],temprgba[2],temprgba[3]);
      Pixel(i,j).Clamp();
    }
  }

  free(spaceWeightArray);
  free(rangeWeightArray);
  free(tempPixels1);
}


void R2Image::
blendOtherImageTranslated(R2Image * otherImage)
{
	// find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
	// compute the matching translation (pixel precision is OK), and blend the translated "otherImage"
	// into this image with a 50% opacity.

  double sigma = 2;

  int windowSide = 6*sigma+1;
  int searchAreaHalfWidth = 0.1*otherImage->Width()+0.5*windowSide;
  int searchAreaHalfHeight = 0.1*otherImage->Height()+0.5*windowSide;

  // Harris

  R2Image harris(*this);


  double alpha = 0.04;

  R2Image Img1(*this);
  R2Image Img2(*this);
  R2Image Img3(*this);

  Img1.SobelX();
  Img2.SobelY();
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Img3.Pixel(i,j) = Img1.Pixel(i,j) * Img2.Pixel(i,j);
      Img1.Pixel(i,j) *= Img1.Pixel(i,j);
      Img2.Pixel(i,j) *= Img2.Pixel(i,j);
    }
  }

  Img1.Blur(sigma);
  Img2.Blur(sigma);
  Img3.Blur(sigma);


  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      harris.Pixel(i,j) = Img1.Pixel(i,j)*Img2.Pixel(i,j) - Img3.Pixel(i,j)*Img3.Pixel(i,j) - alpha*((Img1.Pixel(i,j)+Img2.Pixel(i,j))*(Img1.Pixel(i,j)+Img2.Pixel(i,j)));
      harris.Pixel(i,j).Reset(harris.Pixel(i,j)[0],harris.Pixel(i,j)[1],harris.Pixel(i,j)[2],harris.Pixel(i,j)[3]);
      //harris.Pixel(i,j).Clamp();
      // don't clamp here!
    }
  }

  // feature detection: detect the highest n features
  bool *features;
  features = new bool[width*height];

  std::vector<Feature> featureVec;
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      features[i*height+j] = false;
      featureVec.push_back(Feature(i, j, harris.Pixel(i,j)));
    }
  }

  std::sort (featureVec.begin(), featureVec.end());

  int total = 150;
  int count = total;
  int minDistance = 10;
  int *chosenFeatures;
  chosenFeatures = new int[2*count];
  int *correspondingFeatures;
  correspondingFeatures = new int[2*count];
  while (count > 0 && !featureVec.empty()) {
    bool flag = false;
    Feature feature = featureVec.back();
    featureVec.pop_back();
    for (int i = feature.centerX-minDistance; i <= feature.centerX+minDistance; i++){
      for (int j = feature.centerY-minDistance; j <= feature.centerY+minDistance; j++){
        if (i >= 0 && i < width && j >= 0 && j < height) {
          if (features[i*height+j]) {
            flag = true;
          }
        }
        if (flag) {
          break;
        }
      }
      if (flag) {
        break;
      }
    }

    if (!flag && feature.centerX >= searchAreaHalfWidth && feature.centerX < width - searchAreaHalfWidth && feature.centerY >=  searchAreaHalfHeight && feature.centerY < height - searchAreaHalfHeight) {
      count --;
      features[feature.centerX*height+feature.centerY] = true;
      chosenFeatures[count*2] = feature.centerX;
      chosenFeatures[count*2+1] = feature.centerY;
      //printf("x: %d, y: %d, count: %d", chosenFeatures[count*2], chosenFeatures[count*2+1], count);

    }
  }

  for (int i = 0; i < total; i++) {
    int x = chosenFeatures[2*i];
    int y = chosenFeatures[2*i+1];
    std::vector<SearchSSD> searchVec;

    int searchXStart = x - searchAreaHalfWidth;
    int searchXEnd = x + searchAreaHalfWidth;
    int searchYStart = y - searchAreaHalfHeight;
    int searchYEnd = y + searchAreaHalfHeight;
    if (searchXStart < 0) {
      searchXStart = 0;
    }
    else if (searchXEnd >= otherImage->Width()) {
      searchXEnd = otherImage->Width() - 1;
    }
    if (searchYStart < 0) {
      searchYStart = 0;
    }
    else if (searchYEnd >= otherImage->Width()) {
      searchYEnd = otherImage->Height() - 1;
    }
    for (int m = searchXStart; m <= searchXEnd - windowSide; m++) {
      for (int n = searchYStart; n <= searchYEnd - windowSide; n++) {
        double ssd = 0;
        for (int a = 0; a < windowSide; a++) {
          for (int b = 0; b < windowSide; b++) {
            ssd += pow((Pixel(x-3*sigma+a,y-3*sigma+b)[0]+Pixel(x-3*sigma+a,y-3*sigma+b)[1]+Pixel(x-3*sigma+a,y-3*sigma+b)[2]-otherImage->Pixel(m+a,n+b)[0]-otherImage->Pixel(m+a,n+b)[1]-otherImage->Pixel(m+a,n+b)[2]),2);
          }
        }
        searchVec.push_back(SearchSSD(m+3*sigma, n+3*sigma, ssd));
      }
    }
    std::sort (searchVec.begin(), searchVec.end());
    correspondingFeatures[2*i] = searchVec.front().searchX;
    correspondingFeatures[2*i+1] = searchVec.front().searchY;

  }

  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      otherImage->Pixel(i,j).SetAlpha(0.5);
      Pixel(i,j) = Pixel(i,j)*otherImage->Pixel(i,j);
      Pixel(i,j).Clamp();
    }
  }

  // RANSAC
  int distThreshold = 10;
  int numThreshold = 75;

  std::vector<int> inlierSet;

  for (int k = 0; k < numThreshold; k ++) {
    std::vector<int> tempInlierSet;
    int selected = rand() % total;
    int xSelected = chosenFeatures[2*selected];
    int ySelected = chosenFeatures[2*selected+1];
    int x2Selected = correspondingFeatures[2*selected];
    int y2Selected = correspondingFeatures[2*selected+1];
    int xSelectedMotion = x2Selected - xSelected;
    int ySelectedMotion = y2Selected - ySelected;

    for (int i = 0; i < total; i++) {
      int x = chosenFeatures[2*i];
      int y = chosenFeatures[2*i+1];
      int x2 = correspondingFeatures[2*i];
      int y2 = correspondingFeatures[2*i+1];
      int xMotion = x2 - x;
      int yMotion = y2 - y;
      int diff = (xMotion-xSelectedMotion)*(xMotion-xSelectedMotion) + (yMotion-ySelectedMotion)*(yMotion-ySelectedMotion);

      if (diff < distThreshold) {
        tempInlierSet.push_back(i);
      }

    }

    if (tempInlierSet.size() >= total*0.9) {
      inlierSet = tempInlierSet;
      break;
    }
    else if (inlierSet.size() < tempInlierSet.size()){
      inlierSet = tempInlierSet;
    }
  }

  for (int i = total-1; i >= 0; i--) {
    int x = chosenFeatures[2*i];
    int y = chosenFeatures[2*i+1];
    int x2 = correspondingFeatures[2*i];
    int y2 = correspondingFeatures[2*i+1];
    // printf("x:%d, y:%d, x2:%d, y2:%d", x,y,x2,y2);



    if (i == inlierSet.back()) {
      for (int m = -4; m <= 4; m++) {
        for (int n = -4; n <= 4; n++) {
          if (x+m >= 0 && x+m < width && y+n >= 0 && y+n < width && (m<-2 || m>2 || n<-2 || n>2)) {
            Pixel(x+m,y+n).Reset(0,1,0,0.5);
            Pixel(x+m,y+n).Clamp();
          }
          Pixel(x2+m,y2+n).Reset(0,1,0,0.5);
          Pixel(x2+m,y2+n).Clamp();

        }
      }
      inlierSet.pop_back();
    }
    else {
      for (int m = -4; m <= 4; m++) {
        for (int n = -4; n <= 4; n++) {
          if (x+m >= 0 && x+m < width && y+n >= 0 && y+n < width && (m<-2 || m>2 || n<-2 || n>2)) {
            Pixel(x+m,y+n).Reset(1,0,0,0.5);
            Pixel(x+m,y+n).Clamp();
          }
          Pixel(x2+m,y2+n).Reset(1,0,0,0.5);
          Pixel(x2+m,y2+n).Clamp();

        }
      }
    }

    // int A = 2*(y2-y);
    // int B = A-2*(x2-x);
    // int P = A-(x2-x);
    // int xStart = x2-x>0 ? x : x2;
    // int xEnd = x2-x>0 ? x2 : x;
    // int yTemp = y;

    int xStart;
    int xEnd;
    int yStart;
    int yEnd;
    int A;
    int B;
    int P;
    int yTemp;
    int xTemp;

    if (abs(x2-x)>=abs(y2-y)) {
      if(x2-x > 0) {
        xStart = x;
        xEnd = x2;
        yStart = y;
        yEnd = y2;
        A = 2*abs(y2-y);
        B = A-2*abs(x2-x);
        P = A-abs(x2-x);
        yTemp = y;
      }
      else {
        xStart = x2;
        xEnd = x;
        yStart = y2;
        yEnd = y;
        A = 2*abs(y-y2);
        B = A-2*abs(x-x2);
        P = A-abs(x-x2);
        yTemp = y2;
      }
      // printf("A:%d, B:%d, P:%d", A,B,P);
      for (xTemp = xStart; xTemp <= xEnd; xTemp ++) {
        if (P < 0) {
          Pixel(xTemp,yTemp).Reset(1,0,0,0.5);
          // Pixel(xTemp+1,yTemp).Reset(1,0,0,0.5);
          // Pixel(xTemp-1,yTemp).Reset(1,0,0,0.5);
          // Pixel(xTemp,yTemp+1).Reset(1,0,0,0.5);
          // Pixel(xTemp+1,yTemp+1).Reset(1,0,0,0.5);
          // Pixel(xTemp-1,yTemp+1).Reset(1,0,0,0.5);
          // Pixel(xTemp,yTemp-1).Reset(1,0,0,0.5);
          // Pixel(xTemp+1,yTemp-1).Reset(1,0,0,0.5);
          // Pixel(xTemp-1,yTemp-1).Reset(1,0,0,0.5);
          Pixel(xTemp,yTemp).Clamp();
          // Pixel(xTemp+1,yTemp).Clamp();
          // Pixel(xTemp-1,yTemp).Clamp();
          // Pixel(xTemp,yTemp+1).Clamp();
          // Pixel(xTemp+1,yTemp+1).Clamp();
          // Pixel(xTemp-1,yTemp+1).Clamp();
          // Pixel(xTemp,yTemp-1).Clamp();
          // Pixel(xTemp+1,yTemp-1).Clamp();
          // Pixel(xTemp-1,yTemp-1).Clamp();
          P += A;
          // printf("xtemp: %d, ytemp: %d", xTemp, yTemp);
        }
        else {
          if (yStart < yEnd) {
            Pixel(xTemp,++yTemp).Reset(1,0,0,0.5);
          }
          else {
            Pixel(xTemp,--yTemp).Reset(1,0,0,0.5);
          }

          // Pixel(xTemp+1,yTemp).Reset(1,0,0,0.5);
          // Pixel(xTemp-1,yTemp).Reset(1,0,0,0.5);
          // Pixel(xTemp,yTemp+1).Reset(1,0,0,0.5);
          // Pixel(xTemp+1,yTemp+1).Reset(1,0,0,0.5);
          // Pixel(xTemp-1,yTemp+1).Reset(1,0,0,0.5);
          // Pixel(xTemp,yTemp-1).Reset(1,0,0,0.5);
          // Pixel(xTemp+1,yTemp-1).Reset(1,0,0,0.5);
          // Pixel(xTemp-1,yTemp-1).Reset(1,0,0,0.5);
          Pixel(xTemp,yTemp).Clamp();
          // Pixel(xTemp+1,yTemp).Clamp();
          // Pixel(xTemp-1,yTemp).Clamp();
          // Pixel(xTemp,yTemp+1).Clamp();
          // Pixel(xTemp+1,yTemp+1).Clamp();
          // Pixel(xTemp-1,yTemp+1).Clamp();
          // Pixel(xTemp,yTemp-1).Clamp();
          // Pixel(xTemp+1,yTemp-1).Clamp();
          // Pixel(xTemp-1,yTemp-1).Clamp();
          P += B;
          // printf("xtemp: %d, ytemp: %d", xTemp, yTemp);
        }
      }
    }
    else {
      if(y2-y > 0) {
        yStart = y;
        yEnd = y2;
        xStart = x;
        xEnd = x2;
        A = 2*abs(x2-x);
        B = A-2*abs(y2-y);
        P = A-abs(y2-y);
        xTemp = x;
      }
      else {
        yStart = y2;
        yEnd = y;
        xStart = x2;
        xEnd = x;
        A = 2*abs(x-x2);
        B = A-2*abs(y-y2);
        P = A-abs(y-y2);
        xTemp = x2;
      }
      // printf("A:%d, B:%d, P:%d\n", A,B,P);
      // printf("yStart: %d, yEnd: %d, xStart: %d, xEnd: %d\n", yStart, yEnd, xStart, xEnd);
      // printf("A:%d, B:%d, P:%d", A,B,P);
      for (yTemp = yStart; yTemp <= yEnd; yTemp ++) {
        if (P < 0) {
          Pixel(xTemp,yTemp).Reset(1,0,0,0.5);
          // Pixel(xTemp+1,yTemp).Reset(1,0,0,0.5);
          // Pixel(xTemp-1,yTemp).Reset(1,0,0,0.5);
          // Pixel(xTemp,yTemp+1).Reset(1,0,0,0.5);
          // Pixel(xTemp+1,yTemp+1).Reset(1,0,0,0.5);
          // Pixel(xTemp-1,yTemp+1).Reset(1,0,0,0.5);
          // Pixel(xTemp,yTemp-1).Reset(1,0,0,0.5);
          // Pixel(xTemp+1,yTemp-1).Reset(1,0,0,0.5);
          // Pixel(xTemp-1,yTemp-1).Reset(1,0,0,0.5);
          Pixel(xTemp,yTemp).Clamp();
          // Pixel(xTemp+1,yTemp).Clamp();
          // Pixel(xTemp-1,yTemp).Clamp();
          // Pixel(xTemp,yTemp+1).Clamp();
          // Pixel(xTemp+1,yTemp+1).Clamp();
          // Pixel(xTemp-1,yTemp+1).Clamp();
          // Pixel(xTemp,yTemp-1).Clamp();
          // Pixel(xTemp+1,yTemp-1).Clamp();
          // Pixel(xTemp-1,yTemp-1).Clamp();
          P += A;
          // printf("xtemp: %d, ytemp: %d\n", xTemp, yTemp);
        }
        else {
          if (xStart < xEnd) {
            Pixel(++xTemp,yTemp).Reset(1,0,0,0.5);
          }
          else {
            Pixel(--xTemp,yTemp).Reset(1,0,0,0.5);
          }

          // Pixel(xTemp+1,yTemp).Reset(1,0,0,0.5);
          // Pixel(xTemp-1,yTemp).Reset(1,0,0,0.5);
          // Pixel(xTemp,yTemp+1).Reset(1,0,0,0.5);
          // Pixel(xTemp+1,yTemp+1).Reset(1,0,0,0.5);
          // Pixel(xTemp-1,yTemp+1).Reset(1,0,0,0.5);
          // Pixel(xTemp,yTemp-1).Reset(1,0,0,0.5);
          // Pixel(xTemp+1,yTemp-1).Reset(1,0,0,0.5);
          // Pixel(xTemp-1,yTemp-1).Reset(1,0,0,0.5);
          Pixel(xTemp,yTemp).Clamp();
          // Pixel(xTemp+1,yTemp).Clamp();
          // Pixel(xTemp-1,yTemp).Clamp();
          // Pixel(xTemp,yTemp+1).Clamp();
          // Pixel(xTemp+1,yTemp+1).Clamp();
          // Pixel(xTemp-1,yTemp+1).Clamp();
          // Pixel(xTemp,yTemp-1).Clamp();
          // Pixel(xTemp+1,yTemp-1).Clamp();
          // Pixel(xTemp-1,yTemp-1).Clamp();
          P += B;
          // printf("xtemp: %d, ytemp: %d", xTemp, yTemp);
        }
      }
    }


  }

  free(chosenFeatures);
  free(correspondingFeatures);
  free(features);

	return;
}

void R2Image::
blendOtherImageHomography(R2Image * otherImage)
{
  // find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
	// compute the matching translation (pixel precision is OK), and blend the translated "otherImage"
	// into this image with a 50% opacity.

  double sigma = 2;

  int windowSide = 6*sigma+1;
  int searchAreaHalfWidth = 0.1*otherImage->Width()+0.5*windowSide;
  int searchAreaHalfHeight = 0.1*otherImage->Height()+0.5*windowSide;

  // Harris

  R2Image harris(*this);


  double alpha = 0.04;

  R2Image Img1(*this);
  R2Image Img2(*this);
  R2Image Img3(*this);

  Img1.SobelX();
  Img2.SobelY();
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Img3.Pixel(i,j) = Img1.Pixel(i,j) * Img2.Pixel(i,j);
      Img1.Pixel(i,j) *= Img1.Pixel(i,j);
      Img2.Pixel(i,j) *= Img2.Pixel(i,j);
    }
  }

  Img1.Blur(sigma);
  Img2.Blur(sigma);
  Img3.Blur(sigma);


  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      harris.Pixel(i,j) = Img1.Pixel(i,j)*Img2.Pixel(i,j) - Img3.Pixel(i,j)*Img3.Pixel(i,j) - alpha*((Img1.Pixel(i,j)+Img2.Pixel(i,j))*(Img1.Pixel(i,j)+Img2.Pixel(i,j)));
      harris.Pixel(i,j).Reset(harris.Pixel(i,j)[0],harris.Pixel(i,j)[1],harris.Pixel(i,j)[2],harris.Pixel(i,j)[3]);
      //harris.Pixel(i,j).Clamp();
      // don't clamp here!
    }
  }

  // feature detection: detect the highest n features
  bool *features;
  features = new bool[width*height];

  std::vector<Feature> featureVec;
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      features[i*height+j] = false;
      featureVec.push_back(Feature(i, j, harris.Pixel(i,j)));
    }
  }

  std::sort (featureVec.begin(), featureVec.end());

  int total = 150;
  int count = total;
  int minDistance = 10;
  int *chosenFeatures;
  chosenFeatures = new int[2*count];
  int *correspondingFeatures;
  correspondingFeatures = new int[2*count];
  while (count > 0 && !featureVec.empty()) {
    bool flag = false;
    Feature feature = featureVec.back();
    featureVec.pop_back();
    for (int i = feature.centerX-minDistance; i <= feature.centerX+minDistance; i++){
      for (int j = feature.centerY-minDistance; j <= feature.centerY+minDistance; j++){
        if (i >= 0 && i < width && j >= 0 && j < height) {
          if (features[i*height+j]) {
            flag = true;
          }
        }
        if (flag) {
          break;
        }
      }
      if (flag) {
        break;
      }
    }

    if (!flag && feature.centerX >= searchAreaHalfWidth && feature.centerX < width - searchAreaHalfWidth && feature.centerY >=  searchAreaHalfHeight && feature.centerY < height - searchAreaHalfHeight) {
      count --;
      features[feature.centerX*height+feature.centerY] = true;
      chosenFeatures[count*2] = feature.centerX;
      chosenFeatures[count*2+1] = feature.centerY;
      //printf("x: %d, y: %d, count: %d", chosenFeatures[count*2], chosenFeatures[count*2+1], count);

    }
  }

  for (int i = 0; i < total; i++) {
    int x = chosenFeatures[2*i];
    int y = chosenFeatures[2*i+1];
    std::vector<SearchSSD> searchVec;

    int searchXStart = x - searchAreaHalfWidth;
    int searchXEnd = x + searchAreaHalfWidth;
    int searchYStart = y - searchAreaHalfHeight;
    int searchYEnd = y + searchAreaHalfHeight;
    if (searchXStart < 0) {
      searchXStart = 0;
    }
    else if (searchXEnd >= otherImage->Width()) {
      searchXEnd = otherImage->Width() - 1;
    }
    if (searchYStart < 0) {
      searchYStart = 0;
    }
    else if (searchYEnd >= otherImage->Width()) {
      searchYEnd = otherImage->Height() - 1;
    }
    for (int m = searchXStart; m <= searchXEnd - windowSide; m++) {
      for (int n = searchYStart; n <= searchYEnd - windowSide; n++) {
        double ssd = 0;
        for (int a = 0; a < windowSide; a++) {
          for (int b = 0; b < windowSide; b++) {
            ssd += pow((Pixel(x-3*sigma+a,y-3*sigma+b)[0]+Pixel(x-3*sigma+a,y-3*sigma+b)[1]+Pixel(x-3*sigma+a,y-3*sigma+b)[2]-otherImage->Pixel(m+a,n+b)[0]-otherImage->Pixel(m+a,n+b)[1]-otherImage->Pixel(m+a,n+b)[2]),2);
          }
        }
        searchVec.push_back(SearchSSD(m+3*sigma, n+3*sigma, ssd));
      }
    }
    std::sort (searchVec.begin(), searchVec.end());
    correspondingFeatures[2*i] = searchVec.front().searchX;
    correspondingFeatures[2*i+1] = searchVec.front().searchY;

  }

  // for (int i = 0; i < width; i++) {
  //   for (int j = 0;  j < height; j++) {
  //     otherImage->Pixel(i,j).SetAlpha(0.5);
  //     Pixel(i,j) = Pixel(i,j)*otherImage->Pixel(i,j);
  //     Pixel(i,j).Clamp();
  //   }
  // }

  // RANSAC
  int distThreshold = 3;
  int numThreshold = 1000;

  std::vector<int> inlierSet;
  double homographyMatrix[9];
  int point1;
  int point2;
  int point3;
  int point4;

  for (int k = 0; k < numThreshold; k ++) {
    printf("k=%d\n", k);
    std::vector<int> tempInlierSet;

    // DLT
    // Randomly select a small subset of points

    int idx1 = rand()%total;
    point1 = idx1;
    int idx2 = rand()%total;
    point2 = idx2;
    while (idx2 == idx1) {
      idx2 = rand()%total;
    }
    int idx3 = rand()%total;
    point3 = idx3;
    while (idx3 == idx2 || idx3 == idx1) {
      idx3 = rand()%total;
    }
    int idx4 = rand()%total;
    point4 = idx4;
    while (idx4 == idx3 || idx4 == idx2 || idx4 == idx1) {
      idx4 = rand()%total;
    }

    R2Point p1(chosenFeatures[2*idx1],chosenFeatures[2*idx1+1]);
    R2Point p2(chosenFeatures[2*idx2],chosenFeatures[2*idx2+1]);
    R2Point p3(chosenFeatures[2*idx3],chosenFeatures[2*idx3+1]);
    R2Point p4(chosenFeatures[2*idx4],chosenFeatures[2*idx4+1]);
    R2Point p1p(correspondingFeatures[2*idx1],correspondingFeatures[2*idx1+1]);
    R2Point p2p(correspondingFeatures[2*idx2],correspondingFeatures[2*idx2+1]);
    R2Point p3p(correspondingFeatures[2*idx3],correspondingFeatures[2*idx3+1]);
    R2Point p4p(correspondingFeatures[2*idx4],correspondingFeatures[2*idx4+1]);
    //printf("1\n");

    // build the 8x9 matrix of equations
    double** linEquations = dmatrix(1,8,1,9);

    linEquations[1][1] = 0;
    linEquations[1][2] = 0;
    linEquations[1][3] = 0;
    linEquations[1][4] = -p1[0];
    linEquations[1][5] = -p1[1];
    linEquations[1][6] = -1.0;
    linEquations[1][7] = p1p[1]*p1[0];
    linEquations[1][8] = p1p[1]*p1[1];
    linEquations[1][9] = p1p[1];

    linEquations[2][1] = p1[0];
    linEquations[2][2] = p1[1];
    linEquations[2][3] = 1.0;
    linEquations[2][4] = 0;
    linEquations[2][5] = 0;
    linEquations[2][6] = 0;
    linEquations[2][7] = -p1p[0]*p1[0];
    linEquations[2][8] = -p1p[0]*p1[1];
    linEquations[2][9] = -p1p[0];

    linEquations[3][1] = 0;
    linEquations[3][2] = 0;
    linEquations[3][3] = 0;
    linEquations[3][4] = -p2[0];
    linEquations[3][5] = -p2[1];
    linEquations[3][6] = -1.0;
    linEquations[3][7] = p2p[1]*p2[0];
    linEquations[3][8] = p2p[1]*p2[1];
    linEquations[3][9] = p2p[1];

    linEquations[4][1] = p2[0];
    linEquations[4][2] = p2[1];
    linEquations[4][3] = 1.0;
    linEquations[4][4] = 0;
    linEquations[4][5] = 0;
    linEquations[4][6] = 0;
    linEquations[4][7] = -p2p[0]*p2[0];
    linEquations[4][8] = -p2p[0]*p2[1];
    linEquations[4][9] = -p2p[0];

    linEquations[5][1] = 0;
    linEquations[5][2] = 0;
    linEquations[5][3] = 0;
    linEquations[5][4] = -p3[0];
    linEquations[5][5] = -p3[1];
    linEquations[5][6] = -1.0;
    linEquations[5][7] = p3p[1]*p3[0];
    linEquations[5][8] = p3p[1]*p3[1];
    linEquations[5][9] = p3p[1];

    linEquations[6][1] = p3[0];
    linEquations[6][2] = p3[1];
    linEquations[6][3] = 1.0;
    linEquations[6][4] = 0;
    linEquations[6][5] = 0;
    linEquations[6][6] = 0;
    linEquations[6][7] = -p3p[0]*p3[0];
    linEquations[6][8] = -p3p[0]*p3[1];
    linEquations[6][9] = -p3p[0];

    linEquations[7][1] = 0;
    linEquations[7][2] = 0;
    linEquations[7][3] = 0;
    linEquations[7][4] = -p4[0];
    linEquations[7][5] = -p4[1];
    linEquations[7][6] = -1.0;
    linEquations[7][7] = p4p[1]*p4[0];
    linEquations[7][8] = p4p[1]*p4[1];
    linEquations[7][9] = p4p[1];

    linEquations[8][1] = p4[0];
    linEquations[8][2] = p4[1];
    linEquations[8][3] = 1.0;
    linEquations[8][4] = 0;
    linEquations[8][5] = 0;
    linEquations[8][6] = 0;
    linEquations[8][7] = -p4p[0]*p4[0];
    linEquations[8][8] = -p4p[0]*p4[1];
    linEquations[8][9] = -p4p[0];
    //printf("2\n");

    // compute the SVD
    double singularValues[10]; // 1..9
    double** nullspaceMatrix = dmatrix(1,9,1,9);
    svdcmp(linEquations, 8, 9, singularValues, nullspaceMatrix);

    // find the smallest singular value:
    int smallestIndex = 1;
    for(int i=2;i<10;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;
    //printf("4\n");

    homographyMatrix[0] = nullspaceMatrix[1][smallestIndex];
    homographyMatrix[1] = nullspaceMatrix[2][smallestIndex];
    homographyMatrix[2] = nullspaceMatrix[3][smallestIndex];
    homographyMatrix[3] = nullspaceMatrix[4][smallestIndex];
    homographyMatrix[4] = nullspaceMatrix[5][smallestIndex];
    homographyMatrix[5] = nullspaceMatrix[6][smallestIndex];
    homographyMatrix[6] = nullspaceMatrix[7][smallestIndex];
    homographyMatrix[7] = nullspaceMatrix[8][smallestIndex];
    homographyMatrix[8] = nullspaceMatrix[9][smallestIndex];

    for (int i = 0; i < total; i++) {
      int x = chosenFeatures[2*i];
      int y = chosenFeatures[2*i+1];
      int x2 = correspondingFeatures[2*i];
      int y2 = correspondingFeatures[2*i+1];
      double newx2 = (nullspaceMatrix[1][smallestIndex]*x+nullspaceMatrix[2][smallestIndex]*y+nullspaceMatrix[3][smallestIndex]);
      double newy2 = (nullspaceMatrix[4][smallestIndex]*x+nullspaceMatrix[5][smallestIndex]*y+nullspaceMatrix[6][smallestIndex]);
      double newz2 = (nullspaceMatrix[7][smallestIndex]*x+nullspaceMatrix[8][smallestIndex]*y+nullspaceMatrix[9][smallestIndex]);
      double xp = newx2/newz2;
      double yp = newy2/newz2;

      // int xMotion = x2 - x;
      // int yMotion = y2 - y;
      double diff = sqrt((x2-xp)*(x2-xp) + (y2-yp)*(y2-yp));
      // printf("diff %d\n", diff);

      if (diff < distThreshold) {
        tempInlierSet.push_back(i);
      }

    }

    printf("inlierSet size: %lu\n", inlierSet.size());
    //printf("6\n");
    //printf("tempInlierSet %lu\n", tempInlierSet.size());

    // if (tempInlierSet.size() >= total*0.9) {
    //   //printf("7\n");
    //   inlierSet = tempInlierSet;
    //   break;
    // }
    // else
    if (inlierSet.size() < tempInlierSet.size()){
      //printf("8\n");
      inlierSet = tempInlierSet;
    }
  }
  //printf("hello\n");

  double** linEquations2 = dmatrix(1,inlierSet.size()*2,1,9);
  // printf("size of inlierset: %lu", inlierSet.size());
  int iter = 0;

  for (std::vector<int>::iterator it = inlierSet.begin(); it != inlierSet.end(); ++it) {
    // printf("iter: %d\n", *it);
    int index = *it;
    int x1temp = chosenFeatures[2*index];
    int y1temp = chosenFeatures[2*index+1];
    int x2temp = correspondingFeatures[2*index];
    int y2temp = correspondingFeatures[2*index+1];
    // printf("x1temp %d\n", x1temp);
    // printf("y1temp %d\n", y1temp);
    // printf("x2temp %d\n", x2temp);
    // printf("y2temp %d\n", y1temp);

    linEquations2[2*iter+1][1] = 0;
    linEquations2[2*iter+1][2] = 0;
    linEquations2[2*iter+1][3] = 0;
    linEquations2[2*iter+1][4] = -x1temp;
    linEquations2[2*iter+1][5] = -y1temp;
    linEquations2[2*iter+1][6] = -1.0;
    linEquations2[2*iter+1][7] = y2temp*x1temp;
    linEquations2[2*iter+1][8] = y2temp*y1temp;
    linEquations2[2*iter+1][9] = y2temp;
    // printf("hellohello\n");

    linEquations2[2*iter+2][1] = x1temp;
    linEquations2[2*iter+2][2] = y1temp;
    linEquations2[2*iter+2][3] = 1.0;
    linEquations2[2*iter+2][4] = 0;
    linEquations2[2*iter+2][5] = 0;
    linEquations2[2*iter+2][6] = 0;
    linEquations2[2*iter+2][7] = -x2temp*x1temp;
    linEquations2[2*iter+2][8] = -x2temp*y1temp;
    linEquations2[2*iter+2][9] = -x2temp;
    printf("byebye\n");

    iter ++;
  }

  // compute the SVD
  double singularValues2[10]; // 1..9
  double** nullspaceMatrix2 = dmatrix(1,inlierSet.size()*2+1,1,9);
  svdcmp(linEquations2, inlierSet.size()*2, 9, singularValues2, nullspaceMatrix2);

  // find the smallest singular value:
  int smallestIndex2 = 1;
  for(int i=2;i<10;i++) if(singularValues2[i]<singularValues2[smallestIndex2]) smallestIndex2=i;
  //printf("4\n");

  homographyMatrix[0] = nullspaceMatrix2[1][smallestIndex2];
  homographyMatrix[1] = nullspaceMatrix2[2][smallestIndex2];
  homographyMatrix[2] = nullspaceMatrix2[3][smallestIndex2];
  homographyMatrix[3] = nullspaceMatrix2[4][smallestIndex2];
  homographyMatrix[4] = nullspaceMatrix2[5][smallestIndex2];
  homographyMatrix[5] = nullspaceMatrix2[6][smallestIndex2];
  homographyMatrix[6] = nullspaceMatrix2[7][smallestIndex2];
  homographyMatrix[7] = nullspaceMatrix2[8][smallestIndex2];
  homographyMatrix[8] = nullspaceMatrix2[9][smallestIndex2];


  for (int i = 0; i < width; i ++) {
    for (int j = 0; j < height; j ++) {
      int x = i;
      int y = j;
      double newx2 = (homographyMatrix[0]*x+homographyMatrix[1]*y+homographyMatrix[2]);
      double newy2 = (homographyMatrix[3]*x+homographyMatrix[4]*y+homographyMatrix[5]);
      double newz2 = (homographyMatrix[6]*x+homographyMatrix[7]*y+homographyMatrix[8]);
      double xp = newx2/newz2;
      double yp = newy2/newz2;
      if (xp >= 0 && xp < width && yp >= 0 && yp <= height){
        Pixel(x,y) = Pixel(x,y)*0.5+otherImage->Pixel((int)xp,(int)yp)*0.5;
        //Pixel(x,y) = otherImage->Pixel((int)xp,(int)yp);
        printf("1111: x= %d, y=%d, xp= %d, yp= %d\n", (int)x,(int)y,(int)xp,(int)yp);
      }
      else {
        Pixel(x,y).Reset(0,0,0,0);
        printf("2222: x= %d, y=%d, xp= %d, yp= %d\n", (int)x,(int)y,(int)xp,(int)yp);
      }
    }
  }

  // for (int i = total-1; i >= 0; i--) {
  //   int x = chosenFeatures[2*i];
  //   int y = chosenFeatures[2*i+1];
  //   double newx2 = (homographyMatrix[0]*x+homographyMatrix[1]*y+homographyMatrix[2]);
  //   double newy2 = (homographyMatrix[3]*x+homographyMatrix[4]*y+homographyMatrix[5]);
  //   double newz2 = (homographyMatrix[6]*x+homographyMatrix[7]*y+homographyMatrix[8]);
  //   double xp = newx2/newz2;
  //   double yp = newy2/newz2;
  //   if (xp >= 0 && xp < width && yp >= 0 && yp <= height){
  //     //Pixel(x,y) = Pixel(x,y)*0.5+otherImage->Pixel(xp,yp)*0.5;
  //     Pixel(x,y) = otherImage->Pixel((int)xp,(int)yp);
  //     printf("1111: x= %d, y=%d, xp= %d, yp= %d\n", (int)x,(int)y,(int)xp,(int)yp);
  //   }
  //   else {
  //     Pixel(x,y).Reset(0,0,0,0);
  //     printf("2222: x= %d, y=%d, xp= %d, yp= %d\n", (int)x,(int)y,(int)xp,(int)yp);
  //   }
  //   printf("hahah");
  //
  // }

  // for (int i = total-1; i >= 0; i--) {
  //   printf("1\n");
  //   printf("i %d\n", i);
  //   int x = chosenFeatures[2*i];
  //   printf("x %d\n", x);
  //   int y = chosenFeatures[2*i+1];
  //   printf("y %d\n", y);
  //   int x2 = correspondingFeatures[2*i];
  //   printf("x2 %d\n", x2);
  //   int y2 = correspondingFeatures[2*i+1];
  //   printf("y2 %d\n", y2);
  //
  //   printf("hi\n");
  //   printf("inlierSet size %lu\n", inlierSet.size());
  //   printf("inlierSet back %d\n", inlierSet.back());
  //   bool correct = i == inlierSet.back();
  //   if (correct) {
  //     printf("2\n");
  //     for (int m = -4; m <= 4; m++) {
  //       for (int n = -4; n <= 4; n++) {
  //         if (x+m >= 0 && x+m < width && y+n >= 0 && y+n < width && (m<-2 || m>2 || n<-2 || n>2)) {
  //           Pixel(x+m,y+n).Reset(0,1,0,0.5);
  //           Pixel(x+m,y+n).Clamp();
  //         }
  //         Pixel(x2+m,y2+n).Reset(0,1,0,0.5);
  //         Pixel(x2+m,y2+n).Clamp();
  //
  //       }
  //     }
  //     inlierSet.pop_back();
  //   }
  //   else {
  //     printf("3\n");
  //     for (int m = -4; m <= 4; m++) {
  //       for (int n = -4; n <= 4; n++) {
  //         if (x+m >= 0 && x+m < width && y+n >= 0 && y+n < width && (m<-2 || m>2 || n<-2 || n>2)) {
  //           Pixel(x+m,y+n).Reset(1,0,0,0.5);
  //           Pixel(x+m,y+n).Clamp();
  //         }
  //         Pixel(x2+m,y2+n).Reset(1,0,0,0.5);
  //         Pixel(x2+m,y2+n).Clamp();
  //
  //       }
  //     }
  //   }
  //   printf("4\n");
  //
  //   int xStart;
  //   int xEnd;
  //   int yStart;
  //   int yEnd;
  //   int A;
  //   int B;
  //   int P;
  //   int yTemp;
  //   int xTemp;
  //
  //   if (abs(x2-x)>=abs(y2-y)) {
  //     if(x2-x > 0) {
  //       xStart = x;
  //       xEnd = x2;
  //       yStart = y;
  //       yEnd = y2;
  //       A = 2*abs(y2-y);
  //       B = A-2*abs(x2-x);
  //       P = A-abs(x2-x);
  //       yTemp = y;
  //     }
  //     else {
  //       xStart = x2;
  //       xEnd = x;
  //       yStart = y2;
  //       yEnd = y;
  //       A = 2*abs(y-y2);
  //       B = A-2*abs(x-x2);
  //       P = A-abs(x-x2);
  //       yTemp = y2;
  //     }
  //
  //     for (xTemp = xStart; xTemp <= xEnd; xTemp ++) {
  //       if (P < 0) {
  //         if (correct) {
  //           Pixel(xTemp,yTemp).Reset(0,1,0,0.5);
  //         }
  //         else {
  //           Pixel(xTemp,yTemp).Reset(1,0,0,0.5);
  //         }
  //         Pixel(xTemp,yTemp).Clamp();
  //         P += A;
  //       }
  //       else {
  //         if (yStart < yEnd) {
  //           if (correct) {
  //             Pixel(xTemp,yTemp).Reset(0,1,0,0.5);
  //           }
  //           else {
  //             Pixel(xTemp,yTemp).Reset(1,0,0,0.5);
  //           }
  //         }
  //         else {
  //           if (correct) {
  //             Pixel(xTemp,yTemp).Reset(0,1,0,0.5);
  //           }
  //           else {
  //             Pixel(xTemp,yTemp).Reset(1,0,0,0.5);
  //           }
  //         }
  //
  //         Pixel(xTemp,yTemp).Clamp();
  //         P += B;
  //       }
  //     }
  //   }
  //   else {
  //     if(y2-y > 0) {
  //       yStart = y;
  //       yEnd = y2;
  //       xStart = x;
  //       xEnd = x2;
  //       A = 2*abs(x2-x);
  //       B = A-2*abs(y2-y);
  //       P = A-abs(y2-y);
  //       xTemp = x;
  //     }
  //     else {
  //       yStart = y2;
  //       yEnd = y;
  //       xStart = x2;
  //       xEnd = x;
  //       A = 2*abs(x-x2);
  //       B = A-2*abs(y-y2);
  //       P = A-abs(y-y2);
  //       xTemp = x2;
  //     }
  //
  //     for (yTemp = yStart; yTemp <= yEnd; yTemp ++) {
  //       if (P < 0) {
  //         if (correct) {
  //           Pixel(xTemp,yTemp).Reset(0,1,0,0.5);
  //         }
  //         else {
  //           Pixel(xTemp,yTemp).Reset(1,0,0,0.5);
  //         }
  //         Pixel(xTemp,yTemp).Clamp();
  //         P += A;
  //       }
  //       else {
  //         if (xStart < xEnd) {
  //           if (correct) {
  //             Pixel(xTemp,yTemp).Reset(0,1,0,0.5);
  //           }
  //           else {
  //             Pixel(xTemp,yTemp).Reset(1,0,0,0.5);
  //           }
  //         }
  //         else {
  //           if (correct) {
  //             Pixel(xTemp,yTemp).Reset(0,1,0,0.5);
  //           }
  //           else {
  //             Pixel(xTemp,yTemp).Reset(1,0,0,0.5);
  //           }
  //         }
  //
  //         Pixel(xTemp,yTemp).Clamp();
  //         P += B;
  //       }
  //     }
  //   }
  //
  //
  // }

  free(chosenFeatures);
  free(correspondingFeatures);
  free(features);

	fprintf(stderr, "fit other image using a homography and blend imageB over imageA\n");
	return;
}


void R2Image::test() {
  // // fit a 2D conic to five points
  // R2Point p1(1.2,3.5);
  // R2Point p2(2.1,2.2);
  // R2Point p3(0.2,1.6);
  // R2Point p4(0.0,0.5);
  // R2Point p5(-0.2,4.2);

  // A
  {
  R2Point p1(0.0,0.0);
  R2Point p2(1.0,0.0);
  R2Point p3(1.0,1.0);
  R2Point p4(0.0,1.0);
  R2Point p1p(1.0,0.0);
  R2Point p2p(2.0,0.0);
  R2Point p3p(2.0,1.0);
  R2Point p4p(1.0,1.0);

  // build the 5x6 matrix of equations
  double** linEquations = dmatrix(1,8,1,9);

  linEquations[1][1] = 0;
  linEquations[1][2] = 0;
  linEquations[1][3] = 0;
  linEquations[1][4] = -p1[0];
  linEquations[1][5] = -p1[1];
  linEquations[1][6] = -1.0;
  linEquations[1][7] = p1p[1]*p1[0];
  linEquations[1][8] = p1p[1]*p1[1];
  linEquations[1][9] = p1p[1];

  linEquations[2][1] = p1[0];
  linEquations[2][2] = p1[1];
  linEquations[2][3] = 1.0;
  linEquations[2][4] = 0;
  linEquations[2][5] = 0;
  linEquations[2][6] = 0;
  linEquations[2][7] = -p1p[0]*p1[0];
  linEquations[2][8] = -p1p[0]*p1[1];
  linEquations[2][9] = -p1p[0];

  linEquations[3][1] = 0;
  linEquations[3][2] = 0;
  linEquations[3][3] = 0;
  linEquations[3][4] = -p2[0];
  linEquations[3][5] = -p2[1];
  linEquations[3][6] = -1.0;
  linEquations[3][7] = p2p[1]*p2[0];
  linEquations[3][8] = p2p[1]*p2[1];
  linEquations[3][9] = p2p[1];

  linEquations[4][1] = p2[0];
  linEquations[4][2] = p2[1];
  linEquations[4][3] = 1.0;
  linEquations[4][4] = 0;
  linEquations[4][5] = 0;
  linEquations[4][6] = 0;
  linEquations[4][7] = -p2p[0]*p2[0];
  linEquations[4][8] = -p2p[0]*p2[1];
  linEquations[4][9] = -p2p[0];

  linEquations[5][1] = 0;
  linEquations[5][2] = 0;
  linEquations[5][3] = 0;
  linEquations[5][4] = -p3[0];
  linEquations[5][5] = -p3[1];
  linEquations[5][6] = -1.0;
  linEquations[5][7] = p3p[1]*p3[0];
  linEquations[5][8] = p3p[1]*p3[1];
  linEquations[5][9] = p3p[1];

  linEquations[6][1] = p3[0];
  linEquations[6][2] = p3[1];
  linEquations[6][3] = 1.0;
  linEquations[6][4] = 0;
  linEquations[6][5] = 0;
  linEquations[6][6] = 0;
  linEquations[6][7] = -p3p[0]*p3[0];
  linEquations[6][8] = -p3p[0]*p3[1];
  linEquations[6][9] = -p3p[0];

  linEquations[7][1] = 0;
  linEquations[7][2] = 0;
  linEquations[7][3] = 0;
  linEquations[7][4] = -p4[0];
  linEquations[7][5] = -p4[1];
  linEquations[7][6] = -1.0;
  linEquations[7][7] = p4p[1]*p4[0];
  linEquations[7][8] = p4p[1]*p4[1];
  linEquations[7][9] = p4p[1];

  linEquations[8][1] = p4[0];
  linEquations[8][2] = p4[1];
  linEquations[8][3] = 1.0;
  linEquations[8][4] = 0;
  linEquations[8][5] = 0;
  linEquations[8][6] = 0;
  linEquations[8][7] = -p4p[0]*p4[0];
  linEquations[8][8] = -p4p[0]*p4[1];
  linEquations[8][9] = -p4p[0];

  // compute the SVD
  double singularValues[10]; // 1..9
  double** nullspaceMatrix = dmatrix(1,9,1,9);
  svdcmp(linEquations, 8, 9, singularValues, nullspaceMatrix);

  // get the result
  printf("A\n");
  printf("\n Singular values: %f, %f, %f, %f, %f, %f, %f, %f, %f\n",singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6],singularValues[7],singularValues[8],singularValues[9]);

  // find the smallest singular value:
  int smallestIndex = 1;
  for(int i=2;i<10;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

  // solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
  printf("H: %f, %f, %f, %f, %f, %f, %f, %f, %f\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],nullspaceMatrix[3][smallestIndex],nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex],nullspaceMatrix[7][smallestIndex],nullspaceMatrix[8][smallestIndex],nullspaceMatrix[9][smallestIndex]);

  }
  {
  // B

  R2Point p1(0.0,0.0);
  R2Point p2(1.0,0.0);
  R2Point p3(1.0,1.0);
  R2Point p4(0.0,1.0);
  R2Point p1p(1.0,2.0);
  R2Point p2p(1.0,1.0);
  R2Point p3p(3.0,1.0);
  R2Point p4p(3.0,2.0);

  // build the 5x6 matrix of equations
  double** linEquations = dmatrix(1,8,1,9);

  linEquations[1][1] = 0;
  linEquations[1][2] = 0;
  linEquations[1][3] = 0;
  linEquations[1][4] = -p1[0];
  linEquations[1][5] = -p1[1];
  linEquations[1][6] = -1.0;
  linEquations[1][7] = p1p[1]*p1[0];
  linEquations[1][8] = p1p[1]*p1[1];
  linEquations[1][9] = p1p[1];

  linEquations[2][1] = p1[0];
  linEquations[2][2] = p1[1];
  linEquations[2][3] = 1.0;
  linEquations[2][4] = 0;
  linEquations[2][5] = 0;
  linEquations[2][6] = 0;
  linEquations[2][7] = -p1p[0]*p1[0];
  linEquations[2][8] = -p1p[0]*p1[1];
  linEquations[2][9] = -p1p[0];

  linEquations[3][1] = 0;
  linEquations[3][2] = 0;
  linEquations[3][3] = 0;
  linEquations[3][4] = -p2[0];
  linEquations[3][5] = -p2[1];
  linEquations[3][6] = -1.0;
  linEquations[3][7] = p2p[1]*p2[0];
  linEquations[3][8] = p2p[1]*p2[1];
  linEquations[3][9] = p2p[1];

  linEquations[4][1] = p2[0];
  linEquations[4][2] = p2[1];
  linEquations[4][3] = 1.0;
  linEquations[4][4] = 0;
  linEquations[4][5] = 0;
  linEquations[4][6] = 0;
  linEquations[4][7] = -p2p[0]*p2[0];
  linEquations[4][8] = -p2p[0]*p2[1];
  linEquations[4][9] = -p2p[0];

  linEquations[5][1] = 0;
  linEquations[5][2] = 0;
  linEquations[5][3] = 0;
  linEquations[5][4] = -p3[0];
  linEquations[5][5] = -p3[1];
  linEquations[5][6] = -1.0;
  linEquations[5][7] = p3p[1]*p3[0];
  linEquations[5][8] = p3p[1]*p3[1];
  linEquations[5][9] = p3p[1];

  linEquations[6][1] = p3[0];
  linEquations[6][2] = p3[1];
  linEquations[6][3] = 1.0;
  linEquations[6][4] = 0;
  linEquations[6][5] = 0;
  linEquations[6][6] = 0;
  linEquations[6][7] = -p3p[0]*p3[0];
  linEquations[6][8] = -p3p[0]*p3[1];
  linEquations[6][9] = -p3p[0];

  linEquations[7][1] = 0;
  linEquations[7][2] = 0;
  linEquations[7][3] = 0;
  linEquations[7][4] = -p4[0];
  linEquations[7][5] = -p4[1];
  linEquations[7][6] = -1.0;
  linEquations[7][7] = p4p[1]*p4[0];
  linEquations[7][8] = p4p[1]*p4[1];
  linEquations[7][9] = p4p[1];

  linEquations[8][1] = p4[0];
  linEquations[8][2] = p4[1];
  linEquations[8][3] = 1.0;
  linEquations[8][4] = 0;
  linEquations[8][5] = 0;
  linEquations[8][6] = 0;
  linEquations[8][7] = -p4p[0]*p4[0];
  linEquations[8][8] = -p4p[0]*p4[1];
  linEquations[8][9] = -p4p[0];

  // compute the SVD
  double singularValues[10]; // 1..9
  double** nullspaceMatrix = dmatrix(1,9,1,9);
  svdcmp(linEquations, 8, 9, singularValues, nullspaceMatrix);

  // get the result
  printf("B\n");
  printf("\n Singular values: %f, %f, %f, %f, %f, %f, %f, %f, %f\n",singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6],singularValues[7],singularValues[8],singularValues[9]);

  // find the smallest singular value:
  int smallestIndex = 1;
  for(int i=2;i<10;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

  // solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
  printf("H: %f, %f, %f, %f, %f, %f, %f, %f, %f\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],nullspaceMatrix[3][smallestIndex],nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex],nullspaceMatrix[7][smallestIndex],nullspaceMatrix[8][smallestIndex],nullspaceMatrix[9][smallestIndex]);
  }
  return;
}

////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
  // Initialize everything
  if (pixels) { delete [] pixels; pixels = NULL; }
  npixels = width = height = 0;

  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }

  // Read file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
  else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R2Image::
Write(const char *filename) const
{
  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }

  // Write file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
  else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
  unsigned short int bfType;
  unsigned int bfSize;
  unsigned short int bfReserved1;
  unsigned short int bfReserved2;
  unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
  unsigned int biSize;
  int biWidth;
  int biHeight;
  unsigned short int biPlanes;
  unsigned short int biBitCount;
  unsigned int biCompression;
  unsigned int biSizeImage;
  int biXPelsPerMeter;
  int biYPelsPerMeter;
  unsigned int biClrUsed;
  unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
  unsigned char rgbtBlue;
  unsigned char rgbtGreen;
  unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
  unsigned char rgbBlue;
  unsigned char rgbGreen;
  unsigned char rgbRed;
  unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
  // Read a unsigned short int from a file in little endian format
  unsigned short int lsb, msb;
  lsb = getc(fp);
  msb = getc(fp);
  return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
  // Write a unsigned short int to a file in little endian format
  unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp);
  unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
  // Read a unsigned int word from a file in little endian format
  unsigned int b1 = getc(fp);
  unsigned int b2 = getc(fp);
  unsigned int b3 = getc(fp);
  unsigned int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
  // Write a unsigned int to a file in little endian format
  unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
  unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
  // Read a int word from a file in little endian format
  int b1 = getc(fp);
  int b2 = getc(fp);
  int b3 = getc(fp);
  int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
  // Write a int to a file in little endian format
  char b1 = (x & 0x000000FF); putc(b1, fp);
  char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  /* Read file header */
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = WordReadLE(fp);
  bmfh.bfSize = DWordReadLE(fp);
  bmfh.bfReserved1 = WordReadLE(fp);
  bmfh.bfReserved2 = WordReadLE(fp);
  bmfh.bfOffBits = DWordReadLE(fp);

  /* Check file header */
  assert(bmfh.bfType == BMP_BF_TYPE);
  /* ignore bmfh.bfSize */
  /* ignore bmfh.bfReserved1 */
  /* ignore bmfh.bfReserved2 */
  assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);

  /* Read info header */
  BITMAPINFOHEADER bmih;
  bmih.biSize = DWordReadLE(fp);
  bmih.biWidth = LongReadLE(fp);
  bmih.biHeight = LongReadLE(fp);
  bmih.biPlanes = WordReadLE(fp);
  bmih.biBitCount = WordReadLE(fp);
  bmih.biCompression = DWordReadLE(fp);
  bmih.biSizeImage = DWordReadLE(fp);
  bmih.biXPelsPerMeter = LongReadLE(fp);
  bmih.biYPelsPerMeter = LongReadLE(fp);
  bmih.biClrUsed = DWordReadLE(fp);
  bmih.biClrImportant = DWordReadLE(fp);

  // Check info header
  assert(bmih.biSize == BMP_BI_SIZE);
  assert(bmih.biWidth > 0);
  assert(bmih.biHeight > 0);
  assert(bmih.biPlanes == 1);
  assert(bmih.biBitCount == 24);  /* RGB */
  assert(bmih.biCompression == BI_RGB);   /* RGB */
  int lineLength = bmih.biWidth * 3;  /* RGB */
  if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
  assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);

  // Assign width, height, and number of pixels
  width = bmih.biWidth;
  height = bmih.biHeight;
  npixels = width * height;

  // Allocate unsigned char buffer for reading pixels
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = bmih.biSizeImage;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Read buffer
  fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
  if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
    fprintf(stderr, "Error while reading BMP file %s", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Allocate pixels for image
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double b = (double) *(p++) / 255;
      double g = (double) *(p++) / 255;
      double r = (double) *(p++) / 255;
      R2Pixel pixel(r, g, b, 1);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Compute number of bytes in row
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Write file header
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = BMP_BF_TYPE;
  bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
  bmfh.bfReserved1 = 0;
  bmfh.bfReserved2 = 0;
  bmfh.bfOffBits = BMP_BF_OFF_BITS;
  WordWriteLE(bmfh.bfType, fp);
  DWordWriteLE(bmfh.bfSize, fp);
  WordWriteLE(bmfh.bfReserved1, fp);
  WordWriteLE(bmfh.bfReserved2, fp);
  DWordWriteLE(bmfh.bfOffBits, fp);

  // Write info header
  BITMAPINFOHEADER bmih;
  bmih.biSize = BMP_BI_SIZE;
  bmih.biWidth = width;
  bmih.biHeight = height;
  bmih.biPlanes = 1;
  bmih.biBitCount = 24;       /* RGB */
  bmih.biCompression = BI_RGB;    /* RGB */
  bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
  bmih.biXPelsPerMeter = 2925;
  bmih.biYPelsPerMeter = 2925;
  bmih.biClrUsed = 0;
  bmih.biClrImportant = 0;
  DWordWriteLE(bmih.biSize, fp);
  LongWriteLE(bmih.biWidth, fp);
  LongWriteLE(bmih.biHeight, fp);
  WordWriteLE(bmih.biPlanes, fp);
  WordWriteLE(bmih.biBitCount, fp);
  DWordWriteLE(bmih.biCompression, fp);
  DWordWriteLE(bmih.biSizeImage, fp);
  LongWriteLE(bmih.biXPelsPerMeter, fp);
  LongWriteLE(bmih.biYPelsPerMeter, fp);
  DWordWriteLE(bmih.biClrUsed, fp);
  DWordWriteLE(bmih.biClrImportant, fp);

  // Write image, swapping blue and red in each pixel
  int pad = rowsize - width * 3;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      double r = 255.0 * pixel.Red();
      double g = 255.0 * pixel.Green();
      double b = 255.0 * pixel.Blue();
      if (r >= 255) r = 255;
      if (g >= 255) g = 255;
      if (b >= 255) b = 255;
      fputc((unsigned char) b, fp);
      fputc((unsigned char) g, fp);
      fputc((unsigned char) r, fp);
    }

    // Pad row
    for (int i = 0; i < pad; i++) fputc(0, fp);
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Read PPM file magic identifier
  char buffer[128];
  if (!fgets(buffer, 128, fp)) {
    fprintf(stderr, "Unable to read magic id in PPM file");
    fclose(fp);
    return 0;
  }

  // skip comments
  int c = getc(fp);
  while (c == '#') {
    while (c != '\n') c = getc(fp);
    c = getc(fp);
  }
  ungetc(c, fp);

  // Read width and height
  if (fscanf(fp, "%d%d", &width, &height) != 2) {
    fprintf(stderr, "Unable to read width and height in PPM file");
    fclose(fp);
    return 0;
  }

  // Read max value
  double max_value;
  if (fscanf(fp, "%lf", &max_value) != 1) {
    fprintf(stderr, "Unable to read max_value in PPM file");
    fclose(fp);
    return 0;
  }

  // Allocate image pixels
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for PPM file");
    fclose(fp);
    return 0;
  }

  // Check if raw or ascii file
  if (!strcmp(buffer, "P6\n")) {
    // Read up to one character of whitespace (\n) after max_value
    int c = getc(fp);
    if (!isspace(c)) putc(c, fp);

    // Read raw image data
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
        double r = (double) getc(fp) / max_value;
        double g = (double) getc(fp) / max_value;
        double b = (double) getc(fp) / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }
  else {
    // Read asci image data
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
	// Read pixel values
	int red, green, blue;
	if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
	  fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
	  fclose(fp);
	  return 0;
	}

	// Assign pixel values
	double r = (double) red / max_value;
	double g = (double) green / max_value;
	double b = (double) blue / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
  // Check type
  if (ascii) {
    // Open file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }

    // Print PPM image file
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
        if (((i+1) % 4) == 0) fprintf(fp, "\n");
      }
      if ((width % 4) != 0) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // Close file
    fclose(fp);
  }
  else {
    // Open file
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }

    // Print PPM image file
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%c%c%c", r, g, b);
      }
    }

    // Close file
    fclose(fp);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
  extern "C" {
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
  };
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize decompression info
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, fp);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  // Remember image attributes
  width = cinfo.output_width;
  height = cinfo.output_height;
  npixels = width * height;
  int ncomponents = cinfo.output_components;

  // Allocate pixels for image
  pixels = new R2Pixel [ npixels ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Allocate unsigned char buffer for reading image
  int rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Read scan lines
  // First jpeg pixel is top-left, so read pixels in opposite scan-line order
  while (cinfo.output_scanline < cinfo.output_height) {
    int scanline = cinfo.output_height - cinfo.output_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_read_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  // Close file
  fclose(fp);

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double r, g, b, a;
      if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
        p++;
      }
      else if (ncomponents == 3) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 4) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = (double) *(p++) / 255;
      }
      else {
        fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
        return 0;
      }
      R2Pixel pixel(r, g, b, a);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}




int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize compression info
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);
  cinfo.image_width = width; 	/* image width and height, in pixels */
  cinfo.image_height = height;
  cinfo.input_components = 3;		/* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
  cinfo.dct_method = JDCT_ISLOW;
  jpeg_set_defaults(&cinfo);
  cinfo.optimize_coding = TRUE;
  jpeg_set_quality(&cinfo, 95, TRUE);
  jpeg_start_compress(&cinfo, TRUE);

  // Allocate unsigned char buffer for reading image
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Fill buffer with pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      int r = (int) (255 * pixel.Red());
      int g = (int) (255 * pixel.Green());
      int b = (int) (255 * pixel.Blue());
      if (r > 255) r = 255;
      if (g > 255) g = 255;
      if (b > 255) b = 255;
      *(p++) = r;
      *(p++) = g;
      *(p++) = b;
    }
  }



  // Output scan lines
  // First jpeg pixel is top-left, so write in opposite scan-line order
  while (cinfo.next_scanline < cinfo.image_height) {
    int scanline = cinfo.image_height - cinfo.next_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  // Close file
  fclose(fp);

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return number of bytes written
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}
