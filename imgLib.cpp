#pragma once
#include <cv.h>
#include <highgui.h>
#include <io.h>
#include "imgLib.h"
using namespace cv;


void ImgResize(IplImage* pImgSrc, IplImage* pImgDst)
{
	if(pImgSrc->nChannels != 3 || pImgDst->nChannels != 3)
		return;

	if(pImgSrc->width == pImgDst->width && pImgSrc->height == pImgDst->height)
	{
		cvCopy(pImgSrc, pImgDst);
		return;
	}

	int x, y;
	int* posX = new int[pImgDst->width];
	int* posY = new int[pImgDst->height];
	
	for(x = 0; x < pImgDst->width; x ++)
	{
		float p = (x + 0.5) * (pImgSrc->width - 1) / pImgDst->width;
		posX[x] = p + 0.5;
	}
	for(y = 0; y < pImgDst->height; y ++)
	{
		float p = (y + 0.5) * (pImgSrc->height - 1) / pImgDst->height;
		posY[y] = p + 0.5;
	}

	for(y = 0; y < pImgDst->height; y ++)
	{
		int pos = posY[y];
		unsigned char *pBufDst  = (unsigned char *)pImgDst->imageData + y * pImgDst->widthStep;
		unsigned char *pBufSrc0  = (unsigned char *)pImgSrc->imageData + (pos + 0) * pImgSrc->widthStep;
		unsigned char *pBufSrc1  = (unsigned char *)pImgSrc->imageData + (pos + 1) * pImgSrc->widthStep;
		for(x = 0; x < pImgDst->width; x ++)
		{
			int xOffset = 3 * posX[x];
			pBufDst[0] = pBufSrc0[xOffset + 0];
			pBufDst[1] = pBufSrc0[xOffset + 1];
			pBufDst[2] = pBufSrc0[xOffset + 2];
			pBufDst  += 3;
		}
	}
	delete [] posX;
	delete [] posY;
}




#define GAMMA_ZOOM 256
#define FTBL_SIZE 10240
#define FTBL_BITS 8
#include "RGB2LAB.h"


void RGB2LAB(IplImage *pImgSrc)
{
const int M[] = {
	FTBL_SIZE / GAMMA_ZOOM * 0.412453 / 0.95047, 
	FTBL_SIZE / GAMMA_ZOOM * 0.357580 / 0.95047, 
	FTBL_SIZE / GAMMA_ZOOM * 0.180423 / 0.95047,
	FTBL_SIZE / GAMMA_ZOOM * 0.212671, 
	FTBL_SIZE / GAMMA_ZOOM * 0.715160, 
	FTBL_SIZE / GAMMA_ZOOM * 0.072169,
	FTBL_SIZE / GAMMA_ZOOM * 0.019334 / 1.08883, 
	FTBL_SIZE / GAMMA_ZOOM * 0.119193 / 1.08883, 
	FTBL_SIZE / GAMMA_ZOOM * 0.950227 / 1.08883
};

	for(int y = 0; y <= pImgSrc->height - 1; y ++)
	{
		unsigned char *pBuf = (unsigned char *)pImgSrc->imageData +y * pImgSrc->widthStep;
		for(int x = 0; x <= pImgSrc->width - 1; x ++)
		{
			int B = gammaTbl[pBuf[0]];
			int G = gammaTbl[pBuf[1]];
			int R = gammaTbl[pBuf[2]];

			int FX = fTbl[M[0] * R + M[1] * G + M[2] * B];
			int FY = fTbl[M[3] * R + M[4] * G + M[5] * B];
			int FZ = fTbl[M[6] * R + M[7] * G + M[8] * B];

			pBuf[0] = (116 * FY >> FTBL_BITS);// - 16;		//LÊä³ö·¶Î§ÊÇ16-116
		    pBuf[1] = (500 * (FX - FY) >> FTBL_BITS) + 128;
		    pBuf[2] = (200 * (FY - FZ) >> FTBL_BITS) + 128;

			pBuf += 3;
		}

	}
}


