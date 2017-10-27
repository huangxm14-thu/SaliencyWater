#pragma once
#include <cv.h>
#include <highgui.h>
#include <stdio.h>
#include "imgLib.h"
using namespace cv;
using namespace std;

//please modify header file path of intel ipp library
#include "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2016.2.180/windows/ipp/include/ipp.h"

#define rWater 5
#define wCenter 0.5
#define wEnhance 20
#define wMorph  30
#define wSmooth 0.5
#define wMix 0.7


#define KERSIZE 3   /* Kernel size  */

/* Next two defines are created to simplify code reading and understanding */
#define EXIT_MAIN exitLine:                                  /* Label for Exit */
#define check_sts(st) if((st) != ippStsNoErr) goto exitLine; /* Go to Exit if IPP function returned status different from ippStsNoErr */

/* Results of ippMalloc() are not validated because IPP functions perform bad arguments check and will return an appropriate status  */


int morphSmooth(int nWidth, int nHeight, Mat_<float> &salImg)
{
	Scalar s = mean(salImg);
	int smooth_size = wMorph * sqrt(s.val[0]);

    IppStatus status = ippStsNoErr;
    IppiMorphAdvState *pState = NULL;
    IppiSize roiSize= { nWidth, nHeight };            /* image size */
    IppiSize maskSize={smooth_size,smooth_size};                        /* mask size */
    Ipp8u *pSrc = NULL, *pDst = NULL, *pImg = NULL; /*  initial and working image */
    Ipp8u *pBuf = NULL;
	Ipp8u *pMask= new Ipp8u[smooth_size * smooth_size];
	memset(pMask, 0x1, smooth_size * smooth_size);
    int srcStep = 0;                                /* srcStep, in bytes, through the source image */
    int size = 0, specSize = 0, bufferSize = 0;     /*  working buffers size */
	int x, y;

    pSrc = ippiMalloc_8u_C1(roiSize.width, roiSize.height, &srcStep);
	for(y = 0; y < nHeight; y ++)
	{
		float *pSal = (float *)salImg.ptr(y);
		Ipp8u *pSrcBuf = &pSrc[y * srcStep];
		for(x = 0; x < nWidth; x ++)
		{
			*pSrcBuf ++ = (*pSal ++) * 255.f;
		}
	}

    pDst = ippiMalloc_8u_C1(roiSize.width, roiSize.height, &srcStep);

    check_sts( status = ippiMorphReconstructGetBufferSize(roiSize, ipp8u, 1, &size) )

    pBuf = ippsMalloc_8u(size);

    check_sts( status = ippiMorphAdvGetSize_8u_C1R( roiSize, maskSize, &specSize, &bufferSize ) )

	IppiMorphState* pSpec = NULL;
	Ipp8u* pBuffer = NULL;
	IppiBorderType borderType= ippBorderRepl;
	Ipp16u borderValue = 0;
	status = ippiMorphologyBorderGetSize_8u_C1R( roiSize, maskSize, &specSize, &bufferSize );
	if (status != ippStsNoErr) return status;
	pSpec = (IppiMorphState*)ippsMalloc_8u(specSize);
	pBuffer = (Ipp8u*)ippsMalloc_8u(bufferSize);
	status = ippiMorphologyBorderInit_8u_C1R( roiSize, pMask, maskSize, pSpec, pBuffer );
	if (status != ippStsNoErr) {
		ippsFree(pBuffer);
		ippsFree(pSpec);
		return status;
	}
	status = ippiErodeBorder_8u_C1R( pSrc, srcStep, pDst, srcStep, roiSize, borderType, borderValue,
		pSpec, pBuffer);
    check_sts( status = ippiMorphReconstructDilate_8u_C1IR(pSrc, srcStep, pDst, srcStep, roiSize, pBuf, (IppiNorm)ippiNormL1) )
	status = ippiDilateBorder_8u_C1R( pDst, srcStep, pSrc, srcStep, roiSize, borderType, borderValue,
		pSpec, pBuffer);
    check_sts( status = ippiMorphReconstructErode_8u_C1IR(pDst, srcStep, pSrc, srcStep, roiSize, pBuf, (IppiNorm)ippiNormL1) )


	for(y = 0; y < nHeight; y ++)
	{
		float *pSal = (float *)salImg.ptr(y);
		Ipp8u *pSrcBuf = &pSrc[y * srcStep];
		for(x = 0; x < nWidth; x ++)
		{
			*pSal ++= (*pSrcBuf ++) / 255.f;
		}
	}

EXIT_MAIN
	delete [] pMask;
	ippsFree(pBuffer);
	ippsFree(pSpec);
    ippsFree(pBuf);
    ippiFree(pSrc);
    ippiFree(pDst);
	if(status != ippStsNoErr)
    printf("Exit status %d (%s)\n", (int)status, ippGetStatusString(status));


    return status;

}


int FrameDetect(IplImage* pImgSrc, int& frameWidth)
{
	int i, j, k;

	//frame margin detection
	int wFrame;
	for(i = 0; i < 10; i ++)
	{
		int sumNorm1[3] = {0, 0, 0};
		int sumNorm2[3] = {0, 0, 0};
		uchar *pBufUp = (uchar* )pImgSrc->imageData + i * pImgSrc->widthStep;
		uchar *pBufDn = (uchar* )pImgSrc->imageData + ( pImgSrc->height - 1 - i) * pImgSrc->widthStep;
		for( j = 0; j < pImgSrc->width; j++ )
		{
			sumNorm1[0] += pBufUp[0] + pBufDn[0];
			sumNorm1[1] += pBufUp[1] + pBufDn[1];
			sumNorm1[2] += pBufUp[2] + pBufDn[2];
			sumNorm2[0] += pBufUp[0] * pBufUp[0] + pBufDn[0] * pBufDn[0];
			sumNorm2[1] += pBufUp[1] * pBufUp[1] + pBufDn[1] * pBufDn[1];
			sumNorm2[2] += pBufUp[2] * pBufUp[2] + pBufDn[2] * pBufDn[2];

			pBufUp += 3;
			pBufDn += 3;
		}

		uchar *pBufLt = (uchar* )pImgSrc->imageData + i * pImgSrc->nChannels;
		uchar *pBufRt = (uchar* )pImgSrc->imageData + (pImgSrc->width - 1 - i) * pImgSrc->nChannels;
		for(j = 0; j < pImgSrc->height; j ++ )
		{
			sumNorm1[0] += pBufLt[0] + pBufRt[0];
			sumNorm1[1] += pBufLt[1] + pBufRt[1];
			sumNorm1[2] += pBufLt[2] + pBufRt[2];
			sumNorm2[0] += pBufLt[0] * pBufLt[0] + pBufRt[0] * pBufRt[0];
			sumNorm2[1] += pBufLt[1] * pBufLt[1] + pBufRt[1] * pBufRt[1];
			sumNorm2[2] += pBufLt[2] * pBufLt[2] + pBufRt[2] * pBufRt[2];

			pBufLt += pImgSrc->widthStep;
			pBufRt += pImgSrc->widthStep;
		}

		int N = 2 * (pImgSrc->width + pImgSrc->height);
		double diffSqrSum = 0;
		for(k = 0; k < 3; k ++)
			diffSqrSum += (double)2 * N * sumNorm2[k] - (double)2 * sumNorm1[k] * sumNorm1[k];
		double diffAver = sqrt(diffSqrSum / (N * N * 3));
		if(diffAver > 10)
			break;
	}

	frameWidth = (i > 0) ? (i + 4) : 0;

	return 0;
}


/****************************************************************************************\
*                                       Watershed                                        *
\****************************************************************************************/

typedef struct CvWSNode
{
    struct CvWSNode* next;
    int mask_ofs;
}
CvWSNode;


typedef struct CvWSQueue
{
    CvWSNode* first;
    CvWSNode* last;
}
CvWSQueue;

static CvWSNode*
icvAllocWSNodes( CvMemStorage* storage )
{
    CvWSNode* n = 0;

    int i, count = (storage->block_size - sizeof(CvMemBlock))/sizeof(*n) - 1;

    n = (CvWSNode*)cvMemStorageAlloc( storage, count*sizeof(*n) );
    for( i = 0; i < count-1; i++ )
        n[i].next = n + i + 1;
    n[count-1].next = 0;

    return n;
}


#define m_fg 0x40


void __cvWatershed( Mat* src, Mat* dst, Mat_<float> &salImg, int r)
{
    const int IN_QUEUE = 0x2;
    const int NQ = 256;
    cv::Ptr<CvMemStorage> storage;

    CvSize size;
    CvWSNode* free_node = 0, *node;
    CvWSQueue q[NQ];
    int active_queue;
    int i, j;
    int db, dg, dr;
    char* mask;
    uchar* img;
    int mstep, istep;
	uchar *bdTbl;

	if(src->cols & 0x3)	//need 4 bytes aligned
		return;
	if(src->channels() != 3)
		return;

    #define ws_push(idx,mofs)  \
    {                               \
        if( !free_node )            \
            free_node = icvAllocWSNodes( storage );\
        node = free_node;           \
        free_node = free_node->next;\
        node->next = 0;             \
        node->mask_ofs = mofs;      \
        if( q[idx].last )           \
            q[idx].last->next=node; \
        else                        \
            q[idx].first = node;    \
        q[idx].last = node;         \
    }

    #define ws_pop(idx,mofs)   \
    {                               \
        node = q[idx].first;        \
        q[idx].first = node->next;  \
        if( !node->next )           \
            q[idx].last = 0;        \
        node->next = free_node;     \
        free_node = node;           \
        mofs = node->mask_ofs;      \
    }

    #define c_diff(ptr1,ptr2,diff)      \
    {                                   \
        db = abs((ptr1)[0] - (ptr2)[0]);\
        dg = abs((ptr1)[1] - (ptr2)[1]);\
        dr = abs((ptr1)[2] - (ptr2)[2]);\
        diff = (db + dg + dr) / 3;         \
        assert( 0 <= diff && diff <= 255 ); \
    }

	#define MaxMinUpdate(pBdDst, pBdSrc, ptr, t)	\
	{\
		pBdDst[0] = max(pBdSrc[0], ptr[0]);\
		pBdDst[1] = min(pBdSrc[1], ptr[0]);\
		pBdDst[2] = max(pBdSrc[2], ptr[1]);\
		pBdDst[3] = min(pBdSrc[3], ptr[1]);\
		pBdDst[4] = max(pBdSrc[4], ptr[2]);\
		pBdDst[5] = min(pBdSrc[5], ptr[2]);\
		t = (pBdDst[0] - pBdDst[1] + pBdDst[2] - pBdDst[3] + pBdDst[4] - pBdDst[5]) / 3;	\
	}

	#define MaxMinInit(pBdDst, ptr, ptr2)	\
	{\
		pBdDst[0] = max(ptr[0], ptr2[0]);\
		pBdDst[1] = min(ptr[0], ptr2[0]);\
		pBdDst[2] = max(ptr[1], ptr2[1]);\
		pBdDst[3] = min(ptr[1], ptr2[1]);\
		pBdDst[4] = max(ptr[2], ptr2[2]);\
		pBdDst[5] = min(ptr[2], ptr2[2]);\
	}

	size = cvSize(src->cols, src->rows);
    storage = cvCreateMemStorage();

    istep = src->step;
	img = src->ptr(0);
    mstep = dst->step / sizeof(mask[0]);
    mask = (char *)dst->ptr(0);
	uchar bdTbl2[6];
	bdTbl = new uchar [size.height * mstep * 6];
	memset(bdTbl, 0x0, size.height * mstep * 6);

    memset( q, 0, NQ*sizeof(q[0]) );


	int t;
	uchar *ptr, *ptr2;
	char *m;

	//top
	i = r;
	j = r;
	ptr = src->ptr(i) +  j * 3;
	m = (char *)dst->ptr(i) + j;
	for(; j <= size.width - 1 - r; j ++ )
	{
		ptr2 = ptr - istep;
		c_diff( ptr, ptr2, t );
		salImg.at<float>(i, j) = t;
		ws_push(t, i*mstep + j);
		uchar *pBd = &bdTbl[(i * mstep + j) * 6];
		MaxMinInit(pBd, ptr, ptr2)

		m[0] = IN_QUEUE;
		m ++;
		ptr += 3;
	}

	//bottom
	i = size.height - 1 - r;
	j = r;
	ptr = src->ptr(i) +  j * 3;
	m = (char *)dst->ptr(i) + j;
	for(; j <= size.width - 1 - r; j ++ )
	{
		ptr2 = ptr + istep;
		c_diff( ptr, ptr2, t );
		salImg.at<float>(i, j) = t;
		ws_push(t, i*mstep + j);
		uchar *pBd = &bdTbl[(i * mstep + j) * 6];
		MaxMinInit(pBd, ptr, ptr2)

		m[0] = IN_QUEUE;
		m ++;
		ptr += 3;
	}

	//left
	i = r;
	j = r;
	ptr = src->ptr(i) +  j * 3;
	m = (char *)dst->ptr(i) + j;
	for(; i <= size.height - 1 - r; i ++ )
	{
		ptr2 = ptr - 3;
		c_diff( ptr, ptr2, t );
		salImg.at<float>(i, j) = t;
		ws_push(t, i*mstep + j);
		uchar *pBd = &bdTbl[(i * mstep + j) * 6];
		MaxMinInit(pBd, ptr, ptr2)

		m[0] = IN_QUEUE;
		m += mstep;
		ptr += istep;
	}

	//right
	i = r;
	j = size.width - 1 - r;
	ptr = src->ptr(i) +  j * 3;
	m = (char *)dst->ptr(i) + j;
	for(; i <= size.height - 1 - r; i ++ )
	{
		ptr2 = ptr + 3;
		c_diff( ptr, ptr2, t );
		salImg.at<float>(i, j) = t;
		ws_push(t, i*mstep + j);
		uchar *pBd = &bdTbl[(i * mstep + j) * 6];
		MaxMinInit(pBd, ptr, ptr2)

		m[0] = IN_QUEUE;
		m += mstep;
		ptr += istep;
	}


    // find the first non-empty queue
    for( i = 0; i < NQ; i++ )
        if( q[i].first )
            break;

    // if there is no markers, exit immediately
    if( i == NQ )
        return;

    active_queue = i;
    img = src->ptr(0);
    mask = (char *)dst->ptr(0);

    // recursively fill the basins
    for(;;)
    {
        int mofs;
        char* m;
        uchar* ptr;
		float* pSal;

        if( q[active_queue].first == 0 )
        {
            for( i = active_queue+1; i < NQ; i++ )
                if( q[i].first )
                    break;
            if( i == NQ )
                break;
            active_queue = i;
        }

        ws_pop( active_queue, mofs);

        m = mask + mofs;
        ptr = img + 3 * mofs;
		pSal = (float *)salImg.data + mofs;

		if(m[0] == m_fg)
			continue;
		m[0] = m_fg;


        if( m[-1] == 0 )
        {
			m[-1] = IN_QUEUE;
			uchar *pBdSrc = &bdTbl[mofs * 6];
			uchar *pBd = pBdSrc - 6;
			uchar* pCur = ptr - 3;
			MaxMinUpdate(pBd, pBdSrc, pCur, t);
            ws_push( t, mofs - 1);
			pSal[-1] = t;
        }
        else if( m[-1] == IN_QUEUE && pSal[-1] > pSal[0] + 1)
        {
			uchar *pBd = &bdTbl2[0];
			uchar *pBdSrc = &bdTbl[mofs * 6];
			uchar* pCur = ptr - 3;
			MaxMinUpdate(pBd, pBdSrc, pCur, t);
			if(t < pSal[-1])
			{
//				memcpy(&bdTbl[(mofs - 1) * 6], bdTbl2, 6);
				bdTbl[(mofs - 1) * 6 + 0] = bdTbl2[0];
				bdTbl[(mofs - 1) * 6 + 1] = bdTbl2[1];
				bdTbl[(mofs - 1) * 6 + 2] = bdTbl2[2];
				bdTbl[(mofs - 1) * 6 + 3] = bdTbl2[3];
				bdTbl[(mofs - 1) * 6 + 4] = bdTbl2[4];
				bdTbl[(mofs - 1) * 6 + 5] = bdTbl2[5];
			    ws_push( t, mofs - 1);
				pSal[-1] = t;
			}
        }
        if( m[1] == 0 )
        {
			m[1] = IN_QUEUE;
			uchar *pBdSrc = &bdTbl[mofs * 6];
			uchar *pBd = pBdSrc + 6;
			uchar* pCur = ptr + 3;
			MaxMinUpdate(pBd, pBdSrc, pCur, t);
            ws_push( t, mofs + 1);
			pSal[1] = t;
        }
        else if( m[1] == IN_QUEUE && pSal[1] > pSal[0] + 1)
        {
			uchar *pBd = &bdTbl2[0];
			uchar *pBdSrc = &bdTbl[mofs * 6];
			uchar* pCur = ptr + 3;
			MaxMinUpdate(pBd, pBdSrc, pCur, t);
			if(t < pSal[1])
			{
//				memcpy(&bdTbl[(mofs + 1) * 6], bdTbl2, 6);
				bdTbl[(mofs + 1) * 6 + 0] = bdTbl2[0];
				bdTbl[(mofs + 1) * 6 + 1] = bdTbl2[1];
				bdTbl[(mofs + 1) * 6 + 2] = bdTbl2[2];
				bdTbl[(mofs + 1) * 6 + 3] = bdTbl2[3];
				bdTbl[(mofs + 1) * 6 + 4] = bdTbl2[4];
				bdTbl[(mofs + 1) * 6 + 5] = bdTbl2[5];
			    ws_push( t, mofs + 1);
				pSal[1] = t;
			}
        }
        if( m[-mstep] == 0 )
        {
			m[-mstep] = IN_QUEUE;
			uchar *pBdSrc = &bdTbl[mofs * 6];
			uchar *pBd = pBdSrc  - mstep * 6;
			uchar* pCur = ptr - istep;
			MaxMinUpdate(pBd, pBdSrc, pCur, t);
            ws_push( t, mofs - mstep);
			pSal[-mstep] = t;
        }
        else if( m[-mstep] == IN_QUEUE && pSal[-mstep] > pSal[0] + 1)
        {
			uchar *pBd = &bdTbl2[0];
			uchar *pBdSrc = &bdTbl[mofs * 6];
			uchar* pCur = ptr - istep;
			MaxMinUpdate(pBd, pBdSrc, pCur, t);
			if(t < pSal[-mstep])
			{
//				memcpy(&bdTbl[(mofs - mstep) * 6], bdTbl2, 6);
				bdTbl[(mofs - mstep) * 6 + 0] = bdTbl2[0];
				bdTbl[(mofs - mstep) * 6 + 1] = bdTbl2[1];
				bdTbl[(mofs - mstep) * 6 + 2] = bdTbl2[2];
				bdTbl[(mofs - mstep) * 6 + 3] = bdTbl2[3];
				bdTbl[(mofs - mstep) * 6 + 4] = bdTbl2[4];
				bdTbl[(mofs - mstep) * 6 + 5] = bdTbl2[5];
			    ws_push( t, mofs - mstep);
				pSal[-mstep] = t;
			}
        }
        if( m[mstep] == 0 )
        {
			m[mstep] = IN_QUEUE;
			uchar *pBdSrc = &bdTbl[mofs * 6];
			uchar *pBd = pBdSrc  + mstep * 6;
			uchar* pCur = ptr + istep;
			MaxMinUpdate(pBd, pBdSrc, pCur, t);
            ws_push( t, mofs + mstep);
			pSal[mstep] = t;
        }
        else if( m[mstep] == IN_QUEUE && pSal[mstep] > pSal[0] + 1)
        {
			uchar *pBd = &bdTbl2[0];
			uchar *pBdSrc = &bdTbl[mofs * 6];
			uchar* pCur = ptr + istep;
			MaxMinUpdate(pBd, pBdSrc, pCur, t);
			if(t < pSal[mstep])
			{
//				memcpy(&bdTbl[(mofs + mstep) * 6], bdTbl2, 6);
				bdTbl[(mofs + mstep) * 6 + 0] = bdTbl2[0];
				bdTbl[(mofs + mstep) * 6 + 1] = bdTbl2[1];
				bdTbl[(mofs + mstep) * 6 + 2] = bdTbl2[2];
				bdTbl[(mofs + mstep) * 6 + 3] = bdTbl2[3];
				bdTbl[(mofs + mstep) * 6 + 4] = bdTbl2[4];
				bdTbl[(mofs + mstep) * 6 + 5] = bdTbl2[5];
			    ws_push( t, mofs + mstep);
				pSal[mstep] = t;
			}
        }
    }

	delete [] bdTbl;
}


void SaliencyWaterRaw(IplImage *pImgSrc, Mat_<float> &salImg, int frameWidth)
{
	int64 t1, t2;
	int i, j;
	int r = frameWidth + rWater;	//2;
	Mat markers_bg(pImgSrc->height, pImgSrc->width, CV_8U);
//	markers_bg = Scalar::all(0);

	int bndSum[3] = {0, 0, 0};
	for( i = frameWidth; i < r; i ++ )
	{
		j = frameWidth;
		uchar *pBufUp = (uchar* )pImgSrc->imageData + i * pImgSrc->widthStep + j * pImgSrc->nChannels;
		uchar *pBufDn = (uchar* )pImgSrc->imageData + ( markers_bg.rows - 1 - i) * pImgSrc->widthStep + j * pImgSrc->nChannels;
		for(; j < markers_bg.cols - frameWidth; j++ )
		{
			bndSum[0] += pBufUp[0] + pBufDn[0];
			bndSum[1] += pBufUp[1] + pBufDn[1];
			bndSum[2] += pBufUp[2] + pBufDn[2];
			pBufUp += 3;
			pBufDn += 3;
		}
	}

	for( i = r; i < markers_bg.rows - r; i++ )
	{
		j = frameWidth;
		uchar *pBufLt = (uchar* )pImgSrc->imageData + i * pImgSrc->widthStep + j * pImgSrc->nChannels;
		uchar *pBufRt = (uchar* )pImgSrc->imageData + i * pImgSrc->widthStep + (markers_bg.cols - 1 - j) * pImgSrc->nChannels;
		for(; j < r; j ++ )
		{
			bndSum[0] += pBufLt[0] + pBufRt[0];
			bndSum[1] += pBufLt[1] + pBufRt[1];
			bndSum[2] += pBufLt[2] + pBufRt[2];
			pBufLt += 3;
			pBufRt -= 3;
		}
	}
	int bndCnt = (markers_bg.cols - 2 * frameWidth + markers_bg.rows - 2 * r) * (r - frameWidth) * 2;
	bndSum[0] /= bndCnt;
	bndSum[1] /= bndCnt;
	bndSum[2] /= bndCnt;
	
	

	for( i = 0; i < markers_bg.rows; i++ )
	{
		uchar *pDataBg = markers_bg.data + i * markers_bg.step;
		float *pSal = (float *)salImg.ptr(i);
		for( j = 0; j < markers_bg.cols; j++ )
		{
			if(i < r || i >= markers_bg.rows - r || j < r || j >= markers_bg.cols - r)
			{
				uchar* pBuf = (uchar* )pImgSrc->imageData + i * pImgSrc->widthStep + j * pImgSrc->nChannels;

				*pDataBg = m_fg;
				*pSal = (1 - wMix) * (abs(bndSum[0] - pBuf[0]) + abs(bndSum[1] - pBuf[1]) + abs(bndSum[2] - pBuf[2])) / 3;
				pBuf[0] = wMix * pBuf[0] + (1 - wMix) * bndSum[0];
				pBuf[1] = wMix * pBuf[1] + (1 - wMix) * bndSum[1];
				pBuf[2] = wMix * pBuf[2] + (1 - wMix) * bndSum[2];
			}
			else
				*pDataBg = 0;
			pDataBg ++;
			pSal ++;
		}
	}

	Mat img(pImgSrc); 
    
	t1 = getTickCount();
	__cvWatershed( &img, &markers_bg, salImg, r);
	t2 = getTickCount();

}



void SaliencyCenterPrior(IplImage *pImgSrc, Mat_<float> &salImg)
{
	int i, j;
	int W = pImgSrc->width, H = pImgSrc->height;

	float sigma_space = wCenter * std::max( pImgSrc->width, pImgSrc->height ); 
	sigma_space = sigma_space * sigma_space;
	for( j=0; j<H; j++ )
	{
		int dy = (j - H / 2);
		dy = W * W / 4 + dy * dy;
		float *pSal = (float *)salImg.ptr(j);
		for( i=0; i<W; i++ )
		{
			(*pSal) *= exp(-(i * i - i * W + dy) / sigma_space);
			pSal ++;
		}
	}
}


#define N 24
int colorCntTbl[N*N*N];
float salSmoothTbl[N*N*N];

void SaliencySmooth(IplImage *pImgSrc, Mat_<float> &salImg)
{
	int r = 20;
	int x, y;
	int *pIdxMap = new int[pImgSrc->height * pImgSrc->width];


	memset(colorCntTbl, 0x0, sizeof(colorCntTbl));
	memset(salSmoothTbl, 0x0, sizeof(salSmoothTbl));
	for(y = 0; y < pImgSrc->height; y ++)
	{
		unsigned char *pBuf = (unsigned char *)&pImgSrc->imageData[y * pImgSrc->widthStep];
		float *pSal = (float *)salImg.ptr(y);
		for(x = 0; x < pImgSrc->width; x ++)
		{
			int a = pBuf[0] * N / 256;
			int b = pBuf[1] * N / 256;
			int c = pBuf[2] * N / 256;
			int idx = a * N * N + b * N + c;
			pIdxMap[y * pImgSrc->width + x] = idx;
			colorCntTbl[idx] ++;
			salSmoothTbl[idx] += *pSal;

			if(y < r/* || y > pImgSrc->height * (1 - r)*/ || x < r || x > pImgSrc->width - r)
			{
				colorCntTbl[idx] += 8;
				salSmoothTbl[idx] += (*pSal) * 8;
			}
			pBuf += 3;
			pSal ++;
		}
	}

	for(int i = 0; i < N * N * N; i ++)
	{
		if(colorCntTbl[i] > 0)
		{
			salSmoothTbl[i] /= colorCntTbl[i];
		}
		else
		{
			salSmoothTbl[i] = 0;
		}
	}


	for(y = 0; y < pImgSrc->height; y ++)
	{
		unsigned char *pBuf = (unsigned char *)&pImgSrc->imageData[y * pImgSrc->widthStep];
		float *pSal = (float *)salImg.ptr(y);
		for(x = 0; x < pImgSrc->width; x ++)
		{
			int idx = pIdxMap[y * pImgSrc->width + x];

			(*pSal) = (*pSal) * (1 - wSmooth) + salSmoothTbl[idx] * wSmooth;
			pBuf += 3;
			pSal ++;
		}
	}

	delete pIdxMap;
	return;
}
#undef N


float OtsuThre(int nWidth, int nHeight, Mat_<float> &salImg)
{
	float r = 0.03;
	int pixelCount[256], totalCount;  
    float pixelPro[256];  
    int i, j;
	int threshold = 0; 
 
    //初始化  
    for(i = 0; i < 256; i++)  
    {  
        pixelCount[i] = 0;  
        pixelPro[i] = 0;  
    }  
  
    //统计灰度级中每个像素在整幅图像中的个数  
	totalCount = 0;
    for(i = nHeight * r; i < nHeight * (1 - r); i++)
		for(j = nWidth * r; j < nWidth * (1 - r); j ++)
	    {  
			pixelCount[(int)(salImg.at<float>(i,j) * 255)] ++;  
			totalCount ++;
	    }  
  
    //计算每个像素在整幅图像中的比例
	float utmp = 0, usum = 0;
    for(i = 0; i < 256; i++)  
    {  
        pixelPro[i] = (float)pixelCount[i] / totalCount;  
		utmp += i * pixelPro[i];  
		usum += i * i * pixelPro[i];
    }  
  
    //经典ostu算法,得到前景和背景的分割  
    //遍历灰度级[0,255],计算出方差最大的灰度值,为最佳阈值  
    float w0, w1, u0tmp, u0sum, u1sum, u1tmp, u0, u1, u,deltaTmp, deltaMax = 0;  
	w0 = u0tmp = u0sum = 0;
    for(i = 0; i < 256; i++)  
    {  
		w0 += pixelPro[i];
		u0tmp += i * pixelPro[i];  
		u0sum += i * i * pixelPro[i];

		w1 = 1 - w0;
		u1tmp = utmp - u0tmp;
		u1sum = usum - u0sum;

        u0 = u0tmp / w0;        //第一类的平均灰度  
        u1 = u1tmp / w1;        //第二类的平均灰度  
        u = u0tmp + u1tmp;      //整幅图像的平均灰度  
        //计算类间方差  
        deltaTmp = w0 * (u0 - u)*(u0 - u) + w1 * (u1 - u)*(u1 - u);
		float cov0 = (u0sum / w0 - u0 * u0);
		float cov1 = (u1sum / w1 - u1 * u1);
		float cov = 0.5 * (u0sum / w0 - u0 * u0) + 0.5 * (u1sum / w1 - u1 * u1);
        //找出最大类间方差以及对应的阈值  
        if(deltaTmp > deltaMax)  
        {     
            deltaMax = deltaTmp;  
            threshold = i;  
        }  
    }
    //返回最佳阈值;  
    return (float)threshold / 255;  
}  



void SaliencyEnhance(IplImage *pImgSrc, Mat_<float> &salImg)
{
	int W = pImgSrc->width, H = pImgSrc->height;
	int i, j;
	float th;
	th = OtsuThre(pImgSrc->width, pImgSrc->height, salImg);
	float thCenter = (1 - th) / th;
	for( j=0; j< H; j++ )
	{
		float *pSal = (float *)salImg.ptr(j);
		for( i=0; i< W; i++ )
		{
			*pSal = 1 / (1 + thCenter * exp(-wEnhance * (*pSal - th)));
			pSal ++;
		}
	}
}




int SaliencyWater(IplImage *pImgSrc, Mat_<float> &salImg)
{
	int frameWidth;
	FrameDetect(pImgSrc, frameWidth);
	cvSetImageROI(pImgSrc, cvRect(frameWidth, frameWidth, pImgSrc->width - 2 * frameWidth, pImgSrc->height - 2 * frameWidth));
	cvSmooth(pImgSrc, pImgSrc, CV_BLUR, 3);
	cvResetImageROI(pImgSrc);

	RGB2LAB(pImgSrc);

	SaliencyWaterRaw(pImgSrc, salImg, frameWidth);
	SaliencyCenterPrior(pImgSrc, salImg);

	SaliencySmooth(pImgSrc, salImg);
	normalize(salImg, salImg, 0, 1, NORM_MINMAX);
	SaliencyEnhance(pImgSrc, salImg);

	morphSmooth(pImgSrc->width, pImgSrc->height, salImg);
	normalize(salImg, salImg, 0, 1, NORM_MINMAX);

	return 0;
}
