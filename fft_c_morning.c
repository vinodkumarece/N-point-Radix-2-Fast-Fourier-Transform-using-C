

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_FFT_LEN                   1024

#define ISQRT_2                       (0.70710678118654752440084436210485)

double twiddle[MAX_FFT_LEN];
double in[MAX_FFT_LEN * 2];

//
// getLevels2
//     Get Number of levels for Radix-2
//

int getLevels2(int N)
{
    int k = 1;

    while((1 << k) != N)
        k++;

    return k;
}

//
// breverse
//     Reversal of given input in bits
//

int breverse(int value, int bitcnt)
{
    int i, res = 0;

    for(i = 0; i < bitcnt; i++)
        {
        if((value >> (bitcnt - 1 - i)) & 1)
            res |= (1 << i);
        }

    return res;
}

//
// brev_sort
//     Bit Reversal Sorting
//

void brev_sort(double *in, int n)
{
    int i, j, bits;
    double tmpRe, tmpIm;

    j = (n >> 1);
    bits = getLevels2(n);
    for(i = 0; i < n; i++)
        {
        if((j = breverse(i, bits)) != i && i < j)
            {
            tmpRe = in[2*i];
            tmpIm = in[2*i+1];

            in[2*i] = in[2*j];
            in[2*i+1] = in[2*j+1];

            in[2*j] = tmpRe;
            in[2*j+1] = tmpIm;
            }
        }
}

void genTwiddle2(double *twiddle, int N)
{
    int k, M;

    M = (N >> 1);
    for(k = 0; k < M; k++)
        {
        twiddle[2*k] = cos((3.14159265358979323846 * k)/M);
        twiddle[2*k+1] = -sin((3.14159265358979323846 * k)/M);
        }
}

//
// cfft_dit2_in
//     Radix-2 Decimation in Time Fast Fourier Transform
//

void cfft_dit2_in(double *in, double *twiddle, int N)
{
    int setIdx, bfyIdx, setsCnt, bfyCnt, twIncr;
    double *top, *bottom, *tw, tre, tim, bre, bim, tmp;

    brev_sort(in, N);
    setsCnt = (N >> 1);
    bfyCnt = 1;
    twIncr = N;

    top = in;
    bottom = in + bfyCnt * 2;

    for(setIdx = 0; setIdx < setsCnt; setIdx++)
        {
        tre = top[0];             tim = top[1];
        bre = top[2];             bim = top[3];
        top[0] = tre + bre;       top[1] = tim + bim;
        top[2] = tre - bre;       top[3] = tim - bim;

        top += 4;
        }

    twIncr >>= 1;
    setsCnt >>= 1;
    bfyCnt <<= 1;
    //printf("\nStage 0 passed\n");


    //Stage 1
    top = in;
    bottom = in + bfyCnt * 2;
    double tre1,tre2,tre3,tre4,tim1,tim2,tim3,tim4;
    for(setIdx = 0; setIdx < setsCnt; setIdx++)
    {
        tre1=top[0];                            tim1=top[1];
        tre2=top[2];                            tim2=top[3];
        tre3=top[4];                            tim3=top[5];
        tre4=top[7];                            tim4=top[6];

        top[0]=tre1+tre3;                       top[1]=tim1+tim3;
        top[2]=tre2+tre4;                       top[3]=tim2-tim4;
        top[4]=tre1-tre3;                       top[5]=tim1-tim3;
        top[6]=tre2-tre4;                       top[7]=tim2+tim4;
    
        top += 8;
    }
 
    twIncr >>= 1;
    setsCnt >>= 1;
    bfyCnt <<= 1;

    //printf("Stage 1 Passed\n");

    /*
    One more Approach
    top = in;
    bottom = in + bfyCnt * 2;
    for(setIdx = 0; setIdx < setsCnt; setIdx++)
    {
        tre1=top[0];                            tim1=top[1];
        tre2=top[2];                            tim2=top[3];
        tre3=bottom[0];                            tim3=bottom[1];
        tre4=bottom[2];                            tim4=bottom[3];

        top[0]=tre1+tre3;                       top[1]=tim1+tim3;
        top[2]=tre2+tre4;                       top[3]=tim2-tim4;
        top[4]=tre1-tre3;                       top[5]=tim1-tim3;
        top[6]=tre2-tre4;                       top[7]=tim2+tim4;
    
        top += 8;
        bottom+=8;
    }
    
    
    
    
    */





//Stage 2
top = in;
bottom = in + bfyCnt * 2;;

for (setIdx = 0; setIdx < setsCnt; setIdx++)
{
    double tre1,tre2,tre3,tre4,tre5,tre6,tre7,tre8,tim1,tim2,tim3,tim4,tim5,tim6,tim7,tim8;


    tre1=top[0];                                                tim1=top[1];
    tre2=top[2];                                                tim2=top[3];
    tre3=top[4];                                                tim3=top[5];
    tre4=top[6];                                                tim4=top[7];
    tre5=top[8];                                                tim5=top[9];
    tre6=ISQRT_2*(top[10]+top[11]);                             tim6=ISQRT_2*(top[11]-top[10]);
    tre7=top[13];                                               tim7=top[12];
    tre8=ISQRT_2*(top[15]-top[14]);                             tim8=ISQRT_2*(top[14]+top[15]);




    top[0]=tre1+tre5;                                           top[1]=tim1+tim5;
    top[2]=tre2+tre6;                                           top[3]=tim2+tim6;
    top[4]=tre3+tre7;                                           top[5]=tim3-tim7;
    top[6]=tre4+tre8;                                           top[7]=tim4-tim8;
    top[8]=tre1-tre5;                                           top[9]=tim1-tim5;
    top[10]=tre2-tre6;                                          top[11]=tim2-tim6;
    top[12]=tre3-tre7;                                          top[13]=tim3+tim7;
    top[14]=tre4-tre8;                                          top[15]=tim4+tim8;                                                

    top+=16;



   /*  Another method of doing 
    tre = top[0];                                                       tim = top[1];
    bre = bottom[0];                                                    bim = bottom[1];
    top[0] = tre + bre;                                                 top[1] = tim + bim;
    bottom[0] = tre - bre;                                              bottom[1] = tim - bim;
    
    
    tre = top[2];                                                       tim = top[3];
    bre = ISQRT_2*(bottom[2]+bottom[3]);                                bim = ISQRT_2*(bottom[3]-bottom[2]);
    top[2] = tre + bre;                                                 top[3] = tim + bim;
    bottom[2] = tre - bre;                                              bottom[3] = tim - bim;
    
    
    tre = top[4];                                                       tim = top[5];
    bre = bottom[4];                                                    bim = bottom[5];
    top[4]=tre+bim;                                                     top[5]= tim-bre;
    bottom[4] = tre-bim;                                                bottom[5] = tim+bre;
    

    tre = top[6];                                                       tim = top[7];
    bre = ISQRT_2*(bottom[7]-bottom[6]);                                bim = -1*ISQRT_2*(bottom[7]+bottom[6]);
    top[6] = tre + bre;                                                 top[7] = tim + bim;
    bottom[6] = tre - bre;                                              bottom[7] = tim - bim;
    
    top += 16;
    bottom+=16; */
}


    twIncr >>= 1;
    setsCnt >>= 1;
    bfyCnt <<= 1;
    //printf("Stage 2 Completed");


    //------------------------------------------//



    while(setsCnt)
        {
        top = in;
        bottom = in + bfyCnt * 2;

        for(setIdx = 0; setIdx < setsCnt; setIdx++)
            {
            tw = twiddle;
            for(bfyIdx = 0; bfyIdx < bfyCnt; bfyIdx++)
                {
                bre = bottom[bfyIdx * 2];             bim = bottom[bfyIdx * 2 + 1];
                tre = tw[0];                          tim = tw[1];
                tmp = bre * tre - bim * tim;          bim = bre * tim + bim * tre;
                bre = tmp;

                tre = top[bfyIdx * 2];                tim = top[bfyIdx * 2 + 1];
                top[bfyIdx * 2] = tre + bre;          top[bfyIdx * 2 + 1] = tim + bim;
                bottom[bfyIdx * 2] = tre - bre;       bottom[bfyIdx * 2 + 1] = tim - bim;

                tw += twIncr;
                }

            top += bfyCnt * 4;
            bottom += bfyCnt * 4;
            }

        twIncr >>= 1;
        setsCnt >>= 1;
        bfyCnt <<= 1;
        }
}

//
// main
//     Entry point
//

int main()
{
#if 1
    int i, n;
    double maxError, sqnr;

    for(n = 16; n < MAX_FFT_LEN; n = n * 2)
        {
        // Get twiddle factors
        genTwiddle2(twiddle, n);

        // Generate input
        for(i = 0; i < n; i++)
            {
            in[2*i] = cos((2 * M_PI * i)/n);
            in[2*i+1] = 0;
            }

        // Compute FFT
        cfft_dit2_in(in, twiddle, n);

        // Verify output
        maxError = 0;
        for(i = 0; i <= n/2; i++)
            {
            if(i == 1)
                {
                if(fabs(in[2*i] - n/2) > maxError)
                    maxError = fabs(in[2*i] - n/2);
                }
            else
                {
                if(fabs(in[2*i]) > maxError)
                    maxError = fabs(in[2*i]);
                }
            if(fabs(in[2*i+1]) > maxError)
                maxError = fabs(in[2*i+1]);
            }
        sqnr = 20 * log10f(n/(2 * maxError));
        if(sqnr < 120)
            {
            printf("SQNR below the threshold\n");
            return -1;
            }
        }

    printf("All tests passed\n");
#else
    int i;

    // Generate input
    for(i = 0; i < 8; i++)
        {
        in[2*i] = i;
        in[2*i + 1] = 0;
        }

    // Get twiddle factors
    genTwiddle2(twiddle, 8);

    // Compute FFT
    cfft_dit2_in(in, twiddle, 8);

    for(i = 0; i < 8; i++)
        printf("%f + j%f\n", in[2*i], in[2*i+1]);

#endif
    return 0;
}
