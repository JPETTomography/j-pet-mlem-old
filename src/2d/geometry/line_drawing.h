#ifndef LINE_DRAWING
#define LINE_DRAWING

#if 0
void line_DDA_subpixel1(double x0,double y0,double x1,double y1,int col)    // DDA subpixel -> thick
    {

    // prepare data n-pixels,x1,y1 is line dx,dy step per pixel
    x1-=x0; i=ceil(fabs(x1));
    y1-=y0; n=ceil(fabs(y1));
    if (n<i) n=i;
    if (!n) n=1;
    auto dx = x1/double(n);
    auto dy = y1/double(n);
    n++;
    // rasterize DDA line
    int i,n,x,y,xx,yy;
    for (xx=x0,yy=y0,i=0;i<=n;i++,x0+=dx,y0+=dy) {
        // direct pixel
        pnt(x,y,col);
        // subpixels on change in both axises
        x=x0; y=y0;
        if ((i<n)&&(x!=xx)&&(y!=yy)) { pnt(xx,y,col); pnt(x,yy,col); }
        xx=x; yy=y;
        }
    }

#endif

#endif  // LINE_DRAWING
