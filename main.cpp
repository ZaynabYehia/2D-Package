#if defined(UNICODE) && !defined(_UNICODE)
#define _UNICODE
#elif defined(_UNICODE) && !defined(UNICODE)
#define UNICODE
#endif

#include <tchar.h>
#include <windows.h>
#include <bits/stdc++.h>
#include "menu.h"
using namespace std;

struct Vertex
{
    double x, y;
    Vertex(int xi = 0, int yi = 0)
    {
        x = xi;
        y = yi;
    }

};
bool InLeft(Vertex& v, int edge)
{
    return v.x >= edge;
}
bool InRight(Vertex& v, int edge)
{
    return v.x <= edge;
}
bool InTop(Vertex& v, int edge)
{
    return v.y >= edge;
}
bool InBottom(Vertex& v, int edge)
{
    return v.y <= edge;
}
Vertex VIntersect(Vertex& v1, Vertex& v2, int xedge)
{
    Vertex res;
    res.x = xedge;
    res.y = v1.y + (xedge - v1.x) * (v2.y - v1.y) / (v2.x - v1.x);
    return res;
}
Vertex HIntersect(Vertex& v1, Vertex& v2, int yedge)
{
    Vertex res;
    res.y = yedge;
    res.x = v1.x + (yedge - v1.y) * (v2.x - v1.x) / (v2.y - v1.y);
    return res;
}
typedef vector<Vertex>VertexList;
typedef bool (*IsInFunc)(Vertex& v, int edge);
typedef Vertex (*IntersectFunc)(Vertex& v1, Vertex& v2, int edge);
struct Vector2
{
    double x,y;
    Vector2(double a=0,double b=0)
    {
        x=a;
        y=b;
    }
};
class Vector4
{
    double v[4];
public:
    Vector4(double a=0,double b=0,double c=0,double d=0)
    {
        v[0]=a;
        v[1]=b;
        v[2]=c;
        v[3]=d;
    }
    Vector4(double a[])
    {
        memcpy(v,a,4*sizeof(double));
    }
    double& operator[](int i)
    {
        return v[i];
    }
};
void Draw4Points(HDC hdc,int xc,int yc, int a, int b,COLORREF color)
{
    SetPixel(hdc, xc+a, yc+b, color);
    SetPixel(hdc, xc-a, yc+b, color);
    SetPixel(hdc, xc-a, yc-b, color);
    SetPixel(hdc, xc+a, yc-b, color);

}
class Matrix4
{
    Vector4 M[4];
public:
    Matrix4(double A[])
    {
        memcpy(M,A,16*sizeof(double));
    }
    Vector4& operator[](int i)
    {
        return M[i];
    }
};
Vector4 operator*(Matrix4 M,Vector4& b) // right multiplication of M by b
{
    Vector4 res;
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            res[i]+=M[i][j]*b[j];

    return res;
}
double DotProduct(Vector4& a,Vector4& b) //multiplying a raw vector by a column vector
{
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3];
}
Vector4 GetHermiteCoeff(double x0,double s0,double x1,double s1)
{
    static double H[16]= {2,1,-2,1,-3,-2,3,-1,0,1,0,0,1,0,0,0};
    static Matrix4 basis(H);
    Vector4 v(x0,s0,x1,s1);
    return basis*v;
}

class bezier_curve
{
public:
    int counter = 0;
    Vector2 p0, p1, p2, p3;
    HDC hdc;
    void lbutton(LPARAM lParam, HWND hwnd)
    {

        if(counter == 0)
        {
            p0.x = LOWORD(lParam);
            p0.y = HIWORD(lParam);
        }
        else if(counter == 1)
        {
            p1.x = LOWORD(lParam);
            p1.y = HIWORD(lParam);

        }
        else if(counter == 2)
        {
            p2.x = LOWORD(lParam);
            p2.y = HIWORD(lParam);
        }
        else
        {
            p3.x = LOWORD(lParam);
            p3.y = HIWORD(lParam);
            hdc = GetDC(hwnd);
            DrawBezierCurve(hdc, p0, p1, p2, p3, 1000);
            ReleaseDC(hwnd, hdc);
            counter = -1;
        }
        counter++;

    }

    void DrawHermiteCurve (HDC hdc,Vector2& P0,Vector2& T0,Vector2& P1,Vector2& T1,int numpoints)
    {
        Vector4 xcoeff=GetHermiteCoeff(P0.x,T0.x,P1.x,T1.x);
        Vector4 ycoeff=GetHermiteCoeff(P0.y,T0.y,P1.y,T1.y);
        if(numpoints<2)return;
        double dt=1.0/(numpoints-1);
        for(double t=0; t<=1; t+=dt)
        {
            Vector4 vt;
            vt[3]=1;
            for(int i=2; i>=0; i--)vt[i]=vt[i+1]*t;
            int x=round(DotProduct(xcoeff,vt));
            int y=round(DotProduct(ycoeff,vt));
            if(t==0)MoveToEx(hdc,x,y,NULL);
            else LineTo(hdc,x,y);
        }
    }
    void DrawBezierCurve(HDC hdc,Vector2& P0,Vector2& P1,Vector2& P2,Vector2& P3,int numpoints)
    {
        Vector2 T0(3*(P1.x-P0.x),3*(P1.y-P0.y));
        Vector2 T1(3*(P3.x-P2.x),3*(P3.y-P2.y));
        DrawHermiteCurve(hdc,P0,T0,P3,T1,numpoints);
    }
};
class direct_ellipse
{
public:
    int a = 50, b = 40, xc, yc;
    HDC hdc;
    void lbutton(HWND hwnd, LPARAM lParam)
    {

        xc = LOWORD(lParam);
        yc = HIWORD(lParam);
        hdc = GetDC(hwnd);
        EllipseDirect( hdc, xc, yc,  a,  b,(41, 41, 41));
        ReleaseDC(hwnd, hdc);
    }
    void EllipseDirect(HDC hdc,int xc,int yc, int a, int b,COLORREF color)
    {
        int x=0,y=b;
        double b2 = b*b, a2 = a*a;
        Draw4Points(hdc,xc,yc,x,y,color);
        while(b*b*x<a*a*y)
        {
            x++;
            y=round(sqrt(b2*(double)(1.0-1.0*x*x/a2)));
            Draw4Points(hdc,xc,yc,x,y,color);
        }
        while(y>0)
        {
            y--;
            x=round(sqrt(a2*(double)(1.0-1.0*y*y/b2)));
            Draw4Points(hdc,xc,yc,x,y,color);
        }
    }



};
class cardinal_spline
{
public:
    vector<Vector2>p;
    Vector2 tmp;
    int sz;
    HDC hdc;
    lbutton(LPARAM lParam)
    {
        tmp.x = LOWORD(lParam);
        tmp.y = HIWORD(lParam);
        p.push_back(tmp);



    }
    void rbutton(HWND hwnd)
    {
        hdc = GetDC(hwnd);
        DrawCardinalSpline(hdc, p, p.size(), 0.5, 1000);
        p.clear();

    }
    void DrawHermiteCurve (HDC hdc,Vector2& P0,Vector2& T0,Vector2& P1,Vector2& T1,int numpoints)
    {
        Vector4 xcoeff=GetHermiteCoeff(P0.x,T0.x,P1.x,T1.x);
        Vector4 ycoeff=GetHermiteCoeff(P0.y,T0.y,P1.y,T1.y);
        if(numpoints<2)return;
        double dt=1.0/(numpoints-1);
        for(double t=0; t<=1; t+=dt)
        {
            Vector4 vt;

            vt[3]=1;
            for(int i=2; i>=0; i--)vt[i]=vt[i+1]*t;
            int x=round(DotProduct(xcoeff,vt));
            int y=round(DotProduct(ycoeff,vt));
            // cout<< x << " " << y << endl;
            if(t==0)MoveToEx(hdc,x,y,NULL);
            else LineTo(hdc,x,y);
        }
    }

    void DrawCardinalSpline(HDC hdc, vector<Vector2> P, int n, double c, int numpix)
    {
        double c1 = 1 - c;
        Vector2 T0(c1*(P[2].x - P[0].x), c1*(P[2].y - P[0].y));
        for(int i = 2; i < n - 1; i++)
        {
            Vector2 T1(c1*(P[i + 1].x - P[i - 1].x), c1*(P[i + 1].y - P[i - 1].y));
            // cout << P[i].y << " " <<P[i].x<<endl;
            DrawHermiteCurve(hdc,P[i-1],T0,P[i],T1,numpix);
            T0=T1;
        }
    }


};
class clipping
{
public:
    int counter = 0, x_left, y_top, x_right, y_down;
    VertexList p;
    HDC hdc;
    void rbutton()
    {
        PolygonClip(hdc,p,p.size(),x_left, y_top,x_right,y_down);

    }
    void lbutton(LPARAM lParam, HWND hwnd)
    {
        if(counter == 0)
        {
            x_left = LOWORD(lParam);
            y_top = HIWORD(lParam);
        }
        else if(counter == 1)
        {
            hdc = GetDC(hwnd);
            x_right = LOWORD(lParam);
            y_down = HIWORD(lParam);
            MoveToEx(hdc, x_left, y_top, NULL);
            LineTo(hdc, x_right, y_top);
            MoveToEx(hdc, x_right, y_top, NULL);
            LineTo(hdc, x_right, y_down);
            MoveToEx(hdc, x_right, y_down, NULL);
            LineTo(hdc, x_left, y_down);
            MoveToEx(hdc, x_left, y_down, NULL);
            LineTo(hdc, x_left, y_top);


        }
        else
        {
            Vertex v(LOWORD(lParam), HIWORD(lParam));
            p.push_back(v);

        }

        counter++;


    }
    VertexList ClipWithEdge(VertexList& p, int edge, IsInFunc In, IntersectFunc Intersect)
    {
        VertexList OutList;
        Vertex v1 = p.back();
        bool v1_In = In(v1, edge);
        for(int i = 0; i < p.size(); i++)
        {
            Vertex v2 = p[i];
            bool v2_In = In(v2, edge);
            if(v1_In && v2_In)
                OutList.push_back(v2);
            else if(v1_In && !v2_In)
                OutList.push_back(Intersect(v1, v2, edge));
            else if(!v1_In && v2_In)
            {
                OutList.push_back(Intersect(v1, v2, edge));
                OutList.push_back(v2);
            }
            v1 = v2;
            v1_In = v2_In;
        }

        return OutList;
    }

    void PolygonClip(HDC hdc,VertexList p,int n,int xleft,int ytop,int xright,int ybottom)
    {
        VertexList vlist;
        for(int i = 0; i < n; i++) vlist.push_back(p[i]);
        vlist=ClipWithEdge(vlist, xleft, InLeft, VIntersect);
        vlist=ClipWithEdge(vlist, xright, InRight, VIntersect);
        vlist=ClipWithEdge(vlist, ytop, InTop, HIntersect);
        vlist=ClipWithEdge(vlist, ybottom, InBottom, HIntersect);
        Vertex v1=vlist[vlist.size()-1];
        for(int i=0; i<(int)vlist.size(); i++)
        {
            Vertex v2=vlist[i];
            MoveToEx(hdc,v1.x,round(v1.y),NULL);
            LineTo(hdc,round(v2.x),round(v2.y));
            v1=v2;
        }

    }


};
class polar_ellipse
{
public:
    int a = 50, b = 40, xc, yc;
    HDC hdc;
    void lbutton(LPARAM lParam, HWND hwnd)
    {

        xc = LOWORD(lParam);
        yc = HIWORD(lParam);
        hdc = GetDC(hwnd);
        polarellipse( hdc, xc, yc,  a,  b,(41, 41, 41));
        ReleaseDC(hwnd, hdc);


    }
    void polarellipse(HDC hdc,int xc,int yc,int  a,int  b,COLORREF color)
    {
        int x = a, y = 0;
        Draw4Points(hdc,xc,yc,x,y,color);
        double f = sqrt(a*a - b*b), e = f / a, b2 = b*b;
        for(double theta = 0.001; theta <= 3.14/2; theta+=0.001)
        {
            double R = a*b / sqrt((b*b*cos(theta)*cos(theta)) + (a*a*sin(theta)*sin(theta)));
            x=round(R*cos(theta));
            y=round(R*sin(theta));
            Draw4Points(hdc,xc,yc,x,y,color);
        }
    }



};
class midpoint_ellipse
{
public:

    int a = 50, b = 40, xc, yc;
    HDC hdc;
    void lbutton(LPARAM lParam, HWND hwnd)
    {
        xc = LOWORD(lParam);
        yc = HIWORD(lParam);
        hdc = GetDC(hwnd);
        midptellipse( hdc, xc, yc,  a,  b,(41, 41, 41));
        ReleaseDC(hwnd, hdc);


    }
    void midptellipse(HDC hdc, int xc,int yc, int a,int  b,COLORREF color)
    {
        int ry = b, rx = a;
        float dx, dy, d1, d2, x, y;
        x = 0;
        y = ry;

        // Initial decision parameter of region 1
        d1 = (ry * ry) - (rx * rx * ry) +
             (0.25 * rx * rx);
        dx = 2 * ry * ry * x;
        dy = 2 * rx * rx * y;

        // For region 1
        while (dx < dy)
        {

            // Print points based on 4-way symmetry
            Draw4Points(hdc, xc,yc, x, y, color);

            // Checking and updating value of
            // decision parameter based on algorithm
            if (d1 < 0)
            {
                x++;
                dx = dx + (2 * ry * ry);
                d1 = d1 + dx + (ry * ry);
            }
            else
            {
                x++;
                y--;
                dx = dx + (2 * ry * ry);
                dy = dy - (2 * rx * rx);
                d1 = d1 + dx - dy + (ry * ry);
            }
        }

        // Decision parameter of region 2
        d2 = ((ry * ry) * ((x + 0.5) * (x + 0.5))) +
             ((rx * rx) * ((y - 1) * (y - 1))) -
             (rx * rx * ry * ry);

        // Plotting points of region 2
        while (y >= 0)
        {

            // Print points based on 4-way symmetry
            Draw4Points(hdc, xc,yc, x, y, color);

            // Checking and updating parameter
            // value based on algorithm
            if (d2 > 0)
            {
                y--;
                dy = dy - (2 * rx * rx);
                d2 = d2 + (rx * rx) - dy;
            }
            else
            {
                y--;
                x++;
                dx = dx + (2 * ry * ry);
                dy = dy - (2 * rx * rx);
                d2 = d2 + dx - dy + (rx * rx);
            }
        }
    }




};
class hermite_curve
{
public:
    int counter = 0;
    Vector2 p0, p1, t0, t1;
    HDC hdc;
    void lbutton(LPARAM lParam, HWND hwnd)
    {
        if(counter == 0)
        {
            p0.x = LOWORD(lParam);
            p0.y = HIWORD(lParam);
        }
        else if(counter == 1)
        {
            t0.x = LOWORD(lParam) - p0.x;
            t0.y = HIWORD(lParam) - p0.y;

        }
        else if(counter == 2)
        {
            p1.x = LOWORD(lParam);
            p1.y = HIWORD(lParam);
        }
        else
        {
            t1.x = LOWORD(lParam) - p1.x;
            t1.y = HIWORD(lParam) - p1.y;
            hdc = GetDC(hwnd);
            DrawHermiteCurve (hdc, p0, t0, p1, t1, 1000);
            ReleaseDC(hwnd, hdc);
            counter = -1;
        }
        counter++;


    }
    void DrawHermiteCurve (HDC hdc,Vector2& P0,Vector2& T0,Vector2& P1,Vector2& T1,int numpoints)
    {
        Vector4 xcoeff=GetHermiteCoeff(P0.x,T0.x,P1.x,T1.x);
        Vector4 ycoeff=GetHermiteCoeff(P0.y,T0.y,P1.y,T1.y);
        if(numpoints<2)return;
        double dt=1.0/(numpoints-1);
        for(double t=0; t<=1; t+=dt)
        {
            Vector4 vt;
            vt[3]=1;
            for(int i=2; i>=0; i--)vt[i]=vt[i+1]*t;
            int x=round(DotProduct(xcoeff,vt));
            int y=round(DotProduct(ycoeff,vt));
            if(t==0)MoveToEx(hdc,x,y,NULL);
            else LineTo(hdc,x,y);
        }
    }


};
class Point
{
public:
    double x;
    double y;

    Point(){}

    Point(double x, double y);
};
class parametric
{
public:
    int counter = 0;
    Point p1;
    Point p2;
    void lbutton(LPARAM lParam, HWND hwnd)
    {
        if(counter == 0)
        {
            p1.x = LOWORD(lParam);
            p1.y = HIWORD(lParam);
        }
        else
        {

            p2.x = LOWORD(lParam);
            p2.y = HIWORD(lParam);
            HDC hdc = GetDC(hwnd);
            Parametric(hdc, p1, p2);
            counter = -1;
        }
        ++counter;
    }
    void Parametric(HDC hdc, Point startPoint, Point endPoint)
    {
        int numberOfPoints = max(abs(startPoint.x - endPoint.x), abs(startPoint.y - endPoint.y));
        double dt = 1.0 / numberOfPoints;
        double x = startPoint.x, y = startPoint.y;
        double dx = dt * (endPoint.x - startPoint.x);
        double dy = dt * (endPoint.y - startPoint.y);
        for (double t = 1; t <= numberOfPoints; t++)
        {
            SetPixel(hdc, round(x), round(y), RGB(26, 29, 94));
            x += dx;
            y += dy;
        }
    }
};
class dda
{

public:
    int counter = 0;
    Point p1, p2;
    void lbutton(LPARAM lParam, HWND hwnd)
    {
        if(counter == 0)
        {
            p1.x = LOWORD(lParam);
            p1.y = HIWORD(lParam);
        }
        else
        {

            p2.x = LOWORD(lParam);
            p2.y = HIWORD(lParam);
            HDC hdc = GetDC(hwnd);
            DDA(hdc, p1, p2);
            counter = -1;
        }
        ++counter;
    }
    void DDA(HDC hdc, Point startPoint, Point endPoint)
    {
        double dx = endPoint.x - startPoint.x;
        double dy = endPoint.y - startPoint.y;

        if (abs(dy) < abs(dx)) // slope < 1, x is the independent variable
        {
            double slope = dy / dx;

            if (startPoint.x > endPoint.x)
            {
                swap(startPoint.x, endPoint.x);
                swap(startPoint.y, endPoint.y);
            }

            int x = startPoint.x;
            double y = startPoint.y;
            SetPixel(hdc, x, y, RGB(41, 41, 41));

            while (x < endPoint.x)
            {
                x++;
                y += slope;
                SetPixel(hdc, x, y,RGB(41, 41, 41));
            }
        }
        else   // slope > 1, y is the independent variable
        {
            double slope = dx / dy;

            if (startPoint.y > endPoint.y)
            {
                swap(startPoint.x, endPoint.x);
                swap(startPoint.y, endPoint.y);
            }

            double x = startPoint.x;
            int y = startPoint.y;
            SetPixel(hdc, x, y, RGB(41, 41, 41));

            while (y < endPoint.y)
            {
                y++;
                x += slope;
                SetPixel(hdc, x, y, RGB(41, 41, 41));
            }
        }
    }
};
class bresenham
{
    int counter = 0;
    Point p1, p2;
public:
    void lbutton(LPARAM lParam, HWND hwnd)
    {
        if(counter == 0)
        {
            p1.x = LOWORD(lParam);
            p1.y = HIWORD(lParam);
        }
        else
        {

            p2.x = LOWORD(lParam);
            p2.y = HIWORD(lParam);
            HDC hdc = GetDC(hwnd);
            Bresenham(hdc, p1, p2);
            counter = -1;
        }
        ++counter;
    }
    void Bresenham(HDC hdc, Point startPoint, Point endPoint)
    {
        int ixs = (int) startPoint.x; // converting all double variables to ints, as this algorithm works all integers
        int iys = (int) startPoint.y;
        int ixe = (int) endPoint.x;
        int iye = (int) endPoint.y;

        int deltaX = ixe - ixs;
        int deltaY = iye - iys;

        if (abs(deltaY) <= abs(deltaX))
        {
            if (ixs > ixe)
            {
                swap(ixs, ixe);
                swap(iys, iye);
            }

            deltaX = abs(deltaX);
            deltaY = abs(deltaY);
            int error = 2 * deltaY - deltaX;
            int d1 = 2 * deltaY;
            int d2 = 2 * (deltaY - deltaX);

            int x = ixs;
            int y = iys;

            int increment;
            if (iys < iye)
                increment = 1;
            else
                increment = -1;

            SetPixel(hdc, x, y, RGB(0, 0, 0));
            while (x < ixe)
            {
                if (error < 0)
                    error += d1;
                else
                {
                    error += d2;
                    y += increment;
                }
                x++;
                SetPixel(hdc, x, y, RGB(0, 0, 0));
            }
        }
        else
        {
            if (iys > iye)
            {
                swap(ixs, ixe);
                swap(iys, iye);
            }

            deltaX = abs(deltaX);
            deltaY = abs(deltaY);

            int error = 2 * deltaX - deltaY;
            int d1 = 2 * deltaX;
            int d2 = 2 * (deltaX - deltaY);

            int x = ixs;
            int y = iys;

            int increment;
            if (ixs < ixe)
                increment = 1;
            else
                increment = -1;

            SetPixel(hdc, x, y, RGB(0, 0, 0));
            while (y < iye)
            {
                if (error < 0)
                    error += d1;
                else
                {
                    error += d2;
                    x += increment;
                }
                y++;
                SetPixel(hdc, x, y, RGB(0, 0, 0));
            }
        }
    }
};
/*  Declare Windows procedure  */
LRESULT CALLBACK WindowProcedure (HWND, UINT, WPARAM, LPARAM);

/*  Make the class name into a global variable  */
TCHAR szClassName[ ] = _T("CodeBlocksWindowsApp");

int WINAPI WinMain (HINSTANCE hThisInstance,
                    HINSTANCE hPrevInstance,
                    LPSTR lpszArgument,
                    int nCmdShow)
{
    HWND hwnd;               /* This is the handle for our window */
    MSG messages;            /* Here messages to the application are saved */
    WNDCLASSEX wincl;        /* Data structure for the windowclass */

    /* The Window structure */
    wincl.hInstance = hThisInstance;
    wincl.lpszClassName = szClassName;
    wincl.lpfnWndProc = WindowProcedure;      /* This function is called by windows */
    wincl.style = CS_DBLCLKS;                 /* Catch double-clicks */
    wincl.cbSize = sizeof (WNDCLASSEX);

    /* Use default icon and mouse-pointer */
    wincl.hIcon = LoadIcon (NULL, IDI_APPLICATION);
    wincl.hIconSm = LoadIcon (NULL, IDI_APPLICATION);
    wincl.hCursor = LoadCursor (NULL, IDC_ARROW);
    wincl.lpszMenuName = MAKEINTRESOURCE(MY_MENU);                 /* No menu */
    wincl.cbClsExtra = 0;                      /* No extra bytes after the window class */
    wincl.cbWndExtra = 0;                      /* structure or the window instance */
    /* Use Windows's default colour as the background of the window */
    wincl.hbrBackground = (HBRUSH) COLOR_BACKGROUND;

    /* Register the window class, and if it fails quit the program */
    if (!RegisterClassEx (&wincl))
        return 0;

    /* The class is registered, let's create the program*/
    hwnd = CreateWindowEx (
               0,                   /* Extended possibilites for variation */
               szClassName,         /* Classname */
               _T("Code::Blocks Template Windows App"),       /* Title Text */
               WS_OVERLAPPEDWINDOW, /* default window */
               CW_USEDEFAULT,       /* Windows decides the position */
               CW_USEDEFAULT,       /* where the window ends up on the screen */
               544,                 /* The programs width */
               375,                 /* and height in pixels */
               HWND_DESKTOP,        /* The window is a child-window to desktop */
               NULL,                /* No menu */
               hThisInstance,       /* Program Instance handler */
               NULL                 /* No Window Creation data */
           );

    /* Make the window visible on the screen */
    ShowWindow (hwnd, nCmdShow);

    /* Run the message loop. It will run until GetMessage() returns 0 */
    while (GetMessage (&messages, NULL, 0, 0))
    {
        /* Translate virtual-key messages into character messages */
        TranslateMessage(&messages);
        /* Send message to WindowProcedure */
        DispatchMessage(&messages);
    }

    /* The program return-value is 0 - The value that PostQuitMessage() gave */
    return messages.wParam;
}


/*  This function is called by the Windows function DispatchMessage()  */
direct_ellipse directtellipse;
cardinal_spline cardinalspline;
bezier_curve beziercurve;
hermite_curve hermitecurve;
midpoint_ellipse midpointellipse;
polar_ellipse polarellipse;
clipping clip;
parametric par;
dda d;
bresenham b;
int flag = 10;
int var;
bool bezier,save, load, lineParametric, hermite, bazier, polygonClipping, lineDDA, lineMidPoint, EllipsePolar,  splines,EllipseDirect,EllipseMidPoint;
LRESULT CALLBACK WindowProcedure (HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam)
{

    switch (message)                  /* handle the messages */
    {
    case WM_LBUTTONDOWN:
    {

        switch (flag)
        {
        case 1:
        {
            directtellipse.lbutton(hwnd, lParam);
            break;
        }
        case 2:
        {
            cardinalspline.lbutton(lParam);
            break;
        }
        case 3:
        {
            beziercurve.lbutton(lParam, hwnd);
            break;
        }
        case 4:
        {


            hermitecurve.lbutton(lParam, hwnd);
            break;}
        case 5:
        {
            midpointellipse.lbutton(lParam, hwnd);
            break;
        }
        case 6:
        {
            polarellipse.lbutton(lParam, hwnd);
            break;
        }
        case 7:
        {
            clip.lbutton(lParam, hwnd);
            break;
        }
        case 8:
        {
           par.lbutton(lParam, hwnd);
        }
        case 9:
        {
            d.lbutton(lParam, hwnd);
        }
        case 10:
        {
           b.lbutton(lParam, hwnd);
        }

        }
        break;
    }
    case WM_LBUTTONUP:
    {
        break;
    }
    case WM_RBUTTONDOWN:
    {

        switch (flag)
        {

        case 2:
        {
            cardinalspline.rbutton(hwnd);
            break;
        }


        case 7:
        {
            clip.rbutton();
            break;
        }

        }
        break;
    }
    case WM_RBUTTONUP:
    {
        break;
    }
    case WM_COMMAND:
        {
            var = LOWORD(wParam);
            switch (var) {
                case Save_ID:
                    save = true;
                    break;
                case Load_ID:
                    load = true;
                    break;
                case Exit_ID:
                    PostQuitMessage(0);
                    break;
                case LineParametric_ID:
                    lineParametric = true;
                    flag = 8;
                    break;
                case LineDDA_ID:
                    lineDDA = true;
                    flag = 9;
                    break;
                case LineMidPoint_ID:
                    lineMidPoint = true;
                    flag = 10;
                    break;
                case EllipseDirect_ID:
                    EllipseDirect = true;
                    flag = 1;
                    break;
                case EllipsePolar_ID:
                    EllipsePolar = true;
                    flag = 6;
                    break;
                case EllipseMidPoint_ID:
                    EllipseMidPoint = true;
                    flag = 5;
                    break;
                case PolygonClipping_ID:
                    polygonClipping = true;
                    flag = 7;
                    break;
                case Bezier_ID:
                    bezier = true;
                    flag = 3;
                    break;
                case Hermite_ID:
                    hermite = true;
                    flag = 4;
                    break;
                case Splines_ID:
                    splines = true;
                      flag = 2;
                    break;
                default:
                    return DefWindowProc(hwnd, message, wParam, lParam);


        }
        break;
        }
    case WM_DESTROY:
        PostQuitMessage (0);       /* send a WM_QUIT to the message queue */
        break;
    default:                      /* for messages that we don't deal with */
        return DefWindowProc (hwnd, message, wParam, lParam);
    }

    return 0;
}
