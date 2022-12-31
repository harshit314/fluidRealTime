#include<iostream>
#include<SFML/Graphics.hpp>
#include<vector>
#include<math.h>
#include<array>
#include<string.h> //use memset function

using namespace std;
using namespace sf;

//window size:
Vector2i wSize(1200, 1200);
//number of cells:
int N = 240;
//cell width/height:
float h = (float)wSize.x/(float)N;

inline int IX(int i, int j) {return (i+(N+2)*j);} // i is horizontal, j is vertical indexing; N+2 cells in each direction

inline int cellIndxMap(int i, int j, int N)   //j for vertical, i for horizontal
{
    // return the 1D equivalent index of cell, where each cell has 4 points
    return 4*(i + N*j);
}

inline float lerp(float a, float b, float k) { return a+k*(b-a); }

void setBnd(int b, float *x)
{
    for (int i = 1; i <= N; i++)
    {
        x[IX(0  , i)] = b==1 ? -x[IX(1, i)] : x[IX(1, i)];
        x[IX(N+1, i)] = b==1 ? -x[IX(N, i)] : x[IX(N, i)];
        x[IX(i,   0)] = b==2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N+1)] = b==2 ? -x[IX(i, N)] : x[IX(i, N)];
    }
    //set edges:
    x[IX(0 ,   0 )] = 0.5*(x[IX(1,0 )] + x[IX(0 ,1)]);
    x[IX(0 , N+1 )] = 0.5*(x[IX(1,N+1)] + x[IX(0 ,N )]);
    x[IX(N+1,  0 )] = 0.5*(x[IX(N,0 )] + x[IX(N+1,1)]);
    x[IX(N+1, N+1)] = 0.5*(x[IX(N,N+1)] + x[IX(N+1,N )]);
    
}

void diffuseField(int b, float *x, float *x0, float gamma)   //pass 2d array by reference
{
    int iterMax = 1;
    for (int iter = 0; iter < iterMax; iter++)
    {
       for (int i = 1; i <= N; i++)
       {
            for (int j = 1; j <= N; j++)
            {
               // if(iter==0) x[IX(i,j)] = x0[IX(i,j)];
               // x[IX(i,j)] = (x0[IX(i,j)] + (gamma/4.f)*(x[IX(i-1,j)] + x[IX(i+1,j)] + x[IX(i,j-1)] + x[IX(i,j+1)]) )/(1.f + gamma);
                x[IX(i,j)] = (x0[IX(i,j)] + (gamma/4.f)*(x0[IX(i-1,j)] + x0[IX(i+1,j)] + x0[IX(i,j-1)] + x0[IX(i,j+1)]) )/(1.f + gamma);
            }
       }
       if(b!=-1) setBnd(b, x);
    }
}

void advectField(int b, float *d, float *d0, float *u, float *v, float dt)
{   
    float xPos, yPos, value1, value2;
    
    for (int i = 1; i <= N; i++)
    {
        for (int j = 1; j <= N; j++)
        {
            xPos = i*h - u[IX(i,j)]*dt;
            yPos = j*h - v[IX(i,j)]*dt;
            int ix = floor(xPos/h), iy = floor(yPos/h);  // integer part
            float jx = (xPos/h) - ix, jy = (yPos/h) - iy;    //fractional part
            if(ix<0 || ix>N)  { cout<<"Problem in advection!"<<endl; ix=0; }   
            if(iy<0 || iy>N)    iy=0;   //use zero values at boundaries 
            
            value1 = lerp(d0[IX(ix, iy)], d0[IX(ix + 1, iy)], jx);  // second index for horizontal and first for vertical, be careful...
            value2 = lerp(d0[IX(ix, iy + 1)], d0[IX(ix+1, iy+1)], jx);
            d[IX(i,j)] = lerp(value1, value2, jy);
        }

    }
    if(b!= -1)    setBnd(b, d);
}

void project(float * u, float * v, float * p, float * div)
{
    int iterMax = 15;   // should compare with prev iteration step to terminate iteration loop
    
    //get div field and p field:
    for (int iter = 0; iter < iterMax; iter++)
    {
       for (int i = 1; i <= N; i++)
       {
            for (int j = 1; j <= N; j++)
            {
                if(iter==0)    div[IX(i,j)] = -0.5*h*(u[IX(i+1, j)] - u[IX(i-1, j)] + v[IX(i, j+1)] - v[IX(i, j-1)] ); 
                p[IX(i,j)] = 0.25*(div[IX(i,j)] + (p[IX(i+1, j)] + p[IX(i-1, j)] + p[IX(i, j+1)] + p[IX(i, j-1)] ) );
            }
       }
        setBnd(0, p);
    }

    for (int i = 1; i <= N; i++)
    {
        for (int j = 1; j <= N; j++)
        {
            u[IX(i,j)] -= (p[IX(i+1, j)]-p[IX(i-1,j)])/(2.f*h);
            v[IX(i,j)] -= (p[IX(i, j+1)]-p[IX(i,j-1)])/(2.f*h);
        }
    }
    setBnd(1,u);  setBnd(2,v);
}

void densStep(float * x, float * x0, float * u, float * v, float diff, float dt )
{
    diffuseField(-1, x, x0, diff); 
    advectField(-1, x0, x, u, v, dt);    // advection in x occurs in next iteration
}

void velStep(float * u, float * v, float * u0, float * v0, float * p, float visc, float dt)
{
    advectField(1, u, u0, u0, v0, dt); advectField(2, v, v0, u0, v0, dt);
    project(u, v, p, v0);
}

int main()
{

    //define mesh:
   // int N = 4; // square mesh, size excludes the boundary cells; N is no. of cells in each direction
    int size = (N+2)*(N+2);
    //Velocity scale and time step:
    float Uscale = 200.f;
    float dt = 0.05, h = (float)wSize.x/(float)N;

    //define fields:
    float density[size] = {}; // 1d array of zeroes.
    float densityPrev[size] = {}; // 1d array of zeroes.
    
    float u[size] = {}; // 1d array of zeroes.
    float v[size] = {}; // 1d array of zeroes.
    float uPrev[size] = {}; // 1d array of zeroes.
    float vPrev[size] = {}; // 1d array of zeroes.
    float pField[size] = {};
    
    dt = 2.f*h/Uscale;
    dt = dt>0.08?dt:0.08;
    cout<<"dt: "<<dt<<endl;
    
    RenderWindow w(VideoMode(wSize.x, wSize.y), "improved2Dsolver");
    w.setFramerateLimit(60);

    //**********Create Mesh***********
    //define vertex array (quads) for each cell:
    vector<Vertex> cells; //4 vertices needed to define quad, need to be in a particular order (clockwise or anti-clockwise)  
    for (int j = 0; j < N; j++)  
    {
        for (int i = 0; i < N; i++)
        {
            Vertex temp;
            temp.position = Vector2f(i*h, j*h);
            float tmpVal = 0.25*(density[IX(i,j)]+density[IX(i+1,j)]+density[IX(i,j+1)]+density[IX(i+1,j+1)]);
            temp.color = Color(255*tmpVal);
            cells.push_back(temp);
            
            temp.position = Vector2f((i+1)*h, j*h);
            tmpVal = 0.25*(density[IX(i+1,j)]+density[IX(i+2,j)]+density[IX(i+1,j+1)]+density[IX(i+2,j+1)]);
            temp.color = Color(255*tmpVal);
            cells.push_back(temp);
            
            temp.position = Vector2f((i+1)*h, (j+1)*h);
            tmpVal = 0.25*(density[IX(i+1,j+1)]+density[IX(i+2,j+1)]+density[IX(i+1,j+2)]+density[IX(i+2,j+2)]);
            temp.color = Color(255*tmpVal);
            cells.push_back(temp);
            
            temp.position = Vector2f(i*h, (j+1)*h);
            tmpVal = 0.25*(density[IX(i,j+1)]+density[IX(i+1,j+1)]+density[IX(i,j+2)]+density[IX(i+1,j+2)]);
            temp.color = Color(255*tmpVal);
            cells.push_back(temp);
        }   
    }
    
    cout<<"Total no. of cells by meshSize squared: "<<cells.size()/(4*N*N)<<endl;   //should be 1
    cout<<"cellWidth: "<<h<<endl;
    //Use mouse click to add density
    Vector2i localMousePos(0, 0), localMousePrevPos(0, 0);
    Color currColor = Color::White;
    
    bool begin = false, mouseDown = false, draw = false;
    float dissipation = 0.05;

    while (w.isOpen())
    {
        Event e;
        while (w.pollEvent(e))
        {
            if(e.type==Event::Closed)   w.close();
            if(e.type==Event::KeyPressed)   
            {
                if(e.key.code==Keyboard::C) 
                {
                    for (int i = 0; i < size; i++)
                    {
                        density[i] = 0.f;
                        densityPrev[i] = 0.f;
                        u[i] = 0.f; v[i] = 0.f;
                        uPrev[i] = 0.f; vPrev[i] = 0.f;
                    }   
                }
                if(e.key.code==Keyboard::B) begin=begin?false:true;
                //draw boundaries:
                if(e.key.code==Keyboard::D)
                {
                    begin = false;
                    draw = draw?false:true;
                }
            }
            if(e.type==Event::MouseButtonReleased)
            {
                if(e.mouseButton.button==Mouse::Left)    mouseDown = false;
            }
            if(Mouse::isButtonPressed(Mouse::Left) && !mouseDown) 
            {   
                mouseDown = true;   
                localMousePrevPos = Mouse::getPosition(w); 
               // currColor = Color(255*(float)rand()/(float)RAND_MAX, 255*(float)rand()/(float)RAND_MAX, 255*(float)rand()/(float)RAND_MAX);
            }
        }
        
        //draw boundaries:
        if(draw)
        {
                


        }  

        if(begin)
        {
            //get fields from UI:
            localMousePos = Mouse::getPosition(w);                
            int sourceSize = 4; //radius of source 
            if(mouseDown)
            {
                int localI = floor(localMousePos.x/h), localJ = floor(localMousePos.y/h);

                int localIprev = floor(localMousePrevPos.x/h), localJprev = floor(localMousePrevPos.y/h);
                float delX = localMousePos.x-localMousePrevPos.x, delY = localMousePos.y-localMousePrevPos.y;
                float magn = sqrt(delX*delX + delY*delY);

                for (int jSource = localJ-sourceSize; jSource <= localJ+sourceSize; jSource++)
                {
                    if(jSource<0 || jSource>=N)   continue;
                    for (int iSource = localI-sourceSize; iSource <= localI+sourceSize; iSource++)
                    {
                        if(iSource<0 || iSource>=N)   continue;
                        
                        uPrev[IX(iSource+1, jSource+1)] = magn>0?Uscale*delX/magn:0.0;
                        vPrev[IX(iSource+1, jSource+1)] = magn>0?Uscale*delY/magn:-Uscale;
                        
                    //  vPrev[IX(iSource+1, jSource+1)] = -Uscale;

                        float localD = abs(jSource - localJ)*h + abs(iSource - localI)*h, dMax = sourceSize==0?1.f:sourceSize*(2.f*h);
                    
                        densityPrev[IX(iSource+1, jSource+1)] = densityPrev[IX(iSource+1, jSource+1)]<1.f-localD/(2.f*dMax)?1.f-localD/(2.f*dMax):densityPrev[IX(iSource+1, jSource+1)]; // added ones for the boundary offset
                    }    
                }
            } 

            velStep(u, v, uPrev, vPrev, pField, 0.0, dt);
        
            densStep(density, densityPrev, u, v, 0.01, dt);

            for (int i = 0; i < size; i++)
            {
                uPrev[i] = u[i];
                vPrev[i] = v[i];
                densityPrev[i] -= dissipation*densityPrev[i]*dt;
                density[i]=densityPrev[i];
            }  

            //*****************
            //update density for cells:
            for (int j = 0; j < N; j++)
            {
                for (int i = 0; i < N; i++)
                {
                    int cellCoord = cellIndxMap(i, j, N);
                    
                    float tempVal = 0.25*(densityPrev[IX(i, j)]+densityPrev[IX(i+1, j)]+densityPrev[IX(i, j+1)]+densityPrev[IX(i+1, j+1)]);
                    
                    Color tmpColor = currColor*Color(255.f*tempVal, 255.f*tempVal, 255.f*tempVal);
                    cells[cellCoord].color = tmpColor;
                    

                    tempVal = 0.25*(densityPrev[IX(i+1, j)]+densityPrev[IX(i+2, j)]+densityPrev[IX(i+1, j+1)]+densityPrev[IX(i+2, j+1)]);
                    tmpColor = currColor*Color(255.f*tempVal, 255.f*tempVal, 255.f*tempVal);
                    cells[cellCoord+1].color = tmpColor;
                    
                    
                    tempVal = 0.25*(densityPrev[IX(i+1, j+1)]+densityPrev[IX(i+2, j+1)]+densityPrev[IX(i+1, j+2)]+densityPrev[IX(i+2, j+2)]);
                    tmpColor = currColor*Color(255.f*tempVal, 255.f*tempVal, 255.f*tempVal);
                    cells[cellCoord+2].color = tmpColor;

                    tempVal = 0.25*(densityPrev[IX(i, j+1)]+densityPrev[IX(i+1, j+1)]+densityPrev[IX(i, j+2)]+densityPrev[IX(i+1, j+2)]);
                    tmpColor = currColor*Color(255.f*tempVal, 255.f*tempVal, 255.f*tempVal);
                    cells[cellCoord+3].color = tmpColor;
                                        
                }
            }
            
        }

        w.clear();
        w.draw(cells.data(), cells.size(), Quads);
        w.display();    
    
    
    }

    return 0;
}