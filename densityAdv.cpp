#include<iostream>
#include<SFML/Graphics.hpp>
#include<vector>
#include<math.h>
#include<string.h> //use memset function

using namespace std;
using namespace sf;

//window size:
Vector2i wSize(1200, 1200);

//define mesh:
const int meshSize = 400, meshSizeP2 = 402; // square mesh, size excludes the boundary cells; meshSize is no. of cells in each direction
float meshWidth, meshHeight;

int cellIndxMap(int i, int j)   //i for vertical, j for horizontal
{
    // return the 1D equivalent index of cell, where each cell has 4 points
    return 4*(j + meshSize*i);
}

//velocity function:
int frameRate = 60;
float omega = 6*M_PI/20.f;
float f(float x, float t)
{
    float epsilon = 0.25;
    return epsilon*sin(omega*t)*x*x + (1-2*epsilon*sin(omega*t))*x;
}
float fPrime(float x, float t)
{
    float epsilon = 0.25;
    return 2*epsilon*sin(omega*t)*x + (1-2*epsilon*sin(omega*t));
}
float dt = (1*2.f*M_PI/omega)/(float)frameRate; // half complete cycle in a second


float lerp(float a, float b, float k)
{
    return a+k*(b-a);
}

void diffuseField(float field[][meshSizeP2], float delta)   //pass 2d array by reference
{
    float fieldNext[meshSizeP2][meshSizeP2] = {};
    float del = delta;    //del = 4*kappa*dt/h^2, it controls the diffusion strength
    int iterMax = 8;   // should compare with prev iteration step to terminate iteration loop
    for (int iter = 0; iter <= iterMax; iter++)
    {
       for (int i = 0; i < meshSize; i++)
       {
           for (int j = 0; j < meshSize; j++)
           {
               if(iter==0) fieldNext[i+1][j+1] = field[i+1][j+1];  //start with prev field, boundary cells are already always zero
               else if(iter<iterMax && iter>0)    fieldNext[i+1][j+1] = (fieldNext[i+1][j+1] + del*(1.0/4.0)*(field[i+1+1][j+1]+field[i+1-1][j+1]+field[i+1][j+1+1]+field[i+1][j+1-1] ))/(1+del);
               else field[i+1][j+1] = fieldNext[i+1][j+1];
           }
       }
    }
}

void advectField(float field[][meshSizeP2], float uxField[][meshSizeP2], float uyField[][meshSizeP2])
{
    float xPos, yPos, value1, value2;
    float tempField[meshSizeP2][meshSizeP2] = {};
    for (int i = 0; i < meshSize; i++)
    {
        for (int j = 0; j < meshSize; j++)
        {
            xPos = j*meshWidth - uxField[i+1][j+1]*dt;
            yPos = i*meshHeight - uyField[i+1][j+1]*dt;
            int ix = floor(xPos/(float)meshWidth), iy = floor(yPos/(float)meshHeight);  // integer part
            float jx = (xPos/(float)meshWidth) - ix, jy = (yPos/(float)meshHeight) - iy;    //fractional part
           // ix = jx<=0.5?ix-1:ix;
           // iy = jy<=0.5?iy-1:iy;   //no need to change jx and jy
            if(ix+1<0 || ix+1>meshSize)  { cout<<"Problem in advection!"<<endl; ix=0; }   
            if(iy+1<0 || iy+1>meshSize)    iy=0;   //use zero values at boundaries 
            
            value1 = lerp(field[iy+1][ix+1], field[iy+1][ix+1+1], jx);  // second index for horizontal and first for vertical, be careful...
            value2 = lerp(field[iy+1+1][ix+1], field[iy+1+1][ix+1+1], jx);
            tempField[i+1][j+1] = lerp(value1, value2, jy);
        }
    }
    for (int i = 0; i < meshSize; i++)
    {
        for (int j = 0; j < meshSize; j++)
        {
            field[i+1][j+1] = tempField[i+1][j+1];
        }
    }
    
}

int main()
{
    //define fields:
    float density[meshSizeP2][meshSizeP2] = {}; // 2d array of zeroes.
    float ux[meshSizeP2][meshSizeP2] = {}; // 2d array of zeroes.
    float uy[meshSizeP2][meshSizeP2] = {}; // 2d array of zeroes.


    meshWidth = (float)wSize.x/(float)meshSize;
    meshHeight = (float)wSize.y/(float)meshSize;

    RenderWindow w(VideoMode(wSize.x, wSize.y), "densityAdvection");
    w.setFramerateLimit(frameRate);

    //**********Create Mesh***********
    //define vertex array (quads) for each cell:
    vector<Vertex> cells; //4 vertices needed to define quad, need to be in a particular order (clockwise or anti-clockwise)  
    float Uscale = 50.0;   
    for (int i = 0; i < meshSize; i++)
    {
        for (int j = 0; j < meshSize; j++)
        {
            //set velocity field:
            float x = 2*(j+1)*meshWidth/wSize.x, y = (i+1)*meshHeight/wSize.y;
            ux[i+1][j+1] = Uscale*M_PI*(sin(M_PI*x)*cos(M_PI*y));
            uy[i+1][j+1] = -Uscale*M_PI*(cos(M_PI*x)*sin(M_PI*y));
            
            Vertex temp;
            temp.color = Color(255*density[i+1][j+1]);
            temp.position = Vector2f(j*meshWidth, i*meshHeight);
            cells.push_back(temp);
            
            temp.position = Vector2f((j+1)*meshWidth, i*meshHeight);
            cells.push_back(temp);
            
            temp.position = Vector2f((j+1)*meshWidth, (i+1)*meshHeight);
            cells.push_back(temp);
            
            temp.position = Vector2f(j*meshWidth, (i+1)*meshHeight);
            cells.push_back(temp);
        }   
    }
    
    cout<<"Total no. of cells by meshSize squared: "<<cells.size()/(4*meshSize*meshSize)<<endl;   //should be 1
   
    //Use mouse click to add density
    Vector2i localMousePos(0.f, 0.f);
    float t=0.f;
   // Vector3f currColor(255, 255, 255);
    while (w.isOpen())
    {
        Event e;
        while (w.pollEvent(e))
        {
            if(e.type==Event::Closed)   w.close();
            if(e.type==Event::KeyPressed)   
            {
                if(e.key.code==Keyboard::C) memset(density, 0.f, sizeof(density));
            }
          //  if (Mouse::isButtonPressed(Mouse::Left))    currColor = Vector3f(255*(float)rand()/(float)RAND_MAX, 255*(float)rand()/(float)RAND_MAX, 255*(float)rand()/(float)RAND_MAX);
        }

        //add source:
        localMousePos = Mouse::getPosition(w);
        int sourceSize = 5; //radius of source 
        if(Mouse::isButtonPressed(Mouse::Left))
        {
            int localJ = floor(localMousePos.x/meshWidth), localI = floor(localMousePos.y/meshHeight);
            for (int iSource = localI-sourceSize; iSource <= localI+sourceSize; iSource++)
            {
                if(iSource<0 || iSource>=meshSize)   continue;
                for (int jSource = localJ-sourceSize; jSource <= localJ+sourceSize; jSource++)
                {
                    if(jSource<0 || jSource>=meshSize)   continue;
                    float localD = abs(jSource - localJ)*meshWidth + abs(iSource - localI)*meshHeight, dMax = sourceSize*(meshWidth+meshHeight);
                    
                    density[iSource + 1][jSource + 1] = density[iSource + 1][jSource + 1]<1.f-localD/(2.f*dMax)?1.f-localD/(2.f*dMax):density[iSource + 1][jSource + 1]; // added ones for the boundary offset
                    int cellCoord = cellIndxMap(iSource, jSource);
                    for(int cIndx = 0; cIndx < 4; cIndx++)  cells[cellCoord + cIndx].color = Color::White;
                }   
            }
        }   

        //diffuse:
        diffuseField(density, 0.0001); //second argument characterizes diffusion constant

        //advect:
        advectField(density, ux, uy);
        

        //update velocity and density for cells:
        for (int i = 0; i < meshSize; i++)
        {
            for (int j = 0; j < meshSize; j++)
            {
                float x = 2*(j+1)*meshWidth/wSize.x, y = (i+1)*meshHeight/wSize.y;
                ux[i+1][j+1] = Uscale*M_PI*(sin(M_PI*f(x,t))*cos(M_PI*y));
                uy[i+1][j+1] = -Uscale*M_PI*(cos(M_PI*f(x,t))*(fPrime(x,t))*sin(M_PI*y));
        
                int cellCoord = cellIndxMap(i, j);
                for(int cIndx1 = 0; cIndx1 < 2; cIndx1++)  
                {
                    for (int cIndx2 = 0; cIndx2 < 2; cIndx2++)
                    {
                        float tempRGB = cIndx1==0?255*density[i+1+cIndx1][j+1+cIndx2]:255*density[i+1+cIndx1][j+1+1-cIndx2];
                        Color tmpColor = Color(tempRGB, tempRGB, tempRGB);
                        cells[cellCoord + 2*cIndx1 + cIndx2].color = tmpColor;
                    }
                }
            }
        }
        t += dt;

        w.clear();
        w.draw(cells.data(), cells.size(), Quads);
        w.display();    
    
    
    }
    
    return 0;
}