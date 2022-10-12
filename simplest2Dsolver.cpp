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
const int meshSize = 250, meshSizeP2 = 252; // square mesh, size excludes the boundary cells; meshSize is no. of cells in each direction
float meshWidth, meshHeight;

//Velocity scale and time step:
float Uscale = 200.f;
float dt = 0.05;

int cellIndxMap(int i, int j)   //i for vertical, j for horizontal
{
    // return the 1D equivalent index of cell, where each cell has 4 points
    return 4*(j + meshSize*i);
}

float lerp(float a, float b, float k)
{
    return a+k*(b-a);
}

void diffuseField(float field[meshSizeP2][meshSizeP2], float delta)   //pass 2d array by reference
{
    float fieldNext[meshSizeP2][meshSizeP2] = {};
    float del = delta;    //del = 4*kappa*dt/h^2, it controls the diffusion strength
    int iterMax = 7;   // should compare with prev iteration step to terminate iteration loop
    for (int iter = 0; iter <= iterMax; iter++)
    {
       for (int i = 0; i < meshSize; i++)
       {
           for (int j = 0; j < meshSize; j++)
           {
                if(iter==0) fieldNext[i+1][j+1] = field[i+1][j+1];  //start with prev field, boundary cells are already always zero
                else if(iter>0 && iter<iterMax)    fieldNext[i+1][j+1] = (fieldNext[i+1][j+1] + del*0.25*(field[i+1+1][j+1]+field[i+1-1][j+1]+field[i+1][j+1+1]+field[i+1][j+1-1] ))/(1+del);
                else field[i+1][j+1] = fieldNext[i+1][j+1];
           }
       }
    }
    //memcpy(field, fieldNext, sizeof(field));
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
    //memcpy(field, tempField, sizeof(field));
}

void project(float pField[meshSizeP2][meshSizeP2], float uxField[][meshSizeP2], float uyField[][meshSizeP2])
{
    int iterMax = 20;   // should compare with prev iteration step to terminate iteration loop
    float h = meshWidth; //same as meshHeight
    float pFieldPrev[meshSizeP2][meshSizeP2] = {};
    float divField[meshSizeP2][meshSizeP2] = {};
    for (int iter = 0; iter <= iterMax; iter++)
    {
       for (int i = 0; i < meshSize; i++)
       {
            for (int j = 0; j < meshSize; j++)
            {
                if(iter==0)    divField[i+1][j+1] = h*0.25*(uxField[i+1][j+1+1] - uxField[i+1][j+1-1] + uyField[i+1+1][j+1] - uyField[i+1-1][j+1] )/2.f;    // pField includes delT factor
                pField[i+1][j+1] = -1.0*divField[i+1][j+1] + 0.25*(pField[i+1+1][j+1]+pField[i+1-1][j+1]+pField[i+1][j+1+1]+pField[i+1][j+1-1] );
            }
       }
       //continuity for p at boundaries:
       for (int k = 1; k <= meshSize; k++)
       {
            pField[0][k] = pField[1][k];
            pField[meshSize+1][k] = pField[meshSize][k];
            pField[k][0] = pField[k][1];
            pField[k][meshSize+1] = pField[k][meshSize];
       }
       pField[0][0] = 0.5*(pField[0][1]+pField[1][0]);
       pField[0][meshSize+1] = 0.5*(pField[0][meshSize]+pField[1][meshSize+1]);
       pField[meshSize+1][0] = 0.5*(pField[meshSize+1][1]+pField[meshSize][0]);
       pField[meshSize+1][meshSize+1] = 0.5*(pField[meshSize+1][meshSize]+pField[meshSize][meshSize+1]);
       
       //Check error:
    //    float sum = 0.f, tolerance = 0.05;
    //    for (int i = 0; i < meshSize+2; i++)
    //    {
    //         for (int j = 0; j < meshSize+2; j++)
    //         {
    //             sum += abs(pFieldPrev[i][j] - pField[i][j]);
    //             pFieldPrev[i][j] = pField[i][j];
    //         }
    //    }
    //    if(sum/pow(meshSizeP2, 2.0)<tolerance) break;//{cout<<"iter: "<<iter<<" "<<sum<<"; "; break;} 
       //if(iter==iterMax && sum/pow(meshSizeP2, 2.0)>tolerance)  cout<<sum/pow(meshSizeP2, 2.0)<<endl;
    }
    for (int i = 0; i < meshSize; i++)
    {
        for (int j = 0; j < meshSize; j++)
        {
            uxField[i+1][j+1] -= (pField[i+1][j+1+1]-pField[i+1][j+1-1])/(2.f*h);
            uyField[i+1][j+1] -= (pField[i+1+1][j+1]-pField[i+1-1][j+1])/(2.f*h);
        }
    }

}


int main()
{
    //define fields:
    float density[meshSizeP2][meshSizeP2] = {}; // 2d array of zeroes.
    float ux[meshSizeP2][meshSizeP2] = {}; // 2d array of zeroes.
    float uy[meshSizeP2][meshSizeP2] = {}; // 2d array of zeroes.
    float uxPrev[meshSizeP2][meshSizeP2] = {}; // 2d array of zeroes.
    float uyPrev[meshSizeP2][meshSizeP2] = {}; // 2d array of zeroes.
    
    float pField[meshSizeP2][meshSizeP2] = {};

    meshWidth = (float)wSize.x/(float)meshSize;
    meshHeight = (float)wSize.y/(float)meshSize;

    dt = 2.f*meshWidth/Uscale;
    cout<<"dt: "<<dt<<endl;

    RenderWindow w(VideoMode(wSize.x, wSize.y), "simplest2Dsolver");
    w.setFramerateLimit(60);

    //**********Create Mesh***********
    //define vertex array (quads) for each cell:
    vector<Vertex> cells; //4 vertices needed to define quad, need to be in a particular order (clockwise or anti-clockwise)  
    for (int i = 0; i < meshSize; i++)  //cellsPrev stores previous color info
    {
        for (int j = 0; j < meshSize; j++)
        {
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
    cout<<"meshWidth: "<<meshWidth<<endl;
    //Use mouse click to add density
    Vector2i localMousePos(0, 0), localMousePrevPos(0, 0);
    Color currColor = Color::White;
    
    bool begin = false, mouseDown = false;
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
                    for (int i = 0; i < meshSize; i++)
                    {
                        for (int j = 0; j < meshSize; j++)
                        {
                            density[i+1][j+1] = 0.f;
                            ux[i+1][j+1] = 0.f;
                            uy[i+1][j+1] = 0.f;
                            uxPrev[i+1][j+1] = 0.f;
                            uyPrev[i+1][j+1] = 0.f;
                        }
                    }
                }
                if(e.key.code==Keyboard::B) begin=begin?false:true;

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

        //add source:
        localMousePos = Mouse::getPosition(w);
        int sourceSize = 5; //radius of source 
        if(mouseDown)
        {
            int localJ = floor(localMousePos.x/meshWidth), localI = floor(localMousePos.y/meshHeight);
            for (int iSource = localI-sourceSize; iSource <= localI+sourceSize; iSource++)
            {
                if(iSource<0 || iSource>=meshSize)   continue;
                for (int jSource = localJ-sourceSize; jSource <= localJ+sourceSize; jSource++)
                {
                    if(jSource<0 || jSource>=meshSize)   continue;
                    
                    uy[iSource+1][jSource+1] = -Uscale;
                    // if(localMousePos.x != localMousePrevPos.x && localMousePos.y != localMousePrevPos.y)
                    // {
                    //     ux[iSource+1][jSource+1] = Uscale*tanh(localMousePos.x-localMousePrevPos.x);
                    //     uy[iSource+1][jSource+1] = Uscale*tanh(localMousePos.y-localMousePrevPos.y);
                    // }
                    // else    uy[iSource+1][jSource+1] = -Uscale;
                    float localD = abs(jSource - localJ)*meshWidth + abs(iSource - localI)*meshHeight, dMax = sourceSize*(meshWidth+meshHeight);
                    
                    density[iSource + 1][jSource + 1] = density[iSource + 1][jSource + 1]<1.f-localD/(2.f*dMax)?1.f-localD/(2.f*dMax):density[iSource + 1][jSource + 1]; // added ones for the boundary offset
                    int cellCoord = cellIndxMap(iSource, jSource);
                    for(int cIndx = 0; cIndx < 4; cIndx++)  cells[cellCoord + cIndx].color = currColor;
                }
                
            }
        }   

        if(begin)
        {
            //diffuse:
            diffuseField(density, 0.00001); //second argument characterizes diffusion constant

            //advect:
            advectField(density, ux, uy);

            //*******Calculate velocity for next instant:**********
           // diffuseField(ux, 0.001); // second argument specifies viscosity
           // diffuseField(uy, 0.001);

           // project(pField, ux, uy);    // remove the divergence part
            
            for (int i = 0; i < meshSize; i++)
            {
                for (int j = 0; j < meshSize; j++)
                {
                    uxPrev[i+1][j+1] = ux[i+1][j+1];
                    uyPrev[i+1][j+1] = uy[i+1][j+1];
                }
            }

            advectField(uxPrev, ux, uy);
            advectField(uyPrev, ux, uy);

            for (int i = 0; i < meshSize; i++)
            {
                for (int j = 0; j < meshSize; j++)
                {
                    ux[i+1][j+1] = uxPrev[i+1][j+1];
                    uy[i+1][j+1] = uyPrev[i+1][j+1];
                    density[i+1][j+1] -= dissipation*dt*density[i+1][j+1]; //dissipate density
                }
            }

           // diffuseField(ux, 0.01); // second argument specifies viscosity
           // diffuseField(uy, 0.01);

            project(pField, ux, uy);
            // //*****************
            //update density for cells:
            for (int i = 0; i < meshSize; i++)
            {
                for (int j = 0; j < meshSize; j++)
                {
                    int cellCoord = cellIndxMap(i, j);
                    for(int cIndx1 = 0; cIndx1 < 2; cIndx1++)  
                    {
                        for (int cIndx2 = 0; cIndx2 < 2; cIndx2++)
                        {
                            float tempRGB = cIndx1==0?255.f*density[i+1+cIndx1][j+1+cIndx2]:255.f*density[i+1+cIndx1][j+1+1-cIndx2];
                            Color tmpColor = currColor*Color(tempRGB, tempRGB, tempRGB);
                            cells[cellCoord + 2*cIndx1 + cIndx2].color = tmpColor;
                        }
                    }
                }
            }
        
        }

        
        w.clear();
        w.draw(cells.data(), cells.size(), Quads);
        w.display();    
    
    
    }
    
    return 0;
}